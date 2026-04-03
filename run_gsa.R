# ==============================================================================
# run_gsa.R — Global Sensitivity Analysis (LH-OAT Morris method)
#
# USAGE (RStudio / VSCode):
#   CONFIG_PATH <- "config.yaml"
#   source("run_gsa.R")
#
# USAGE (terminal):
#   Rscript run_gsa.R config.yaml
# ==============================================================================

# Preserve variables set by main.R before clearing the environment
.bak_cfg   <- if (exists("CONFIG_PATH"))         CONFIG_PATH         else NULL
.bak_state <- if (exists(".CALIBRATOR_STATE")) .CALIBRATOR_STATE else NULL
rm(list = ls())
gc()
if (!is.null(.bak_cfg))   CONFIG_PATH         <- .bak_cfg
if (!is.null(.bak_state)) .CALIBRATOR_STATE   <- .bak_state
rm(.bak_cfg, .bak_state)

library(SWATrunR)
library(yaml)
library(dplyr)
library(tibble)
library(lhs)
library(purrr)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. CONFIG
# ------------------------------------------------------------------------------
config_path <- if (exists("CONFIG_PATH")) {
  CONFIG_PATH
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) args[1] else "config.yaml"
}

cfg <- yaml::read_yaml(config_path)
cat("================================================================\n")
cat(" GSA — LH-OAT Morris screening\n")
cat(" Config:", config_path, "\n")
cat("================================================================\n")

# ------------------------------------------------------------------------------
# 2. MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/lhoat_engine.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/morris_classify.R")
source("R/plots.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. OUTPUT PATH + VERSIONED FOLDER
# ------------------------------------------------------------------------------
out_path <- if (!is.null(cfg$output_path)) cfg$output_path else cfg$project_path
ver      <- create_versioned_folder(out_path, "GSA_v")
cat("Output folder:", ver$path, "\n\n")

# ------------------------------------------------------------------------------
# 4. PARAMETER INFO
#    Default: use bounds from config.yaml
#    Iteration mode: override bounds from new_ranges.csv of a previous run
# ------------------------------------------------------------------------------
yaml_params <- tibble(
  short    = sapply(cfg$parameters, `[[`, "short"),
  par_name = sapply(cfg$parameters, `[[`, "par_name"),
  min      = as.numeric(sapply(cfg$parameters, `[[`, "min")),
  max      = as.numeric(sapply(cfg$parameters, `[[`, "max"))
)

if (isTRUE(cfg$iteration$enabled) && !is.null(cfg$iteration$ranges_file)) {
  cat("Iteration mode: loading ranges from\n ", cfg$iteration$ranges_file, "\n")
  prev_ranges <- read.csv(cfg$iteration$ranges_file, stringsAsFactors = FALSE)

  param_info <- yaml_params %>%
    inner_join(prev_ranges %>% select(parameter, new_min, new_max),
               by = c("short" = "parameter")) %>%
    mutate(
      min = new_min,
      max = new_max
    ) %>%
    select(short, par_name, min, max)

  cat("Ranges updated from previous calibration.\n")
} else {
  param_info <- yaml_params
  cat("Using parameter bounds from config.yaml.\n")
}

p <- nrow(param_info)
cat("\nParameters (", p, "):\n", sep = "")
print(param_info)
cat("\n")

# ------------------------------------------------------------------------------
# 5. GENERATE LH-OAT TRAJECTORIES
# ------------------------------------------------------------------------------
design_norm <- generate_lhoat_trajectories(
  m_traj      = cfg$gsa$m_traj,
  p           = p,
  delta       = cfg$gsa$delta,
  param_names = param_info$short
)

scaled <- scale_parameters(design_norm, param_info)
cat("Trajectories:", cfg$gsa$m_traj, " | Total runs:", nrow(scaled$swat_params), "\n\n")

# ------------------------------------------------------------------------------
# 6. CLEAN SWAT+ CACHE
# ------------------------------------------------------------------------------
cache_dirs <- list.files(cfg$project_path, pattern = "^(run_|lhoat_|gsa_)",
                         full.names = TRUE, include.dirs = TRUE)
cache_dirs <- cache_dirs[file.info(cache_dirs)$isdir]
if (length(cache_dirs) > 0) {
  unlink(cache_dirs, recursive = TRUE, force = TRUE)
  cat("Cache cleared (", length(cache_dirs), " folder(s)).\n", sep = "")
}

# ------------------------------------------------------------------------------
# 7. DEFINE SWAT+ OUTPUTS
# ------------------------------------------------------------------------------
outputs <- lapply(cfg$gsa_outputs, function(o) {
  define_output(file = o$file, variable = o$variable, unit = o$unit)
})

# ------------------------------------------------------------------------------
# 8. RUN SWAT+
# ------------------------------------------------------------------------------
cat("Running SWAT+ (", cfg$threads, " threads)...\n", sep = "")

sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = scaled$swat_params,
  start_date   = cfg$simulation$start_date,
  end_date     = cfg$simulation$end_date,
  years_skip   = cfg$simulation$warmup_years,
  n_thread     = cfg$threads,
  save_file    = ver$tag,
  keep_folder  = TRUE
)

# ------------------------------------------------------------------------------
# 9. VALID RUNS
# ------------------------------------------------------------------------------
valid_runs <- get_valid_runs(sim$simulation, names(outputs))
cat("Valid runs:", length(valid_runs), "/", nrow(scaled$swat_params), "\n")

# ------------------------------------------------------------------------------
# 10. METRICS
# ------------------------------------------------------------------------------
out_nm    <- names(cfg$gsa_outputs)[1]
suffix    <- if (cfg$temporal_scale == "daily") "_day" else "_mon"
floor_mon <- cfg$temporal_scale == "monthly"

obs_file <- if (cfg$temporal_scale == "daily") {
  cfg$observed$daily_file
} else {
  cfg$observed$monthly_file
}

obs_data <- read.csv(file.path(cfg$project_path, obs_file), stringsAsFactors = FALSE) %>%
  mutate(Date = as.Date(Date))

if (floor_mon) obs_data <- obs_data %>% mutate(Date = floor_date(Date, "month"))

obs_data <- obs_data %>%
  filter(Date >= as.Date(cfg$observed$start_date) &
           Date <= as.Date(cfg$simulation$end_date)) %>%
  arrange(Date)

comparativo <- build_comparativo(sim$simulation[[out_nm]], obs_data, valid_runs, floor_mon)
filt        <- filter_invalid_runs(comparativo, valid_runs)
valid_runs  <- filt$runs_ok
comparativo <- comparativo %>% select(Date, Flow, all_of(valid_runs))

metrics <- calc_hard_metrics(comparativo, valid_runs, suffix = suffix) %>%
  mutate(run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_id)))

# ------------------------------------------------------------------------------
# 11. ATTACH TRAJECTORY DATA + COMPUTE ELEMENTARY EFFECTS
# ------------------------------------------------------------------------------
metrics <- design_norm %>%
  select(run_index, traj_id, step_id, all_of(param_info$short)) %>%
  inner_join(metrics, by = "run_index") %>%
  arrange(traj_id, step_id)

ee_all <- calc_all_ee(metrics, cfg$gsa_metrics, param_info$short)
sens   <- compute_sensitivity(ee_all)

cat("\n=== PARAMETER RANKING (obj_total) ===\n")
print(sens$ranking)

# ------------------------------------------------------------------------------
# 12. MORRIS CLASSIFICATION + PLOT
# ------------------------------------------------------------------------------
ranking_classified <- classify_morris(sens$ranking)

iter_label <- if (isTRUE(cfg$iteration$enabled)) "Iteration 2+" else "Iteration 1"
g_morris   <- plot_morris(ranking_classified, iter_label)
print(g_morris)
save_tiff(g_morris, "Morris_sensitivity_screening.tif", ver$path)

# ------------------------------------------------------------------------------
# 13. SAVE RESULTS
# ------------------------------------------------------------------------------
saveRDS(
  list(
    config      = cfg,
    sim         = sim,
    ranking     = sens$ranking,
    sensitivity = sens$sensitivity,
    metrics     = metrics,
    ee_all      = ee_all,
    param_info  = param_info,
    design_norm = design_norm,
    swat_params = scaled$swat_params
  ),
  file.path(ver$path, "resultado_gsa.rds")
)

write.csv(sens$ranking,  file.path(ver$path, "ranking_global_gsa.csv"), row.names = FALSE)
write.csv(metrics,       file.path(ver$path, "metricas_gsa.csv"),       row.names = FALSE)
write.csv(ee_all,        file.path(ver$path, "ee_all_gsa.csv"),         row.names = FALSE)
write.csv(param_info,    file.path(ver$path, "param_info_gsa.csv"),     row.names = FALSE)

cat("\nGSA completed. Results saved to:\n ", ver$path, "\n")
