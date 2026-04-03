# ==============================================================================
# run_cal.R — Monte Carlo calibration using Morris-informed parameter bounds
#
# USAGE (RStudio / VSCode):
#   CONFIG_PATH <- "config.yaml"
#   source("run_cal.R")
#
# USAGE (terminal):
#   Rscript run_cal.R config.yaml
#
# REQUIRES: a GSA_v* folder must exist in out_path (produced by run_gsa.R)
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
library(tidyr)
library(purrr)
library(hydroGOF)
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
cat(" Calibration — Monte Carlo sampling\n")
cat(" Config:", config_path, "\n")
cat("================================================================\n")

# ------------------------------------------------------------------------------
# 2. MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/morris_classify.R")
source("R/plots.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. OUTPUT PATH + VERSIONED FOLDER
# ------------------------------------------------------------------------------
out_path <- if (!is.null(cfg$output_path)) cfg$output_path else cfg$project_path
ver      <- create_versioned_folder(out_path, "CAL_v")
cat("Output folder:", ver$path, "\n\n")

# ------------------------------------------------------------------------------
# 4. LOAD GSA RESULTS (param_info + Morris ranking)
#    Always reads the latest GSA_v* folder so GSA and CAL are paired by version.
# ------------------------------------------------------------------------------
gsa_folder     <- find_latest_folder(out_path, "GSA_v")
param_info     <- read.csv(file.path(gsa_folder$path, "param_info_gsa.csv"),
                           stringsAsFactors = FALSE) %>% as_tibble()
ranking_morris <- read.csv(file.path(gsa_folder$path, "ranking_global_gsa.csv"),
                           stringsAsFactors = FALSE) %>% as_tibble()

cat("GSA source:", gsa_folder$path, "\n")

# ------------------------------------------------------------------------------
# 4b. EXCLUDE INSENSITIVE PARAMETERS
#     Two modes (can be combined):
#       1. auto_exclude: true → exclude parameters with mu_star <= threshold
#       2. exclude_parameters: [...] → manual list of parameter short names
#     Excluded parameters are fixed at their GSA median value.
# ------------------------------------------------------------------------------

# --- Auto-exclude based on Morris mu_star ---
auto_exclude <- isTRUE(cfg$calibration$auto_exclude)
auto_threshold <- if (!is.null(cfg$calibration$auto_exclude_threshold)) {
  cfg$calibration$auto_exclude_threshold
} else {
  0
}

exclude_names <- character(0)

if (auto_exclude) {
  insensitive <- ranking_morris %>%
    filter(metric == "obj_total" & mu_star <= auto_threshold) %>%
    pull(parameter)
  exclude_names <- union(exclude_names, insensitive)
  cat("Auto-exclude (mu_star <=", auto_threshold, "):",
      if (length(insensitive) > 0) paste(insensitive, collapse = ", ") else "none", "\n")
}

# --- Manual exclude list (additive) ---
manual_list <- cfg$calibration$exclude_parameters
if (!is.null(manual_list) && length(manual_list) > 0) {
  manual_names <- unlist(manual_list)
  exclude_names <- union(exclude_names, manual_names)
  cat("Manual exclude:", paste(manual_names, collapse = ", "), "\n")
}

# --- Apply exclusion ---
if (length(exclude_names) > 0) {
  # Load GSA metrics to compute median values for excluded parameters
  gsa_rds <- readRDS(file.path(gsa_folder$path, "resultado_gsa.rds"))
  gsa_scaled <- gsa_rds$design_norm
  gsa_pi     <- gsa_rds$param_info

  # Scale design to physical values for median computation
  for (k in seq_len(nrow(gsa_pi))) {
    nm <- gsa_pi$short[k]
    gsa_scaled[[nm]] <- gsa_pi$min[k] + gsa_scaled[[nm]] * (gsa_pi$max[k] - gsa_pi$min[k])
  }

  excluded_params <- param_info %>% filter(short %in% exclude_names)
  param_info      <- param_info %>% filter(!short %in% exclude_names)

  cat("\nExcluded from calibration (fixed at GSA median):\n")
  for (i in seq_len(nrow(excluded_params))) {
    med_val <- median(gsa_scaled[[excluded_params$short[i]]], na.rm = TRUE)
    excluded_params$fixed_value[i] <- med_val
    cat("  ", excluded_params$short[i], "=", round(med_val, 4), "\n")
  }
  cat("\n")
} else {
  excluded_params <- NULL
  cat("No parameters excluded.\n\n")
}

cat("Calibration parameters:", nrow(param_info), "\n")
print(param_info[, c("short", "min", "max")])
cat("\n")

# ------------------------------------------------------------------------------
# 5. SAMPLE PARAMETERS (uniform random within GSA bounds)
# ------------------------------------------------------------------------------
n_sim      <- cfg$calibration$n_simulations
cal_params <- map_dfc(seq_len(nrow(param_info)), function(i) {
  runif(n_sim, min = param_info$min[i], max = param_info$max[i])
})
names(cal_params) <- param_info$par_name
cal_params        <- as_tibble(cal_params)

# Add excluded parameters as fixed columns
if (!is.null(excluded_params) && nrow(excluded_params) > 0) {
  for (i in seq_len(nrow(excluded_params))) {
    cal_params[[excluded_params$par_name[i]]] <- excluded_params$fixed_value[i]
  }
}

cat("Simulations:", n_sim, "\n")

# ------------------------------------------------------------------------------
# 6. DEFINE SWAT+ OUTPUTS
# ------------------------------------------------------------------------------
outputs <- lapply(cfg$cal_outputs, function(o) {
  define_output(file = o$file, variable = o$variable, unit = o$unit)
})

# ------------------------------------------------------------------------------
# 7. CLEAN SWAT+ CACHE
# ------------------------------------------------------------------------------
cache_dirs <- list.files(cfg$project_path, pattern = "^(run_|cal_)",
                         full.names = TRUE, include.dirs = TRUE)
cache_dirs <- cache_dirs[file.info(cache_dirs)$isdir]
if (length(cache_dirs) > 0) {
  unlink(cache_dirs, recursive = TRUE, force = TRUE)
  cat("Cache cleared (", length(cache_dirs), " folder(s)).\n", sep = "")
}

# ------------------------------------------------------------------------------
# 8. RUN SWAT+
# ------------------------------------------------------------------------------
cat("Running SWAT+ (", cfg$threads, " threads)...\n", sep = "")

sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = cal_params,
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
cat("Valid runs:", length(valid_runs), "/", n_sim, "\n")

# ------------------------------------------------------------------------------
# 10. OBSERVED DATA + COMPARISON TABLE
# ------------------------------------------------------------------------------
out_nm    <- names(cfg$cal_outputs)[1]
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

comp       <- build_comparativo(sim$simulation[[out_nm]], obs_data, valid_runs, floor_mon)
filt       <- filter_invalid_runs(comp, valid_runs)
valid_runs <- filt$runs_ok
comp       <- comp %>% select(Date, Flow, all_of(valid_runs))

# ------------------------------------------------------------------------------
# 11. METRICS + CLASSIFICATION
# ------------------------------------------------------------------------------
metrics <- calc_hard_metrics(comp, valid_runs, suffix = suffix)
metrics <- classify_hard_runs(
  metrics,
  cfg$calibration$threshold_nse,
  cfg$calibration$threshold_kge,
  cfg$calibration$threshold_abs_pbias,
  suffix = suffix
)

# Join parameter values to results
params_idx <- cal_params %>% mutate(run_index = row_number())
results    <- metrics %>% left_join(params_idx, by = "run_index")
sat        <- results %>% filter(satisfactory_run)

n_sat <- nrow(sat)
cat("\nSatisfactory runs:", n_sat, "/", length(valid_runs),
    "(", round(100 * n_sat / max(length(valid_runs), 1), 1), "%)\n\n")

# ------------------------------------------------------------------------------
# 11b. WATER BALANCE DIAGNOSTIC (optional — if targets are defined)
# ------------------------------------------------------------------------------
balance_diag <- NULL
targets <- cfg$calibration$targets
has_balance_outputs <- all(c("precip", "et", "wateryld") %in% names(sim$simulation))

if (!is.null(targets) && has_balance_outputs) {
  balance_diag <- calc_balance_diagnostic(sim, valid_runs, targets)

  # Attach to results
  results <- results %>%
    left_join(balance_diag %>% select(run_id, ET_P, WYLD_P, err_ET_P_pct,
                                       err_WYLD_P_pct, agg_rel_error_pct),
              by = "run_id")
  sat <- results %>% filter(satisfactory_run)

  # Summary for satisfactory runs
  if (n_sat > 0) {
    sat_balance <- balance_diag %>% filter(run_id %in% sat$run_id)
    cat("=== WATER BALANCE DIAGNOSTIC (satisfactory runs) ===\n")
    cat(sprintf("  Target ET/P  = %.4f  |  Sim median = %.4f  (range: %.4f - %.4f)\n",
                targets$et_rto,
                median(sat_balance$ET_P), min(sat_balance$ET_P), max(sat_balance$ET_P)))
    cat(sprintf("  Target WYLD/P = %.4f  |  Sim median = %.4f  (range: %.4f - %.4f)\n",
                targets$wyld_rto,
                median(sat_balance$WYLD_P), min(sat_balance$WYLD_P), max(sat_balance$WYLD_P)))
    cat(sprintf("  Median aggregate error = %.1f%%\n\n",
                median(sat_balance$agg_rel_error_pct)))
  }
} else if (!is.null(targets) && !has_balance_outputs) {
  cat("NOTE: Water balance targets defined but precip/et/wateryld outputs not found.\n")
  cat("  Add precip, et, wateryld to cal_outputs in config.yaml.\n\n")
}

# ------------------------------------------------------------------------------
# 12. PLOTS
# ------------------------------------------------------------------------------
nse_col <- paste0("NSE", suffix)
kge_col <- paste0("KGE", suffix)

g_scatter <- plot_performance_scatter(
  metrics, nse_col, kge_col,
  cfg$calibration$threshold_nse,
  cfg$calibration$threshold_kge,
  paste(cfg$temporal_scale, "calibration — performance space")
)
print(g_scatter)
save_tiff(g_scatter, "scatter_performance.tif", ver$path)

# Behavioural band + P-factor / R-factor
band_result <- NULL
if (n_sat >= 2) {
  band_result <- calc_behavioural_band(comp, sat$run_id)
  cat(sprintf("Behavioural band  P-factor = %.3f  |  R-factor = %.3f\n",
              band_result$p_factor, band_result$r_factor))
}

g_hydro <- plot_hydrograph(
  comp,
  satisfactory_runs = if (n_sat > 0) sat$run_id else NULL,
  band              = if (!is.null(band_result)) band_result$band     else NULL,
  p_factor          = if (!is.null(band_result)) band_result$p_factor else NULL,
  r_factor          = if (!is.null(band_result)) band_result$r_factor else NULL,
  title_label       = paste(cfg$temporal_scale, "streamflow calibration")
)
print(g_hydro)
save_tiff(g_hydro, "hydrograph.tif", ver$path, width = 20, height = 10)

# Uncertainty envelope plot (SWAT-CUP style: band + observed only)
if (!is.null(band_result)) {
  obs_trim <- obs_data %>% filter(Date %in% comp$Date)
  g_envelope <- plot_uncertainty_envelope(
    band        = band_result$band,
    obs_data    = obs_trim,
    p_factor    = band_result$p_factor,
    r_factor    = band_result$r_factor,
    title_label = paste(cfg$temporal_scale, "— behavioural uncertainty band")
  )
  print(g_envelope)
  save_tiff(g_envelope, "uncertainty_envelope.tif", ver$path, width = 20, height = 10)
}

# FDC envelope (Flow Duration Curve with behavioural band)
if (!is.null(band_result)) {
  fdc_data <- calc_fdc_band(comp, sat$run_id)
  if (!is.null(fdc_data)) {
    g_fdc <- plot_fdc_envelope(
      fdc_data,
      title_label = paste(cfg$temporal_scale, "— Flow Duration Curve")
    )
    print(g_fdc)
    save_tiff(g_fdc, "fdc_envelope.tif", ver$path, width = 20, height = 12)
    write.csv(fdc_data, file.path(ver$path, "fdc_band.csv"), row.names = FALSE)
  }
}

if (n_sat > 0) {
  g_box <- plot_param_boxplot(
    sat, param_info,
    "Satisfactory parameter distributions"
  )
  print(g_box)
  save_tiff(g_box, "boxplot_params.tif", ver$path, width = 25, height = 10)
}

# ------------------------------------------------------------------------------
# 13. NEW PARAMETER RANGES (for next iteration)
# ------------------------------------------------------------------------------
if (n_sat > 0) {
  new_ranges <- calc_new_ranges(sat, param_info$par_name, param_info, ranking_morris)

  cat("=== NEW PARAMETER RANGES ===\n")
  print(new_ranges)

  ranges_file <- file.path(ver$path, "new_ranges.csv")
  write.csv(new_ranges, ranges_file, row.names = FALSE)

  cat("\nNew ranges saved to:\n ", ranges_file, "\n\n")
  cat("To run another iteration, update config.yaml:\n")
  cat("  iteration:\n")
  cat("    enabled: true\n")
  cat("    ranges_file: \"", ranges_file, "\"\n", sep = "")
} else {
  cat("WARNING: No satisfactory runs found.\n")
  cat("  Consider relaxing calibration thresholds or increasing n_simulations.\n")
}

# ------------------------------------------------------------------------------
# 14. SAVE RESULTS
# ------------------------------------------------------------------------------
saveRDS(
  list(
    config               = cfg,
    gsa_source           = gsa_folder$path,
    sim                  = sim,
    param_info           = param_info,
    excluded_params      = excluded_params,
    results              = results,
    satisfactory_results = sat,
    behavioural_band     = band_result,
    balance_diagnostic   = balance_diag
  ),
  file.path(ver$path, "resultado_cal.rds")
)

write.csv(results, file.path(ver$path, "results_cal.csv"), row.names = FALSE)

if (n_sat > 0) {
  write.csv(sat, file.path(ver$path, "satisfactory_results.csv"), row.names = FALSE)
}

if (!is.null(band_result)) {
  write.csv(band_result$band, file.path(ver$path, "behavioural_band.csv"), row.names = FALSE)
  write.csv(
    data.frame(p_factor = band_result$p_factor, r_factor = band_result$r_factor),
    file.path(ver$path, "behavioural_band_factors.csv"), row.names = FALSE
  )
}

if (!is.null(balance_diag)) {
  write.csv(balance_diag, file.path(ver$path, "balance_diagnostic.csv"), row.names = FALSE)
}

cat("\nCalibration completed. Results saved to:\n ", ver$path, "\n")
