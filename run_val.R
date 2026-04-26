# ==============================================================================
# run_val.R — Split-sample validation of a calibrated ensemble
#
# Loads a set of parameter combinations produced by a previous calibration
# (typically the satisfactory or robust subset) and propagates them through
# SWAT+ over an independent period. Computes NSE/KGE/PBIAS and the behavioural
# band (P-factor / R-factor) on the validation window.
#
# USAGE (RStudio / VSCode):
#   CONFIG_PATH <- "config.yaml"
#   source("run_val.R")
#
# USAGE (terminal):
#   Rscript run_val.R config.yaml
#
# REQUIRES: validation block in config.yaml (see config.yaml for schema).
# ==============================================================================

# Preserve variables set by main.R before clearing the environment
.bak_cfg   <- if (exists("CONFIG_PATH"))       CONFIG_PATH       else NULL
.bak_state <- if (exists(".CALIBRATOR_STATE")) .CALIBRATOR_STATE else NULL
rm(list = ls())
gc()
if (!is.null(.bak_cfg))   CONFIG_PATH       <- .bak_cfg
if (!is.null(.bak_state)) .CALIBRATOR_STATE <- .bak_state
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
cat(" Validation — split-sample test on independent period\n")
cat(" Config:", config_path, "\n")
cat("================================================================\n")

if (is.null(cfg$validation)) {
  stop("No 'validation' block found in ", config_path)
}
if (is.null(cfg$validation$param_sets_csv)) {
  stop("'validation.param_sets_csv' must point to a CSV with calibrated parameter sets.")
}

# Resolve validation-specific overrides with fall-back to top-level blocks
val_cfg     <- cfg$validation
sim_cfg     <- if (!is.null(val_cfg$simulation))     val_cfg$simulation     else cfg$simulation
obs_cfg     <- if (!is.null(val_cfg$observed))       val_cfg$observed       else cfg$observed
temporal_sc <- if (!is.null(val_cfg$temporal_scale)) val_cfg$temporal_scale else cfg$temporal_scale
val_outputs_cfg <- if (!is.null(val_cfg$outputs)) val_cfg$outputs else cfg$cal_outputs

# Thresholds used to classify "still satisfactory in validation" (temporal stability)
val_thr <- val_cfg$thresholds
if (is.null(val_thr)) val_thr <- cfg$calibration[c("threshold_nse",
                                                   "threshold_kge",
                                                   "threshold_abs_pbias")]

# ------------------------------------------------------------------------------
# 2. MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/plots.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. OUTPUT PATH + VERSIONED FOLDER
# ------------------------------------------------------------------------------
out_path <- if (!is.null(cfg$output_path)) cfg$output_path else cfg$project_path
ver      <- create_versioned_folder(out_path, "VAL_v")
cat("Output folder:", ver$path, "\n\n")

# ------------------------------------------------------------------------------
# 4. LOAD PARAMETER SETS
#    Accept any CSV that contains at least one column named with SWATrunR
#    notation (e.g. "cn2::cn2.hru | change = pctchg").
# ------------------------------------------------------------------------------
param_sets_path <- val_cfg$param_sets_csv
if (!file.exists(param_sets_path)) {
  stop("Parameter-sets CSV not found: ", param_sets_path)
}

param_sets <- read.csv(param_sets_path, stringsAsFactors = FALSE,
                       check.names = FALSE) %>% as_tibble()

par_cols <- grep("::", names(param_sets), value = TRUE, fixed = TRUE)
if (length(par_cols) == 0) {
  stop("No SWATrunR-style parameter columns (containing '::') in: ", param_sets_path)
}

val_params <- param_sets[, par_cols, drop = FALSE]

# Preserve original run ids for traceability between calibration and validation
orig_run_ids <- if ("run_id" %in% names(param_sets)) {
  as.character(param_sets$run_id)
} else {
  sprintf("run_%03d", seq_len(nrow(param_sets)))
}

n_sim <- nrow(val_params)
cat("Parameter sets file :", param_sets_path, "\n")
cat("Ensemble size       :", n_sim, "parameter sets\n")
cat("Parameters per set  :", length(par_cols), "\n\n")

# ------------------------------------------------------------------------------
# 5. DEFINE SWAT+ OUTPUTS
# ------------------------------------------------------------------------------
outputs <- lapply(val_outputs_cfg, function(o) {
  define_output(file = o$file, variable = o$variable, unit = o$unit)
})

# ------------------------------------------------------------------------------
# 6. CLEAN SWAT+ CACHE
# ------------------------------------------------------------------------------
cache_dirs <- list.files(cfg$project_path, pattern = "^(run_|val_)",
                         full.names = TRUE, include.dirs = TRUE)
cache_dirs <- cache_dirs[file.info(cache_dirs)$isdir]
if (length(cache_dirs) > 0) {
  unlink(cache_dirs, recursive = TRUE, force = TRUE)
  cat("Cache cleared (", length(cache_dirs), " folder(s)).\n", sep = "")
}

# ------------------------------------------------------------------------------
# 7. RUN SWAT+ OVER THE VALIDATION PERIOD
# ------------------------------------------------------------------------------
cat("Running SWAT+ (", cfg$threads, " threads) on ",
    sim_cfg$start_date, " to ", sim_cfg$end_date,
    " | warm-up: ", sim_cfg$warmup_years, " yr\n", sep = "")

sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = val_params,
  start_date   = sim_cfg$start_date,
  end_date     = sim_cfg$end_date,
  years_skip   = sim_cfg$warmup_years,
  n_thread     = cfg$threads,
  save_file    = ver$tag,
  keep_folder  = TRUE
)

# ------------------------------------------------------------------------------
# 8. VALID RUNS
# ------------------------------------------------------------------------------
valid_runs <- get_valid_runs(sim$simulation, names(outputs))
cat("Valid runs:", length(valid_runs), "/", n_sim, "\n")

# Map SWATrunR run names (run_1, run_2, ...) to original calibration run ids
run_map <- tibble(
  run_id       = valid_runs,
  orig_run_id  = orig_run_ids[as.integer(sub("^run_0*", "", valid_runs))]
)

# ------------------------------------------------------------------------------
# 9. OBSERVED DATA + COMPARISON TABLE (validation window only)
# ------------------------------------------------------------------------------
out_nm    <- names(val_outputs_cfg)[1]
suffix    <- if (temporal_sc == "daily") "_day" else "_mon"
floor_mon <- temporal_sc == "monthly"

obs_file <- if (temporal_sc == "daily") {
  if (!is.null(obs_cfg$daily_file))   obs_cfg$daily_file   else cfg$observed$daily_file
} else {
  if (!is.null(obs_cfg$monthly_file)) obs_cfg$monthly_file else cfg$observed$monthly_file
}
obs_path <- file.path(cfg$project_path, obs_file)

obs_data <- read.csv(obs_path, stringsAsFactors = FALSE) %>%
  mutate(Date = as.Date(Date))

if (floor_mon) obs_data <- obs_data %>% mutate(Date = floor_date(Date, "month"))

val_obs_start <- as.Date(if (!is.null(obs_cfg$start_date)) obs_cfg$start_date else cfg$observed$start_date)
val_sim_end   <- as.Date(sim_cfg$end_date)

obs_data <- obs_data %>%
  filter(Date >= val_obs_start & Date <= val_sim_end) %>%
  arrange(Date)

comp       <- build_comparativo(sim$simulation[[out_nm]], obs_data, valid_runs, floor_mon)
filt       <- filter_invalid_runs(comp, valid_runs)
valid_runs <- filt$runs_ok
comp       <- comp %>% select(Date, Flow, all_of(valid_runs))

cat("Validation window   :", format(min(comp$Date)), "->", format(max(comp$Date)),
    " (", nrow(comp), " steps)\n", sep = "")

# ------------------------------------------------------------------------------
# 10. METRICS + TEMPORAL-STABILITY CLASSIFICATION
# ------------------------------------------------------------------------------
metrics <- calc_hard_metrics(comp, valid_runs, suffix = suffix)
metrics <- classify_hard_runs(
  metrics,
  val_thr$threshold_nse,
  val_thr$threshold_kge,
  val_thr$threshold_abs_pbias,
  suffix = suffix
) %>%
  rename(still_satisfactory = satisfactory_run)

# Attach parameter values and original run ids
val_params_idx <- val_params %>% mutate(run_id = sprintf("run_%d", row_number()))
results <- metrics %>%
  left_join(run_map,         by = "run_id") %>%
  left_join(val_params_idx,  by = "run_id")

n_stable <- sum(results$still_satisfactory, na.rm = TRUE)
cat("\n=== TEMPORAL STABILITY ===\n")
cat("Runs still satisfactory in validation:", n_stable, "/", length(valid_runs),
    sprintf(" (%.1f%%)\n", 100 * n_stable / max(length(valid_runs), 1)))

# ------------------------------------------------------------------------------
# 11. BEHAVIOURAL BAND — FULL ENSEMBLE (propagated calibration uncertainty)
#     p-/r-factor are computed over ALL valid runs; the calibrated ensemble is
#     not re-filtered by validation performance (GLUE/SUFI-2 convention).
# ------------------------------------------------------------------------------
band_full <- NULL
if (length(valid_runs) >= 2) {
  band_full <- calc_behavioural_band(comp, valid_runs)
  cat(sprintf("Full ensemble band      P-factor = %.3f  |  R-factor = %.3f\n",
              band_full$p_factor, band_full$r_factor))
}

# Optional: also compute the band over the still-satisfactory subset (diagnostic)
band_stable <- NULL
if (n_stable >= 2) {
  stable_runs <- results %>% filter(still_satisfactory) %>% pull(run_id)
  band_stable <- calc_behavioural_band(comp, stable_runs)
  cat(sprintf("Stable subset band      P-factor = %.3f  |  R-factor = %.3f\n",
              band_stable$p_factor, band_stable$r_factor))
}

# ------------------------------------------------------------------------------
# 12. WATER BALANCE DIAGNOSTIC (if targets + wb outputs are available)
# ------------------------------------------------------------------------------
balance_diag <- NULL
targets <- if (!is.null(val_cfg$targets)) val_cfg$targets else cfg$calibration$targets
has_balance_outputs <- all(c("precip", "et", "wateryld") %in% names(sim$simulation))

if (!is.null(targets) && has_balance_outputs) {
  balance_diag <- calc_balance_diagnostic(sim, valid_runs, targets)

  results <- results %>%
    left_join(balance_diag %>% select(run_id, ET_P, WYLD_P, err_ET_P_pct,
                                       err_WYLD_P_pct, agg_rel_error_pct),
              by = "run_id")

  cat("\n=== WATER BALANCE DIAGNOSTIC (validation window) ===\n")
  cat(sprintf("  Target ET/P  = %.4f  |  Sim median = %.4f  (range: %.4f - %.4f)\n",
              targets$et_rto,
              median(balance_diag$ET_P), min(balance_diag$ET_P), max(balance_diag$ET_P)))
  cat(sprintf("  Target WYLD/P = %.4f  |  Sim median = %.4f  (range: %.4f - %.4f)\n",
              targets$wyld_rto,
              median(balance_diag$WYLD_P), min(balance_diag$WYLD_P), max(balance_diag$WYLD_P)))
  cat(sprintf("  Median aggregate error = %.1f%%\n\n",
              median(balance_diag$agg_rel_error_pct)))
} else if (!is.null(targets) && !has_balance_outputs) {
  cat("NOTE: Water balance targets defined but precip/et/wateryld outputs not found.\n")
  cat("  Add precip, et, wateryld to cal_outputs (or validation.outputs) in config.yaml.\n\n")
}

# ------------------------------------------------------------------------------
# 13. PLOTS
# ------------------------------------------------------------------------------
nse_col <- paste0("NSE", suffix)
kge_col <- paste0("KGE", suffix)

# Scatter needs a column called 'satisfactory_run'
metrics_plot <- metrics %>% mutate(satisfactory_run = still_satisfactory)

g_scatter <- plot_performance_scatter(
  metrics_plot, nse_col, kge_col,
  val_thr$threshold_nse,
  val_thr$threshold_kge,
  paste("Validation — performance space (", temporal_sc, ")")
)
print(g_scatter)
save_tiff(g_scatter, "scatter_performance_val.tif", ver$path)

# Hydrograph uses the FULL-ensemble band so the calibration uncertainty is visible
stable_runs <- results %>% filter(still_satisfactory) %>% pull(run_id)

g_hydro <- plot_hydrograph(
  comp,
  satisfactory_runs = if (length(stable_runs) > 0) stable_runs else NULL,
  band              = if (!is.null(band_full)) band_full$band     else NULL,
  p_factor          = if (!is.null(band_full)) band_full$p_factor else NULL,
  r_factor          = if (!is.null(band_full)) band_full$r_factor else NULL,
  title_label       = paste("Validation hydrograph (", temporal_sc, ")")
)
print(g_hydro)
save_tiff(g_hydro, "hydrograph_val.tif", ver$path, width = 20, height = 10)

# Uncertainty envelope (SWAT-CUP style) over the full ensemble
if (!is.null(band_full)) {
  obs_trim <- obs_data %>% filter(Date %in% comp$Date)
  g_envelope <- plot_uncertainty_envelope(
    band        = band_full$band,
    obs_data    = obs_trim,
    p_factor    = band_full$p_factor,
    r_factor    = band_full$r_factor,
    title_label = paste("Validation — full-ensemble uncertainty band (", temporal_sc, ")")
  )
  print(g_envelope)
  save_tiff(g_envelope, "uncertainty_envelope_val.tif", ver$path, width = 20, height = 10)
}

# FDC envelope
if (!is.null(band_full)) {
  fdc_data <- calc_fdc_band(comp, valid_runs)
  if (!is.null(fdc_data)) {
    g_fdc <- plot_fdc_envelope(
      fdc_data,
      title_label = paste("Validation — Flow Duration Curve (", temporal_sc, ")")
    )
    print(g_fdc)
    save_tiff(g_fdc, "fdc_envelope_val.tif", ver$path, width = 20, height = 12)
    write.csv(fdc_data, file.path(ver$path, "fdc_band_val.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 14. SAVE RESULTS
# ------------------------------------------------------------------------------
saveRDS(
  list(
    config             = cfg,
    param_sets_source  = param_sets_path,
    run_map            = run_map,
    sim                = sim,
    results            = results,
    band_full          = band_full,
    band_stable        = band_stable,
    balance_diagnostic = balance_diag
  ),
  file.path(ver$path, "resultado_val.rds")
)

write.csv(results,  file.path(ver$path, "results_val.csv"), row.names = FALSE)
write.csv(run_map,  file.path(ver$path, "run_map.csv"),     row.names = FALSE)

if (!is.null(band_full)) {
  write.csv(band_full$band, file.path(ver$path, "behavioural_band_full.csv"),
            row.names = FALSE)
  write.csv(
    data.frame(
      subset   = "full_ensemble",
      n_runs   = length(valid_runs),
      p_factor = band_full$p_factor,
      r_factor = band_full$r_factor
    ),
    file.path(ver$path, "behavioural_band_factors_full.csv"), row.names = FALSE
  )
}

if (!is.null(band_stable)) {
  write.csv(band_stable$band, file.path(ver$path, "behavioural_band_stable.csv"),
            row.names = FALSE)
  write.csv(
    data.frame(
      subset   = "still_satisfactory",
      n_runs   = n_stable,
      p_factor = band_stable$p_factor,
      r_factor = band_stable$r_factor
    ),
    file.path(ver$path, "behavioural_band_factors_stable.csv"), row.names = FALSE
  )
}

if (!is.null(balance_diag)) {
  write.csv(balance_diag, file.path(ver$path, "balance_diagnostic_val.csv"),
            row.names = FALSE)
}

cat("\nValidation completed. Results saved to:\n ", ver$path, "\n")
