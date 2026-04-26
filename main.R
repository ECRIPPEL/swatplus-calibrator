# ==============================================================================
# main.R — SWAT+ Calibrator: GSA → Calibration → (optional) Validation
#
# USAGE (RStudio / VSCode):
#   source("main.R")                   # uses config.yaml in working directory
#
#   CONFIG_PATH <- "config.yaml"       # or set path explicitly
#   source("main.R")
#
# USAGE (terminal):
#   Rscript main.R                     # uses config.yaml
#   Rscript main.R path/to/config.yaml # custom config path
#
# WORKFLOW:
#   Step 1 — run_gsa.R  : LH-OAT Morris sensitivity analysis
#                          → GSA_v1/ (ranking, plots, param_info_gsa.csv)
#   Step 2 — run_cal.R  : Monte Carlo calibration
#                          → CAL_v1/ (results, plots, new_ranges.csv)
#   Step 3 — run_val.R  : Split-sample validation (only if validation.enabled)
#                          → VAL_v1/ (results, band, plots, diagnostics)
#
# ITERATION:
#   After step 2, set in config.yaml:
#     iteration:
#       enabled: true
#       ranges_file: "path/to/CAL_v1/new_ranges.csv"
#   Then re-run main.R → creates GSA_v2/ and CAL_v2/
#
# VALIDATION:
#   Step 3 is skipped unless validation.enabled = true.
#   It can also be run in isolation: source("run_val.R") after CAL_v*/ exists.
# ==============================================================================

if (!exists("CONFIG_PATH")) {
  args        <- commandArgs(trailingOnly = TRUE)
  CONFIG_PATH <- if (length(args) >= 1) args[1] else "config.yaml"
}

cat("================================================================\n")
cat(" SWAT+ Calibrator\n")
cat(" Config :", CONFIG_PATH, "\n")
cat(" Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

.CALIBRATOR_STATE <- list(t_total = proc.time())

# ---- Step 1: GSA ----
cat("--- Step 1: Global Sensitivity Analysis ---\n")
.CALIBRATOR_STATE$t_step <- proc.time()
source("run_gsa.R")
elapsed_gsa <- round((proc.time() - .CALIBRATOR_STATE$t_step)["elapsed"] / 60, 1)
cat("GSA elapsed:", elapsed_gsa, "min\n\n")

# ---- Step 2: Calibration ----
cat("--- Step 2: Calibration ---\n")
.CALIBRATOR_STATE$t_step <- proc.time()
source("run_cal.R")
elapsed_cal <- round((proc.time() - .CALIBRATOR_STATE$t_step)["elapsed"] / 60, 1)
cat("Calibration elapsed:", elapsed_cal, "min\n\n")

# ---- Step 3: Validation (optional) ----
# Re-read the YAML because run_gsa/run_cal clear the environment.
.cfg_main <- yaml::read_yaml(CONFIG_PATH)
if (isTRUE(.cfg_main$validation$enabled)) {
  cat("--- Step 3: Validation ---\n")
  .CALIBRATOR_STATE$t_step <- proc.time()
  source("run_val.R")
  elapsed_val <- round((proc.time() - .CALIBRATOR_STATE$t_step)["elapsed"] / 60, 1)
  cat("Validation elapsed:", elapsed_val, "min\n\n")
} else {
  cat("Step 3 (validation) skipped (validation.enabled = false).\n\n")
}
rm(.cfg_main)

# ---- Summary ----
elapsed_total <- round((proc.time() - .CALIBRATOR_STATE$t_total)["elapsed"] / 60, 1)
cat("================================================================\n")
cat(" COMPLETED at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(" Total time :", elapsed_total, "min\n")
cat("================================================================\n")
