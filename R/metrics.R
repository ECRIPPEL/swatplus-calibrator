# ==============================================================================
# metrics.R — Performance metrics for streamflow calibration
# ==============================================================================

library(dplyr)
library(tibble)
library(purrr)
library(hydroGOF)

#' Compute NSE, KGE, PBIAS, and obj_total for each run
#'
#' @param comparativo data.frame with Date, Flow (observed), and run_* columns
#' @param valid_runs  Character vector of run names
#' @param suffix      Metric column suffix (e.g., "_mon" or "_day")
#' @return Tibble with metrics per run
calc_hard_metrics <- function(comparativo, valid_runs, suffix = "_mon") {
  map_dfr(valid_runs, function(run_nm) {
    sim_vals <- comparativo[[run_nm]]
    obs_vals <- comparativo$Flow

    if (all(is.na(sim_vals)) || length(na.omit(sim_vals)) < 2 ||
        sd(sim_vals, na.rm = TRUE) == 0) {
      out <- tibble(
        run_id    = run_nm,
        NSE       = NA_real_,
        KGE       = NA_real_,
        PBIAS     = NA_real_,
        absPBIAS  = NA_real_,
        err_nse   = NA_real_,
        err_kge   = NA_real_,
        err_pbias = NA_real_,
        obj_total = NA_real_,
        run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_nm))
      )
      names(out) <- gsub("^(NSE|KGE|PBIAS)$", paste0("\\1", suffix), names(out))
      return(out)
    }

    nse_n   <- as.numeric(NSE(sim = sim_vals,   obs = obs_vals))
    kge_n   <- as.numeric(KGE(sim = sim_vals,   obs = obs_vals))
    pbias_n <- as.numeric(pbias(sim = sim_vals, obs = obs_vals))

    out <- tibble(
      run_id    = run_nm,
      NSE       = nse_n,
      KGE       = kge_n,
      PBIAS     = pbias_n,
      absPBIAS  = abs(pbias_n),
      err_nse   = pmax(0, 1 - pmin(nse_n, 1)),
      err_kge   = pmax(0, 1 - pmin(kge_n, 1)),
      err_pbias = abs(pbias_n) / 100,
      obj_total = pmax(0, 1 - pmin(nse_n, 1)) + pmax(0, 1 - pmin(kge_n, 1)) + abs(pbias_n) / 100,
      run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_nm))
    )

    names(out) <- gsub("^(NSE|KGE|PBIAS)$", paste0("\\1", suffix), names(out))
    out
  })
}

#' Classify runs as satisfactory based on NSE / KGE / PBIAS thresholds
#'
#' @param metrics             Tibble of metrics (output of calc_hard_metrics)
#' @param threshold_nse       Minimum NSE
#' @param threshold_kge       Minimum KGE
#' @param threshold_abs_pbias Maximum |PBIAS|
#' @param suffix              Column suffix (e.g., "_mon" or "_day")
#' @return Tibble with added column satisfactory_run
classify_hard_runs <- function(metrics, threshold_nse, threshold_kge,
                               threshold_abs_pbias, suffix = "_mon") {
  nse_col <- paste0("NSE", suffix)
  kge_col <- paste0("KGE", suffix)

  metrics %>%
    mutate(
      satisfactory_run =
        !is.na(.data[[nse_col]]) &
        !is.na(.data[[kge_col]]) &
        .data[[nse_col]] >= threshold_nse &
        .data[[kge_col]] >= threshold_kge &
        absPBIAS        <= threshold_abs_pbias
    ) %>%
    arrange(desc(satisfactory_run), desc(.data[[nse_col]]), desc(.data[[kge_col]]), absPBIAS)
}

#' Compute the behavioural band (min–max envelope) and P-factor / R-factor
#'
#' The band covers the full range of satisfactory simulations at each timestep.
#' P-factor = fraction of observed points bracketed by the band.
#' R-factor = mean band width / sd(observed).
#'
#' @param comparativo      data.frame with Date, Flow, and run_* columns
#' @param satisfactory_runs Character vector of satisfactory run IDs
#' @return List with $band (Date, lower, upper), $p_factor, $r_factor.
#'         Returns NULL if fewer than 2 satisfactory runs.
calc_behavioural_band <- function(comparativo, satisfactory_runs) {
  if (length(satisfactory_runs) < 2) {
    warning("Behavioural band requires at least 2 satisfactory runs.")
    return(NULL)
  }

  sat_mat <- as.matrix(comparativo[, satisfactory_runs, drop = FALSE])

  lower <- apply(sat_mat, 1, min, na.rm = TRUE)
  upper <- apply(sat_mat, 1, max, na.rm = TRUE)

  band <- data.frame(
    Date  = comparativo$Date,
    lower = lower,
    upper = upper
  )

  obs <- comparativo$Flow
  n_valid   <- sum(!is.na(obs))
  bracketed <- sum(obs >= lower & obs <= upper, na.rm = TRUE)

  p_factor <- bracketed / n_valid
  r_factor <- mean(upper - lower, na.rm = TRUE) / sd(obs, na.rm = TRUE)

  list(band = band, p_factor = p_factor, r_factor = r_factor)
}

#' Compute Flow Duration Curve (FDC) envelope from satisfactory simulations
#'
#' Sorts each satisfactory run and observed data by exceedance probability,
#' then builds the min-max envelope across runs at each probability level.
#'
#' @param comparativo      data.frame with Date, Flow, and run_* columns
#' @param satisfactory_runs Character vector of satisfactory run IDs
#' @return data.frame with columns: exceedance (0-100%), obs_flow, lower, upper.
#'         Returns NULL if fewer than 2 satisfactory runs.
calc_fdc_band <- function(comparativo, satisfactory_runs) {
  if (length(satisfactory_runs) < 2) {
    warning("FDC band requires at least 2 satisfactory runs.")
    return(NULL)
  }

  obs <- sort(comparativo$Flow[!is.na(comparativo$Flow)], decreasing = TRUE)
  n   <- length(obs)
  exc <- 100 * seq_len(n) / (n + 1)

  sat_mat <- as.matrix(comparativo[, satisfactory_runs, drop = FALSE])
  sat_sorted <- apply(sat_mat, 2, function(x) sort(x[!is.na(x)], decreasing = TRUE)[seq_len(n)])

  lower <- apply(sat_sorted, 1, min, na.rm = TRUE)
  upper <- apply(sat_sorted, 1, max, na.rm = TRUE)

  data.frame(
    exceedance = exc,
    obs_flow   = obs,
    lower      = lower,
    upper      = upper
  )
}

#' Compute water balance diagnostic (ET/P and WYLD/P) for valid runs
#'
#' Compares simulated ratios against observed targets.
#' Requires "precip", "et", and "wateryld" outputs in sim$simulation.
#'
#' @param sim         SWATrunR simulation object
#' @param valid_runs  Character vector of run names
#' @param targets     List with et_rto and wyld_rto
#' @return Tibble with run_id, ET_P, WYLD_P, target values, and relative errors
calc_balance_diagnostic <- function(sim, valid_runs, targets) {
  tibble(
    run_id   = valid_runs,
    P_sim    = as.numeric(sim$simulation$precip[nrow(sim$simulation$precip),     valid_runs]),
    ET_sim   = as.numeric(sim$simulation$et[nrow(sim$simulation$et),             valid_runs]),
    WYLD_sim = as.numeric(sim$simulation$wateryld[nrow(sim$simulation$wateryld), valid_runs])
  ) %>%
    mutate(
      ET_P             = ET_sim / P_sim,
      WYLD_P           = WYLD_sim / P_sim,
      target_ET_P      = targets$et_rto,
      target_WYLD_P    = targets$wyld_rto,
      err_ET_P_pct     = 100 * abs(ET_P - target_ET_P) / target_ET_P,
      err_WYLD_P_pct   = 100 * abs(WYLD_P - target_WYLD_P) / target_WYLD_P,
      agg_rel_error_pct = err_ET_P_pct + err_WYLD_P_pct
    )
}
