# ==============================================================================
# morris_classify.R — Morris screening classification + quantile-based range narrowing
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

#' Classify parameters by Morris group (mu_star vs sigma)
#'
#' Groups:
#'   - High importance + interaction : mu_star > mean AND sigma > mean
#'   - High importance               : mu_star > mean AND sigma <= mean
#'   - Low importance                : mu_star <= mean
#'
#' @param ranking Tibble with columns: parameter, mu_star, sigma
#' @return Tibble with added column 'grupo'
classify_morris <- function(ranking) {
  mu_med    <- mean(ranking$mu_star, na.rm = TRUE)
  sigma_med <- mean(ranking$sigma,   na.rm = TRUE)

  ranking %>%
    mutate(
      grupo = case_when(
        mu_star > mu_med & sigma > sigma_med  ~ "High importance + interaction",
        mu_star > mu_med & sigma <= sigma_med ~ "High importance",
        TRUE                                  ~ "Low importance"
      ),
      grupo = factor(
        grupo,
        levels = c("Low importance", "High importance", "High importance + interaction")
      )
    )
}

#' Quantile narrowing rules per Morris group
#'
#' @return Tibble with grupo, q_low, q_high
morris_quantile_rules <- function() {
  tibble(
    grupo  = c("High importance + interaction", "High importance", "Low importance", "Not screened"),
    q_low  = c(0.05, 0.10, 0.25, 0.05),
    q_high = c(0.95, 0.90, 0.75, 0.95)
  )
}

#' Compute new parameter ranges from satisfactory posterior + Morris classification
#'
#' @param satisfactory_results Tibble of satisfactory runs (with parameter columns)
#' @param param_names_swat     Character vector of SWATrunR par_name columns
#' @param param_info           Tibble with short and par_name columns
#' @param ranking_morris       Tibble with parameter, mu_star, sigma (from GSA)
#' @return Tibble with new_min, new_max per parameter
calc_new_ranges <- function(satisfactory_results, param_names_swat,
                            param_info, ranking_morris) {
  quantile_rules <- morris_quantile_rules()

  mu_med    <- mean(ranking_morris$mu_star, na.rm = TRUE)
  sigma_med <- mean(ranking_morris$sigma,   na.rm = TRUE)

  morris_df <- ranking_morris %>%
    transmute(
      parameter,
      mu_star,
      sigma,
      grupo = case_when(
        mu_star > mu_med & sigma > sigma_med  ~ "High importance + interaction",
        mu_star > mu_med & sigma <= sigma_med ~ "High importance",
        TRUE                                  ~ "Low importance"
      )
    )

  param_long_sat <- satisfactory_results %>%
    select(all_of(param_names_swat)) %>%
    pivot_longer(
      cols      = everything(),
      names_to  = "par_name",
      values_to = "value"
    ) %>%
    left_join(param_info %>% select(short, par_name), by = "par_name") %>%
    rename(parameter = short)

  new_ranges <- param_long_sat %>%
    left_join(morris_df, by = "parameter") %>%
    mutate(grupo = ifelse(is.na(grupo), "Not screened", grupo)) %>%
    left_join(quantile_rules, by = "grupo") %>%
    group_by(parameter, grupo, q_low, q_high) %>%
    summarise(
      n_runs   = n(),
      post_min = round(min(value,                                  na.rm = TRUE), 3),
      post_max = round(max(value,                                  na.rm = TRUE), 3),
      new_min  = round(quantile(value, probs = first(q_low),  na.rm = TRUE, type = 8), 3),
      new_max  = round(quantile(value, probs = first(q_high), na.rm = TRUE, type = 8), 3),
      mediana  = round(median(value,                               na.rm = TRUE), 3),
      .groups  = "drop"
    ) %>%
    select(parameter, grupo, n_runs, post_min, post_max, q_low, q_high, new_min, new_max, mediana) %>%
    arrange(
      factor(grupo, levels = c("High importance + interaction", "High importance",
                               "Low importance", "Not screened")),
      desc(new_max - new_min)
    )

  new_ranges
}
