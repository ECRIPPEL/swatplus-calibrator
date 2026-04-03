# ==============================================================================
# plots.R — Publication-ready scientific plots (Morris, hydrograph, scatter, boxplot)
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

if (.Platform$OS.type == "windows") {
  windowsFonts(TimesNR = windowsFont("Times New Roman"))
}
BASE_FAMILY <- "TimesNR"

# ==============================================================================
# MORRIS SCREENING PLOT — mu* vs sigma
# ==============================================================================

#' Morris sensitivity screening plot
#'
#' @param ranking     Classified tibble with columns: parameter, mu_star, sigma, grupo
#' @param title_label Plot title suffix (e.g., "Iteration 1")
#' @return ggplot object
plot_morris <- function(ranking, title_label = "") {
  mu_med    <- mean(ranking$mu_star, na.rm = TRUE)
  sigma_med <- mean(ranking$sigma,   na.rm = TRUE)

  ranking <- ranking %>%
    mutate(plot_label = if_else(mu_star > 0 | sigma > 0, parameter, ""))

  ggplot(ranking, aes(x = mu_star, y = sigma)) +
    geom_hline(yintercept = sigma_med, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = mu_med,   linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(
      aes(fill = grupo, shape = grupo),
      color = "black", size = 3.8, stroke = 0.5, alpha = 0.8
    ) +
    geom_text_repel(
      aes(label = plot_label),
      family = BASE_FAMILY, size = 4, color = "black",
      box.padding = 0.5, point.padding = 0.3,
      segment.color = "grey70", segment.size = 0.4,
      min.segment.length = 0.2, max.overlaps = Inf
    ) +
    scale_fill_manual(values = c(
      "Low importance"                = "grey85",
      "High importance"               = "grey55",
      "High importance + interaction" = "black"
    )) +
    scale_shape_manual(values = c(
      "Low importance"                = 21,
      "High importance"               = 22,
      "High importance + interaction" = 24
    )) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      title = paste0("Morris sensitivity screening",
                     if (nchar(title_label) > 0) paste0(" — ", title_label) else ""),
      x = expression(mu^"*"),
      y = expression(sigma),
      fill = NULL, shape = NULL
    ) +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      plot.title      = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      axis.title      = element_text(size = 14, color = "black"),
      axis.text       = element_text(size = 12, color = "black"),
      legend.position = "top",
      legend.text     = element_text(size = 11),
      axis.line       = element_line(color = "black", linewidth = 0.5),
      axis.ticks      = element_line(color = "black")
    )
}

# ==============================================================================
# HYDROGRAPH — observed vs simulated (all runs + satisfactory highlighted)
# ==============================================================================

#' Hydrograph: observed vs simulated streamflow with behavioural band
#'
#' @param comparativo      data.frame with Date, Flow, and run_* columns
#' @param satisfactory_runs Character vector of satisfactory run IDs (NULL = no highlight)
#' @param band             data.frame with Date, lower, upper (NULL = no band)
#' @param p_factor         Numeric P-factor value (NULL = no annotation)
#' @param r_factor         Numeric R-factor value (NULL = no annotation)
#' @param title_label      Plot title
#' @return ggplot object
plot_hydrograph <- function(comparativo, satisfactory_runs = NULL,
                            band = NULL, p_factor = NULL, r_factor = NULL,
                            title_label = "Observed vs simulated streamflow") {
  sim_long <- comparativo %>%
    select(-Flow) %>%
    pivot_longer(cols = -Date, names_to = "run_id", values_to = "sim")

  g <- ggplot()

  # Layer 1: behavioural band (behind everything)
  if (!is.null(band)) {
    g <- g +
      geom_ribbon(
        data    = band,
        mapping = aes(x = Date, ymin = lower, ymax = upper),
        fill    = "#AED6F1", alpha = 0.45
      )
  }

  # Layer 2-3: simulation lines
  if (!is.null(satisfactory_runs)) {
    sim_long <- sim_long %>%
      mutate(group = ifelse(run_id %in% satisfactory_runs, "Satisfactory", "Non-satisfactory"))

    g <- g +
      geom_line(
        data    = sim_long %>% filter(group == "Non-satisfactory"),
        mapping = aes(x = Date, y = sim, group = run_id),
        colour  = "grey85", alpha = 0.60, linewidth = 0.22
      ) +
      geom_line(
        data    = sim_long %>% filter(group == "Satisfactory"),
        mapping = aes(x = Date, y = sim, group = run_id),
        colour  = "#5DADE2", alpha = 0.75, linewidth = 0.30
      )
  } else {
    g <- g +
      geom_line(
        data    = sim_long,
        mapping = aes(x = Date, y = sim, group = run_id),
        color   = "lightblue", linewidth = 0.4, alpha = 0.7
      )
  }

  # Layer 4: observed (on top)
  g <- g +
    geom_line(
      data    = comparativo,
      mapping = aes(x = Date, y = Flow),
      colour  = "black", linewidth = 0.45
    )

  # Layer 5: P-factor / R-factor annotation
  if (!is.null(p_factor) && !is.null(r_factor)) {
    g <- g +
      annotate(
        "label",
        x     = min(comparativo$Date) + 0.02 * as.numeric(diff(range(comparativo$Date))),
        y     = max(comparativo$Flow, na.rm = TRUE) * 0.95,
        label = paste0("P-factor = ", round(p_factor, 2), "\n",
                       "R-factor = ", round(r_factor, 2)),
        hjust = 0, vjust = 1,
        size  = 4, family = BASE_FAMILY,
        fill  = "white", alpha = 0.8,
        label.size = 0.3
      )
  }

  g +
    labs(title = title_label, x = NULL, y = expression(Streamflow ~ (m^3/s))) +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11, colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      plot.margin      = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# SCATTER — NSE vs KGE performance space
# ==============================================================================

#' Performance space scatter plot (NSE vs KGE)
#'
#' @param metrics       Tibble with NSE_*, KGE_*, satisfactory_run columns
#' @param nse_col       Name of the NSE column (e.g., "NSE_mon")
#' @param kge_col       Name of the KGE column (e.g., "KGE_mon")
#' @param threshold_nse NSE reference line
#' @param threshold_kge KGE reference line
#' @param title_label   Plot title
#' @return ggplot object
plot_performance_scatter <- function(metrics, nse_col, kge_col,
                                     threshold_nse, threshold_kge,
                                     title_label = "Performance space") {
  ggplot() +
    geom_point(
      data    = metrics %>% filter(!satisfactory_run),
      mapping = aes(x = .data[[nse_col]], y = .data[[kge_col]]),
      shape   = 21, size = 2.4, stroke = 0.20, color = "black", fill = "grey85", alpha = 0.75
    ) +
    geom_point(
      data    = metrics %>% filter(satisfactory_run),
      mapping = aes(x = .data[[nse_col]], y = .data[[kge_col]]),
      shape   = 21, size = 2.8, stroke = 0.22, color = "black", fill = "#6BAED6", alpha = 0.90
    ) +
    geom_vline(xintercept = threshold_nse, linetype = "dashed", color = "grey35") +
    geom_hline(yintercept = threshold_kge, linetype = "dashed", color = "grey35") +
    labs(title = title_label, x = "NSE", y = "KGE") +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text       = element_text(size = 11, colour = "black"),
      legend.position = "none"
    )
}

# ==============================================================================
# BOXPLOT — Normalized satisfactory parameter distributions
# ==============================================================================

#' Normalized boxplot of satisfactory parameter distributions
#'
#' @param satisfactory_results Tibble of satisfactory runs
#' @param param_info           Tibble with short, par_name, min, max
#' @param title_label          Plot title
#' @return ggplot object
plot_param_boxplot <- function(satisfactory_results, param_info,
                               title_label = "Satisfactory parameter distributions") {
  prior_df <- param_info %>%
    transmute(parameter = short, par_name, prior_min = min, prior_max = max)

  param_long <- satisfactory_results %>%
    select(all_of(param_info$par_name)) %>%
    pivot_longer(cols = everything(), names_to = "par_name", values_to = "value") %>%
    left_join(param_info %>% select(short, par_name), by = "par_name") %>%
    rename(parameter = short) %>%
    filter(!is.na(parameter), !is.na(value))

  param_norm <- param_long %>%
    left_join(prior_df, by = c("parameter", "par_name")) %>%
    mutate(
      value_norm = case_when(
        is.na(prior_min) | is.na(prior_max) ~ NA_real_,
        prior_max <= prior_min              ~ NA_real_,
        TRUE ~ (value - prior_min) / (prior_max - prior_min)
      )
    ) %>%
    filter(!is.na(value_norm))

  param_norm$parameter <- factor(param_norm$parameter, levels = param_info$short)

  ggplot(param_norm, aes(x = parameter, y = value_norm)) +
    geom_jitter(width = 0.22, alpha = 0.28, color = "black", size = 1.1) +
    geom_boxplot(
      width = 0.28, fill = "white", color = "black",
      linewidth = 0.40, outlier.shape = NA, alpha = 0.85
    ) +
    geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "grey45", linewidth = 0.45) +
    geom_hline(yintercept = 0.5,     linetype = "dotted", color = "grey65", linewidth = 0.40) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title_label, x = NULL, y = "Normalized parameter value") +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, face = "italic", size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.line    = element_line(color = "black", linewidth = 0.5),
      axis.ticks   = element_line(color = "black")
    )
}

# ==============================================================================
# UNCERTAINTY ENVELOPE — behavioural band + observed (SWAT-CUP style)
# ==============================================================================

#' Uncertainty envelope plot (behavioural band + observed)
#'
#' Clean plot showing only the shaded uncertainty band and the observed line,
#' similar to SWAT-CUP 95PPU plots but using the behavioural (satisfactory) band.
#'
#' @param band       data.frame with Date, lower, upper
#' @param obs_data   data.frame with Date, Flow
#' @param p_factor   Numeric P-factor value (NULL = no annotation)
#' @param r_factor   Numeric R-factor value (NULL = no annotation)
#' @param title_label Plot title
#' @return ggplot object
plot_uncertainty_envelope <- function(band, obs_data,
                                      p_factor = NULL, r_factor = NULL,
                                      title_label = "Behavioural uncertainty band") {

  g <- ggplot() +
    geom_ribbon(
      data    = band,
      mapping = aes(x = Date, ymin = lower, ymax = upper),
      fill    = "#2ECC71", alpha = 0.55
    ) +
    geom_line(
      data    = obs_data,
      mapping = aes(x = Date, y = Flow),
      colour  = "blue", linewidth = 0.50
    )

  # P-factor / R-factor annotation
  if (!is.null(p_factor) && !is.null(r_factor)) {
    g <- g +
      annotate(
        "label",
        x     = min(band$Date) + 0.02 * as.numeric(diff(range(band$Date))),
        y     = max(band$upper, obs_data$Flow, na.rm = TRUE) * 0.95,
        label = paste0("P-factor = ", round(p_factor, 2), "\n",
                       "R-factor = ", round(r_factor, 2)),
        hjust = 0, vjust = 1,
        size  = 4, family = BASE_FAMILY,
        fill  = "white", alpha = 0.8,
        label.size = 0.3
      )
  }

  g +
    labs(title = title_label, x = NULL, y = expression(Streamflow ~ (m^3/s))) +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11, colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      plot.margin      = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# FDC ENVELOPE — Flow Duration Curve with behavioural band
# ==============================================================================

#' Flow Duration Curve with uncertainty envelope
#'
#' Shows the observed FDC overlaid on the min-max envelope of satisfactory
#' simulations, plotted against exceedance probability (%).
#'
#' @param fdc_data    data.frame from calc_fdc_band (exceedance, obs_flow, lower, upper)
#' @param title_label Plot title
#' @return ggplot object
plot_fdc_envelope <- function(fdc_data,
                              title_label = "Flow Duration Curve — behavioural band") {

  ggplot(fdc_data) +
    geom_ribbon(
      aes(x = exceedance, ymin = lower, ymax = upper),
      fill = "#2ECC71", alpha = 0.55
    ) +
    geom_line(
      aes(x = exceedance, y = obs_flow),
      colour = "blue", linewidth = 0.55
    ) +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    labs(
      title = title_label,
      x     = "Exceedance probability (%)",
      y     = expression(Streamflow ~ (m^3/s))
    ) +
    theme_bw(base_family = BASE_FAMILY) +
    theme(
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11, colour = "black"),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# HELPER — Save plot as TIFF (publication-ready, 300 dpi, LZW compression)
# ==============================================================================

#' Save a ggplot as TIFF
#'
#' @param plot      ggplot object
#' @param filename  File name (e.g., "morris.tif")
#' @param save_path Destination folder
#' @param width     Width in cm  (default 18)
#' @param height    Height in cm (default 15)
save_tiff <- function(plot, filename, save_path, width = 18, height = 15) {
  ggsave(
    filename    = file.path(save_path, filename),
    plot        = plot,
    width       = width,
    height      = height,
    units       = "cm",
    dpi         = 300,
    compression = "lzw",
    device      = "tiff"
  )
}
