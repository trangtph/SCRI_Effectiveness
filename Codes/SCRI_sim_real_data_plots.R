if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(ggplot2, 
               dplyr,
               readr,
               tidyr,
               here,
               rio
)
options(scipen = 999)

source(here("Codes", "SCRI_sim_main_core_functions.R"))
source(here("Codes", "SCRI_sim_main_workflow_functions.R"))
source(here("Codes", "SCRI_helper_functions.R"))

# Get data ----

# Run first the script "SCRI_sim_real_data_execution.R"



# Some data manipulation ----
results_real <- results_real %>% 
  mutate(log_IRR_V_hat = log(IRR_V_hat),
         log_IRR_V_low_CI = log(IRR_V_hat_low_CI),
         log_IRR_V_up_CI = log(IRR_V_hat_up_CI),
         abs_relative_bias_estV = abs(relative_bias_est_V),
         abs_relative_bias_estV_low_CI = abs(relative_bias_est_V_low_CI),
         abs_relative_bias_estV_up_CI = abs(relative_bias_est_V_up_CI))

results_real <- results_real %>% mutate(
  cohort_size_lab = factor(cohort_size,
                           levels = c(20000, 50000, 100000, 500000),
                           labels = c("Cohort size: 20000","Cohort size: 50000","Cohort size: 100000", "Cohort size: 500000")),
  risk_end_id = case_when(risk_end == 21 ~ "days 15-21",
                        TRUE ~ "days 15-35")
) %>% mutate(across(c(risk_end_id, methods), as.factor)) %>%
  mutate(risk_end_id_num = as.numeric(risk_end_id), 
         method_offset = as.numeric(methods) * 0.15 - 0.3,   # Manual dodge: create offsets for each method
         y_dodged = risk_end_id_num + method_offset)


lollipop_plot_real <- function(data,
                           aes_x, aes_x_low_ci, aes_x_up_ci,
                           mode = c("est", "irr"),   # "est" = linear scale, "irr" = log scale
                           refline = 0,               # 0 for bias, 1 or 2 for IRR
                           irr_breaks = NULL,         # vector of IRR ticks, e.g. c(0.5,1,2,4)
                           irr_limits = NULL,         # two-element vector, IRR scale
                           x_break = NULL,            # for 'est' mode
                           x_limits = NULL,           # for 'est' mode
                           xlabel = "",
                           plot_name) {
  
  mode <- match.arg(mode)
  
  png(here("Plots", paste0(plot_name, ".png")), 
      width = 32, height = 18, units = "cm", res = 300)
  
  p <- ggplot(
    data,
    aes(x = .data[[aes_x]], y = y_dodged, color = risk_end_id)
  ) +
    geom_segment(aes(
      x = if (mode == "est") refline else log(refline),
      xend = .data[[aes_x]],
      y = y_dodged,
      yend = y_dodged
    ),
    linewidth = 1) +
    
    geom_point(aes(shape = methods), size = 3) +
    
    geom_errorbar(
      aes(xmin = .data[[aes_x_low_ci]], xmax = .data[[aes_x_up_ci]]),
      width = 0.2, alpha = 0.8,
      orientation = "y") +
    
    scale_y_continuous(
      breaks = unique(data$risk_end_id_num),
      labels = levels(data$risk_end_id)
    ) +
    
    scale_color_manual(values = c(
      "#bf3729","#2f357c", "#b0799a", "#e69b00", "#355828",
      "#6c5d9e", "#e48171", "#f5bb50",
      "#9d9cd5", "#17154f", "#f6b3b0", "#ada43b",
      "#1b9e77", "#4d4d4d", "#8c6d31")) + 
    guides(color = "none") + 
    
    scale_shape_manual(
      name = "Method",
      values = c(
        "no_calendar"  = 16,
        "calendar_7d"  = 17,
        "calendar_30d" = 15, 
        "calendar_7df3" = 23,
        "calendar_7df5" = 11
      ),
      labels = c(
        "no_calendar"  = "No calendar adjustment",
        "calendar_7d"  = "Calendar adj (7-day bin)",
        "calendar_30d" = "Calendar adj (30-day bin)",
        "calendar_7df3"= "Calendar adj (7-day from d3)",
        "calendar_7df5"= "Calendar adj (7-day from d5)"
      )
    ) +
    
    facet_wrap(~ cohort_size_lab, ncol = 4) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 11),
      strip.text   = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    labs(y = "Specification of protection window", x = xlabel) +
    coord_flip()
  
  # ---- SCALE LOGIC ----------------------------------------------------------
  
  if (mode == "est") {
    # linear scale
    p <- p +
      scale_x_continuous(
        breaks = x_break,
        limits = x_limits
      ) +
      geom_vline(xintercept = refline, linewidth = 1.2)
    
  } else if (mode == "irr") {
    # log scale, but axis shows IRR values
    p <- p +
      scale_x_continuous(
        breaks = log(irr_breaks),
        labels = irr_breaks,
        limits = log(irr_limits)
      ) +
      geom_vline(xintercept = log(refline), linewidth = 1.2)
  }
  
  print(p)
  dev.off()
}

#### Relative bias of est_V ----
lollipop_plot_real(data = results_real, aes_x ="relative_bias_est_V", 
               aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
               mode = "est",
               refline = 0,
               x_break = round(seq(from = 0, to = 7.5, by = 0.5),1),
               x_limits = c(0, 7.5),
               xlabel = "Relative bias of est_V",
               plot_name = "real_dat_bias_estV_relative_3models")

#### Estimated VE ----
# This is the VE corresponding to the average coefficient across replicates and its MCSE
lollipop_plot_real(data = results_real, aes_x ="VE_mean_est", 
               aes_x_low_ci ="VE_mean_est_low_CI", aes_x_up_ci = "VE_mean_est_up_CI",
               mode = "est",
               refline = 0.67,
               x_break = round(seq(from = 0.65, to = 1, by = 0.05),2),
               x_limits = c(0.65, 1),
               xlabel = "Estimated VE and 95% Monte Carlo CI",
               plot_name = "real_dat_VE_avg_estV")

#### Absolute bias of VE ----
lollipop_plot_real(data = results_real, aes_x ="bias_VE", 
               aes_x_low_ci ="bias_VE_low_CI", aes_x_up_ci = "bias_VE_up_CI",
               mode = "est",
               refline = 0,
               x_break = round(seq(from = -0.4, to = 0.2, by = 0.1),1),
               x_limits = c(-0.4, 0.2),
               xlabel = "Bias of VE",
               plot_name = "real_dat_bias_VE_3models")

#### Mean number of events ----

mean_events_plot <- function(data,
                             aes_x,
                             x_break = NULL,
                             x_limits = NULL,
                             xlabel = "",
                             plot_name) {
  
  png(here("Plots", paste0(plot_name, ".png")),
      width = 25, height = 18, units = "cm", res = 300)
  
  p <- ggplot(
    data,
    aes(x = .data[[aes_x]],
        y = risk_end_id,
        fill = risk_end_id)
  ) +
    
    geom_bar(stat = "identity") +
    
    geom_text(
      aes(label = round(.data[[aes_x]], 0)),
      vjust = 0.1,
      size = 5) + 
    
    scale_y_discrete(
      breaks = unique(data$risk_end_id),
      labels = levels(data$risk_end_id)
    ) +
    
    scale_fill_manual(values = c(
      "#bf3729", "#2f357c", "#b0799a", "#e69b00", "#355828",
      "#6c5d9e",  "#e48171", "#f5bb50",
      "#9d9cd5", "#17154f", "#f6b3b0", "#ada43b",
      "#1b9e77", "#4d4d4d", "#8c6d31")) +
    
    guides(fill = "none") +
    
    facet_wrap(~ cohort_size_lab, ncol = 4) +
    
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 11),
      strip.text   = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    
    labs(
      y = "Specification of protection window",
      x = xlabel
    ) +
    coord_flip()
  
  print(p)
  dev.off()
}

mean_events_plot(data = results_real[results_real$methods=="no_calendar",], aes_x ="mean_n_event", 
                 x_break = seq(from = 0, to = 1000, by = 50),
                 x_limits = c(0, 1000),
                 xlabel = "Mean number of events",
                 plot_name = "real_dat_mean_nr_events")

# -----------------------------------------------------------------------------
# Repeat the analysis for the test scenario of no seasonality of infection ----
# -----------------------------------------------------------------------------

## Some data manipulation ----
results_real_test <- results_real_test %>% 
  mutate(log_IRR_V_hat = log(IRR_V_hat),
         log_IRR_V_low_CI = log(IRR_V_hat_low_CI),
         log_IRR_V_up_CI = log(IRR_V_hat_up_CI),
         abs_relative_bias_estV = abs(relative_bias_est_V),
         abs_relative_bias_estV_low_CI = abs(relative_bias_est_V_low_CI),
         abs_relative_bias_estV_up_CI = abs(relative_bias_est_V_up_CI))

results_real_test <- results_real_test %>% mutate(
  cohort_size_lab = factor(cohort_size,
                           levels = c(20000, 50000, 100000, 500000),
                           labels = c("Cohort size: 20000","Cohort size: 50000","Cohort size: 100000", "Cohort size: 500000")),
  risk_end_id = case_when(risk_end == 21 ~ "days 15-21",
                          TRUE ~ "days 15-35")
) %>% mutate(across(c(risk_end_id, methods), as.factor)) %>%
  mutate(risk_end_id_num = as.numeric(risk_end_id), 
         method_offset = as.numeric(methods) * 0.15 - 0.3,   # Manual dodge: create offsets for each method
         y_dodged = risk_end_id_num + method_offset)


#### Relative bias of est_V ----
lollipop_plot_real(data = results_real_test, aes_x ="relative_bias_est_V", 
                   aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
                   mode = "est",
                   refline = 0,
                   x_break = round(seq(from = -0.5, to = 9.5, by = 0.5),1),
                   x_limits = c(-0.25, 9.6),
                   xlabel = "Relative bias of est_V",
                   plot_name = "real_test_bias_estV_relative_3models")

#### Estimated VE ----
# This is the VE corresponding to the average coefficient across replicates and its MCSE
lollipop_plot_real(data = results_real_test, aes_x ="VE_mean_est", 
                   aes_x_low_ci ="VE_mean_est_low_CI", aes_x_up_ci = "VE_mean_est_up_CI",
                   mode = "est",
                   refline = 0.67,
                   x_break = round(seq(from = 0.5, to = 1, by = 0.05),2),
                   x_limits = c(0.5, 1),
                   xlabel = "Estimated VE and 95% Monte Carlo CI",
                   plot_name = "real_test_VE_avg_estV")

#### Absolute bias of VE ----
lollipop_plot_real(data = results_real_test, aes_x ="bias_VE", 
                   aes_x_low_ci ="bias_VE_low_CI", aes_x_up_ci = "bias_VE_up_CI",
                   mode = "est",
                   refline = 0,
                   x_break = round(seq(from = -0.5, to = 0.1, by = 0.05),2),
                   x_limits = c(-0.55, 0.1),
                   xlabel = "Bias of VE",
                   plot_name = "real_test_bias_VE_3models")

#### Mean number of events ----
mean_events_plot(data = results_real_test[results_real_test$methods=="no_calendar",], aes_x ="mean_n_event", 
                 x_break = seq(from = 0, to = 1800, by = 100),
                 x_limits = c(0, 1800),
                 xlabel = "Mean number of events",
                 plot_name = "real_test_mean_nr_events")