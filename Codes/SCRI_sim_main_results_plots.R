##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Script 4: Plots summarizing ###
### simulation results          ###
###################################

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

##############################
# 1 - Base case scenarios ----
##############################

## 1.1. Table of the scenarios' information ----
base_case <- data.frame(sample_size = c(6000, 10000, 20000))
base_case <-base_case %>% mutate(
  scen = "base_case",
  base_infect_shape = 2.5,
  base_infect_mode = 100,
  base_infect_min = 0.0002,
  base_infect_max = 0.002,
  vacc_season = "uniform",
  vacc_mean_d = NA,
  vacc_sd_d = NA, 
  control_start = 3,
  control_end = 7,
  risk_start = 16,
  risk_end = 35,
  scen_id = seq(1,3)
) %>% mutate(
  scen_name = paste0("S", scen_id, "_", scen, "_size", sample_size)
  
)

## 1.2. Only SCRI model without adjustment for seasonality was fitted ----
methods <- c("no_calendar")

## 1.3. Summarize the results ----
base_case_results <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                            scenario_table = base_case),
                                                  nsim = n_sim,
                                                  true_VE = 0.6,
                                                  results_dir = here("Results","Raw_results_base_case"),
                                                  summary_dir = file.path(here("Results"), "Summary"),
                                                  summary_file_name = "Summary_base_case_20251223")


## 1.4. Exploratory analysis ----

### 1.4.1. Extract data -----
load_raw_data <- function(method_scen_ = method_scen(method_table = as.data.frame(methods),
                                                        scenario_table = base_case),
                             results_dir = here("Results", "Raw_results_base_case")){
  scenarios <- unique(method_scen_$scen_name)
  methods <- unique(method_scen_$methods)
  
  all_data <- list()
  
  for (i in seq_len(nrow(method_scen_))) {
    method_i <- method_scen_$methods[i]
    scen_file <- paste0(method_scen_$scen_name[i], ".csv")
    
    file_path <- file.path(results_dir, method_i, scen_file)
    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    df <- tryCatch(read.csv(file_path),
                   error = function(e) NULL)
    
    if (is.null(df) ) {
      warning("No data in: ", file_path)
      next
    }
    
    all_data[[length(all_data) + 1]] <- data.frame(
      scen_name = method_scen_$scen_name[i],
      cohort_size = method_scen_$sample_size[i],
      est_V = df$est_V,
      se_est_V = df$se_V,
      n_event = df$n_event
    )
  }
  
  data_all <- bind_rows(all_data)
  return(data_all)
}

data_explor <- load_raw_data()
non_convergence <- data_explor %>% filter(est_V <= -10 | est_V >= 10 )


### 1.4.2. Histogram of est_V, n_event and se_V ---------------
#est_V: coefficient of vaccine effect
#se_V: SE of coefficient of vaccine effect

plot_exploratory <- function(
    data,
    output_dir = file.path(here(),"Plots"),
    width = 7,
    height = 6,
    dpi = 300
){
  scenarios <- unique(data$scen_name)
  
  output_dir1 = file.path(output_dir,"Hist_est_V")
  output_dir2 = file.path(output_dir,"Hist_se_V")
  output_dir3 = file.path(output_dir,"Hist_n_event")
  create_directory(output_dir1)
  create_directory(output_dir2)
  create_directory(output_dir3)
  
  
  for (sc in scenarios) {
    df_sc <- data %>% filter(scen_name == sc)
    
    if (nrow(df_sc) == 0) {
      warning("No data for scenario: ", sc)
      next
    }
    
    p <- ggplot(df_sc, aes(x = cohort_size, y = est_V, fill = cohort_size)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      labs(
        title = paste("Scenario:", sc),
        x = "Cohort Size",
        y = "est_V"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(output_dir1, paste0("plot_", sc, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
    
    p2 <- ggplot(df_sc, aes(x = cohort_size, y = se_est_V, fill = cohort_size)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      labs(
        title = paste("Scenario:", sc),
        x = "Cohort Size",
        y = "se_est_V"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(output_dir2, paste0("plot_", sc, ".png")),
      plot = p2,
      width = width,
      height = height,
      dpi = dpi
    )
    
    p3 <- ggplot(df_sc, aes(x = cohort_size, y = n_event, fill = cohort_size)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      labs(
        title = paste("Scenario:", sc),
        x = "Cohort Size",
        y = "n_event"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(output_dir3, paste0("plot_", sc, ".png")),
      plot = p3,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  message("All plots saved to: ", output_dir)
}

plot_exploratory(data = data_explor)

### 1.4.5. Agreement between empirical se_gamma_1 and model-based se_gamma_1 ---

se_agreement <- ggplot(data = base_case_results, aes(x = mod_se_est_V_hat, y = se_est_V_hat)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
se_agreement

## 1.5. Visualizing results ------

### 1.5.1. Some data manipulation for plotting IRR ----
base_case_results1 <- base_case_results %>% 
  mutate(log_IRR_V_hat = log(IRR_V_hat),
         log_IRR_V_low_CI = log(IRR_V_hat_low_CI),
         log_IRR_V_up_CI = log(IRR_V_hat_up_CI))

base_case_results2 <- base_case_results1 %>% mutate(across(c(sample_size), as.factor))

### 1.5.2 Lollipop plot ----

#### 1.5.2.1. Plot for bias ---
lollipop_plot <- function(data,
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
      width = 15, height = 8, units = "cm", res = 300)
  
  p <- ggplot(
    data,
    aes(x = .data[[aes_x]], y = sample_size, color = sample_size)
  ) +
    geom_segment(aes(
      x = if (mode == "est") refline else log(refline),
      xend = .data[[aes_x]],
      y = sample_size,
      yend = sample_size
    ),
    linewidth = 0.6) +
    
    geom_point(size = 1.5) +
    
    geom_errorbar(
      aes(xmin = .data[[aes_x_low_ci]], xmax = .data[[aes_x_up_ci]]),
      width = 0.3, alpha = 0.6,
      orientation = "y"
    ) +
    
    scale_y_discrete(
      breaks = unique(data$sample_size),
      labels = levels(data$sample_size)
    ) +
    
    scale_color_manual(
      values = c("#bf3729","#2f357c","#e69b00")
    ) +
    
    guides(color = "none") +
    
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12)
    ) +
    labs(y = "Cohort size", x = xlabel) +
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

# Absolute bias of est_V
lollipop_plot(data = base_case_results2, aes_x ="bias_est_V", 
              aes_x_low_ci ="bias_est_V_low_CI", aes_x_up_ci = "bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.1, by = 0.02),
              x_limits = c(-0.1, 0.1),
              xlabel = "Bias of est_V",
              plot_name = "basecase_bias_estV")

# Relative bias of est_V
lollipop_plot(data = base_case_results2, aes_x ="relative_bias_est_V", 
              aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.1, by = 0.02),
              x_limits = c(-0.1, 0.1),
              xlabel = "Relative bias of est_V",
              plot_name = "basecase_bias_estV_relative")

## Plot for variance ------------------------

lollipop_plot_var <- function(data, aes_x, aes_x_low_ci, aes_x_up_ci, plot_name){
  
  png(here("Plots", paste0(plot_name, ".png")), width = 15, height = 8, units = "cm", res = 300)
  
  p <- ggplot(data,
              aes(x = .data[[aes_x]],
                  y = sample_size,
                  color = sample_size)) +
    
    geom_point(size = 1.5) +
    
    geom_errorbar(aes(xmin = .data[[aes_x_low_ci]],
                       xmax = .data[[aes_x_up_ci]]),
                   width = 0.1,
                   alpha = 0.4,
                  orientation = "y") +
    
    scale_y_discrete(
      breaks = unique(data$sample_size),
      labels = levels(data$sample_size)
    ) +
    
    scale_x_continuous(
      limits = c(0, 0.8)) +
    
    scale_color_manual(
      values = c("#bf3729","#2f357c","#e69b00")
    ) +
    guides(color = "none") +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12)
    ) +
    labs(y = "Cohort size", x = "Empirical standard error of est_V ") +
    coord_flip()
  
  print(p)
  dev.off()
}

lollipop_plot_var(data = base_case_results2, aes_x ="se_est_V_hat", 
                  aes_x_low_ci ="se_est_V_hat_low_CI", aes_x_up_ci = "se_est_V_hat_up_CI",
                  plot_name = "basecase_empirical_se_est_V")



## Plot for mean number of events ----------------------------

nr_event_plot <- function(data, aes_x, plot_name){
  
  png(here("Plots", paste0(plot_name, ".png")), width = 8, height = 8, units = "cm", res = 300)
  
  p <- ggplot(data,
              aes(y = .data[[aes_x]],
                  x = sample_size,
                  fill = sample_size)) +
    
    geom_bar(stat = "identity") +
    
    scale_x_discrete(
      breaks = unique(data$sample_size),
      labels = levels(data$sample_size)
    ) +
    
    scale_y_continuous(breaks = seq(from = 0, to = 500, by = 50)) +
    
    scale_fill_manual(
      values = c("#bf3729","#2f357c","#e69b00")
    ) +
    guides(fill = "none") + 
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12)
    ) +
    labs(x = "Cohort size", y = "Mean number of events") #+
  #    coord_flip()
  
  print(p)
  dev.off()
}

nr_event_plot(data = base_case_results2, aes_x ="mean_n_event", 
              plot_name = "basecase_nr_event")

##############################
# 2 - Scenarios for bias quantification ----
##############################

## 2.1. Table of all scenarios ----
scen_table <- scenarios(vacc_seasonality = c("uniform", "beta"),
                        vacc_mean = c(180, 80),
                        vacc_sd = c(60, 20),
                        baseline_infection_shape = c(2.5,10, 20),
                        baseline_infection_mode = c(100, 300, 200),
                        baseline_infection_min = c(0.0002, 0.0002, 0.0002),
                        baseline_infection_max = c(0.002, 0.003, 0.003),
                        control_start_d = c(3),
                        control_end_d = c(15, 7),
                        risk_start_d = c(16, 8),
                        risk_end_d = c(35,77),
                        cohort_size = c(6000, 10000, 20000))

## 2.2. Three SCRI models  ----
methods <- c("no_calendar", "calendar_30d", "calendar_7d")

## 2.3.Summarize the results ---- 
results_all_scens <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                            scenario_table = scen_table),
                                                  nsim = n_sim,
                                                  true_VE = 0.6,
                                                  results_dir = here("Results","Raw_results_all_scens"),
                                                  summary_dir = file.path(here("Results"), "Summary"),
                                                  summary_file_name = "Summary_all_scens_20260412")

## 2.4. Plots ----

### 2.4.1 Some data manipulation ----
results_all_scens1 <- results_all_scens %>% 
  mutate(log_IRR_V_hat = log(IRR_V_hat),
         log_IRR_V_low_CI = log(IRR_V_hat_low_CI),
         log_IRR_V_up_CI = log(IRR_V_hat_up_CI),
         abs_relative_bias_estV = abs(relative_bias_est_V),
         abs_relative_bias_estV_low_CI = abs(relative_bias_est_V_low_CI),
         abs_relative_bias_estV_up_CI = abs(relative_bias_est_V_up_CI))

results_all_scens1 <- results_all_scens1 %>% mutate(across(c(sample_size), as.factor))

### 2.4.2. Misspecifying control window ----

misspecify_control <- results_all_scens1 %>% filter(scen == "misspecify_control", methods == "no_calendar")

### Absolute bias of est_V ----
lollipop_plot(data = misspecify_control, aes_x ="bias_est_V", 
              aes_x_low_ci ="bias_est_V_low_CI", aes_x_up_ci = "bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.3, by = 0.05),
              x_limits = c(-0.1, 0.3),
              xlabel = "Bias of est_V",
              plot_name = "mis_control_bias_estV")

#### Relative bias of est_V ---
lollipop_plot(data = misspecify_control, aes_x ="relative_bias_est_V", 
              aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.3, to = 0, by = 0.05),
              x_limits = c(-0.3, 0),
              xlabel = "Relative bias of est_V",
              plot_name = "mis_control_bias_estV_relative")

### Absolute bias of IRR ----
lollipop_plot(data = misspecify_control, aes_x ="bias_IRR_V", 
              aes_x_low_ci ="bias_IRR_V_low_CI", aes_x_up_ci = "bias_IRR_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.3, by = 0.05),
              x_limits = c(-0.1, 0.3),
              xlabel = "Bias of IRR",
              plot_name = "mis_control_bias_IRRV")

### Estimated VE ----
# This is the VE corresponding to the average coefficient across replicates and its MCSE
lollipop_plot(data = misspecify_control, aes_x ="VE_mean_est", 
              aes_x_low_ci ="VE_mean_est_low_CI", aes_x_up_ci = "VE_mean_est_up_CI",
              mode = "est",
              refline = 0.6,
              x_break = round(seq(from = 0.50, to = 0.6, by = 0.01),2),
              x_limits = c(0.50, 0.6),
              xlabel = "Estimated VE and \n95% Monte Carlo CI",
              plot_name = "mis_control_VE_avg_estV")


### Absolute bias of VE ----
lollipop_plot(data = misspecify_control, aes_x ="bias_VE", 
              aes_x_low_ci ="bias_VE_low_CI", aes_x_up_ci = "bias_VE_up_CI",
              mode = "est",
              refline = 0,
              x_break = round(seq(from = -0.3, to = 0.1, by = 0.05),2),
              x_limits = c(-0.3, 0.1),
              xlabel = "Bias of VE",
              plot_name = "mis_control_bias_VE")

### Variance ----
lollipop_plot_var(data = misspecify_control, aes_x ="se_est_V_hat", 
                  aes_x_low_ci ="se_est_V_hat_low_CI", aes_x_up_ci = "se_est_V_hat_up_CI",
                  plot_name = "mis_control_empirical_se_est_V")

# Reduction of SE compared to basecase
misspecify_control$se_est_V_hat/base_case_results2$se_est_V_hat

### Coverage ----
lollipop_plot(data = misspecify_control, aes_x ="coverage_irr_V", 
               aes_x_low_ci ="coverage_irr_V_low_CI", aes_x_up_ci = "coverage_irr_V_up_CI",
               mode = "est",
               refline = 0.95,
               x_break = round(seq(from = 0.2, to = 1, by = 0.1),1),
               x_limits = c(0.2, 1),
               xlabel = "Coverage of the IRR estimates",
               plot_name = "mis_control_coverage")

### 2.4.3. Misspecifying start of risk window ----

misspecify_risk_sta <- results_all_scens1 %>% filter(scen == "misspecify_risk_sta", methods == "no_calendar")

### Absolute bias of est_V ----
lollipop_plot(data = misspecify_risk_sta, aes_x ="bias_est_V", 
              aes_x_low_ci ="bias_est_V_low_CI", aes_x_up_ci = "bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.3, by = 0.05),
              x_limits = c(-0.1, 0.3),
              xlabel = "Bias of est_V",
              plot_name = "mis_risksta_bias_estV")

#### Relative bias of est_V ---
lollipop_plot(data = misspecify_risk_sta, aes_x ="relative_bias_est_V", 
              aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.3, to = 0, by = 0.05),
              x_limits = c(-0.3, 0),
              xlabel = "Relative bias of est_V",
              plot_name = "mis_risksta_bias_estV_relative")

### Absolute bias of IRR ----
lollipop_plot(data = misspecify_risk_sta, aes_x ="bias_IRR_V", 
              aes_x_low_ci ="bias_IRR_V_low_CI", aes_x_up_ci = "bias_IRR_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.3, by = 0.05),
              x_limits = c(-0.1, 0.3),
              xlabel = "Bias of IRR",
              plot_name = "mis_risksta_bias_IRRV")

### Estimated VE ----
# This is the VE corresponding to the average coefficient across replicates and its MCSE
lollipop_plot(data = misspecify_risk_sta, aes_x ="VE_mean_est", 
              aes_x_low_ci ="VE_mean_est_low_CI", aes_x_up_ci = "VE_mean_est_up_CI",
              mode = "est",
              refline = 0.6,
              x_break = round(seq(from = 0.50, to = 0.6, by = 0.01),2),
              x_limits = c(0.50, 0.6),
              xlabel = "Estimated VE and \n95% Monte Carlo CI",
              plot_name = "mis_risksta_VE_avg_estV")



### Absolute bias of VE ----
lollipop_plot(data = misspecify_risk_sta, aes_x ="bias_VE", 
              aes_x_low_ci ="bias_VE_low_CI", aes_x_up_ci = "bias_VE_up_CI",
              mode = "est",
              refline = 0,
              x_break = round(seq(from = -0.3, to = 0.1, by = 0.05),2),
              x_limits = c(-0.3, 0.1),
              xlabel = "Bias of VE",
              plot_name = "mis_risksta_bias_VE")

### Coverage ----
lollipop_plot(data = misspecify_risk_sta, aes_x ="coverage_irr_V", 
              aes_x_low_ci ="coverage_irr_V_low_CI", aes_x_up_ci = "coverage_irr_V_up_CI",
              mode = "est",
              refline = 0.95,
              x_break = round(seq(from = 0.2, to = 1, by = 0.1),1),
              x_limits = c(0.2, 1),
              xlabel = "Coverage of the IRR estimates",
              plot_name = "mis_risksta_coverage")
### Variance ----
lollipop_plot_var(data = misspecify_risk_sta, aes_x ="se_est_V_hat", 
                  aes_x_low_ci ="se_est_V_hat_low_CI", aes_x_up_ci = "se_est_V_hat_up_CI",
                  plot_name = "mis_risksta_empirical_se_est_V")

# Reduction of SE compared to basecase
1 - misspecify_risk_sta$se_est_V_hat/base_case_results2$se_est_V_hat

### 2.4.4. Misspecifying end of risk window ----

misspecify_risk_end <- results_all_scens1 %>% filter(scen == "misspecify_risk_end", methods == "no_calendar")

### Absolute bias of est_V ----
lollipop_plot(data = misspecify_risk_end, aes_x ="bias_est_V", 
              aes_x_low_ci ="bias_est_V_low_CI", aes_x_up_ci = "bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.3, by = 0.05),
              x_limits = c(-0.1, 0.3),
              xlabel = "Bias of est_V",
              plot_name = "mis_riskend_bias_estV")

#### Relative bias of est_V ---
lollipop_plot(data = misspecify_risk_end, aes_x ="relative_bias_est_V", 
              aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.3, to = 0, by = 0.05),
              x_limits = c(-0.3, 0),
              xlabel = "Relative bias of est_V",
              plot_name = "mis_riskend_bias_estV_relative")

### Estimated VE ----
# This is the VE corresponding to the average coefficient across replicates and its MCSE
lollipop_plot(data = misspecify_risk_end, aes_x ="VE_mean_est", 
              aes_x_low_ci ="VE_mean_est_low_CI", aes_x_up_ci = "VE_mean_est_up_CI",
              mode = "est",
              refline = 0.6,
              x_break = round(seq(from = 0.50, to = 0.6, by = 0.01),2),
              x_limits = c(0.50, 0.6),
              xlabel = "Estimated VE and \n95% Monte Carlo CI",
              plot_name = "mis_riskend_VE_avg_estV")


### Absolute bias of IRR ----
lollipop_plot(data = misspecify_risk_end, aes_x ="bias_IRR_V", 
              aes_x_low_ci ="bias_IRR_V_low_CI", aes_x_up_ci = "bias_IRR_V_up_CI",
              mode = "est",
              refline = 0,
              x_break = seq(from = -0.1, to = 0.3, by = 0.05),
              x_limits = c(-0.1, 0.3),
              xlabel = "Bias of VE",
              plot_name = "mis_riskend_bias_IRRV")

### Absolute bias of VE ----
lollipop_plot(data = misspecify_risk_end, aes_x ="bias_VE", 
              aes_x_low_ci ="bias_VE_low_CI", aes_x_up_ci = "bias_VE_up_CI",
              mode = "est",
              refline = 0,
              x_break = round(seq(from = -0.3, to = 0.1, by = 0.05),2),
              x_limits = c(-0.3, 0.1),
              xlabel = "Bias of VE",
              plot_name = "mis_riskend_bias_VE")

### Coverage ----
lollipop_plot(data = misspecify_risk_end, aes_x ="coverage_irr_V", 
              aes_x_low_ci ="coverage_irr_V_low_CI", aes_x_up_ci = "coverage_irr_V_up_CI",
              mode = "est",
              refline = 0.95,
              x_break = round(seq(from = 0.2, to = 1, by = 0.1),1),
              x_limits = c(0.2, 1),
              xlabel = "Coverage of the IRR estimates",
              plot_name = "mis_riskend_coverage")

### Variance ----
lollipop_plot_var(data = misspecify_risk_end, aes_x ="se_est_V_hat", 
                  aes_x_low_ci ="se_est_V_hat_low_CI", aes_x_up_ci = "se_est_V_hat_up_CI",
                  plot_name = "mis_riskend_empirical_se_est_V")

# Reduction of SE compared to basecase
1 - misspecify_risk_end$se_est_V_hat/base_case_results2$se_est_V_hat


### 2.4.5. Time-varying confounding ----------------------------------------

time_var <- results_all_scens1[results_all_scens1$scen == "seasonality",]

time_var <- time_var %>% mutate(
  cohort_size_lab = factor(sample_size,
                           levels = c(6000, 10000, 20000),
                           labels = c("Cohort size: 6000","Cohort size: 10000","Cohort size: 20000")),
  infect_dist_id = case_when(base_infect_shape == 2.5 ~ 1,
                             base_infect_shape == 20 ~ 2,
                             TRUE ~ 3),
  vacc_dist_id = case_when(vacc_mean_d == 180 & vacc_sd_d == 60 ~ 1,
                           vacc_mean_d == 180 & vacc_sd_d == 20 ~ 2,
                           vacc_mean_d == 80 & vacc_sd_d == 60 ~ 3,
                           vacc_mean_d == 80 & vacc_sd_d == 20 ~ 4, 
                           vacc_season == "uniform" ~ 0)
)

time_var <- time_var %>% mutate(
  season_id = paste0(infect_dist_id,"-", vacc_dist_id)) %>% arrange(season_id) %>%
  mutate(across(c(season_id), as.factor))

time_var <- time_var %>% mutate(
  season_id = paste0(infect_dist_id,"-", vacc_dist_id)) %>% arrange(season_id) %>%
  mutate(across(c(season_id, methods), as.factor)) %>%
  mutate(season_id_num = as.numeric(season_id), 
  method_offset = as.numeric(methods) * 0.15 - 0.30,   # Manual dodge: create offsets for each method
  y_dodged = season_id_num + method_offset)


lollipop_plot3 <- function(data,
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
      width = 30, height = 18, units = "cm", res = 300)
  
  p <- ggplot(
    data,
    aes(x = .data[[aes_x]], y = y_dodged, color = season_id)
  ) +
    geom_segment(aes(
      x = if (mode == "est") refline else log(refline),
      xend = .data[[aes_x]],
      y = y_dodged,
      yend = y_dodged
    ),
    linewidth = 0.6) +
    
    geom_point(aes(shape = methods), size = 1.5) +
    
    geom_errorbar(
      aes(xmin = .data[[aes_x_low_ci]], xmax = .data[[aes_x_up_ci]]),
      width = 0.3, alpha = 0.6,
      orientation = "y") +
    
    scale_y_continuous(
      breaks = unique(data$season_id_num),
      labels = levels(data$season_id)
    ) +
    
    scale_color_manual(values = c(
      "#2f357c", "#b0799a", "#e69b00", "#355828",
      "#6c5d9e", "#bf3729", "#e48171", "#f5bb50",
      "#9d9cd5", "#17154f", "#f6b3b0", "#ada43b",
      "#1b9e77", "#4d4d4d", "#8c6d31")) + 
    guides(color = "none") + 
    
    scale_shape_manual(
      name = "Method",
      values = c(
        "no_calendar"  = 16,
        "calendar_7d"  = 17,
        "calendar_30d" = 15
      ),
      labels = c(
        "no_calendar"  = "No calendar adjustment",
        "calendar_7d"  = "Calendar adjustment (7-day bin)",
        "calendar_30d" = "Calendar adjustment (30-day bin)"
      )
    ) +
    
    facet_wrap(~ cohort_size_lab, ncol = 2) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal"
      ) +
    labs(y = "Scenario of varying seasonality", x = xlabel) +
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

#### Absolute bias of est_V ----
lollipop_plot3(data = time_var, aes_x ="bias_est_V", 
               aes_x_low_ci ="bias_est_V_low_CI", aes_x_up_ci = "bias_est_V_up_CI",
               mode = "est",
               refline = 0,
               x_break = round(seq(from = -2.4, to = 0.7, by = 0.2),1),
               x_limits = c(-2.4, 0.7),
               xlabel = "Bias of est_V",
               plot_name = "seasonality_bias_estV_3models")

#### Relative bias of est_V ---
lollipop_plot3(data = time_var, aes_x ="relative_bias_est_V", 
               aes_x_low_ci ="relative_bias_est_V_low_CI", aes_x_up_ci = "relative_bias_est_V_up_CI",
               mode = "est",
               refline = 0,
               x_break = round(seq(from = -0.6, to = 2.6, by = 0.2),1),
               x_limits = c(-0.7, 2.6),
               xlabel = "Relative bias of est_V",
               plot_name = "seasonality_bias_estV_relative_3models")

#### Estimated VE ----
# This is the VE corresponding to the average coefficient across replicates and its MCSE
lollipop_plot3(data = time_var, aes_x ="VE_mean_est", 
              aes_x_low_ci ="VE_mean_est_low_CI", aes_x_up_ci = "VE_mean_est_up_CI",
              mode = "est",
              refline = 0.6,
              x_break = round(seq(from = 0.2, to = 1, by = 0.05),2),
              x_limits = c(0.2, 1),
              xlabel = "Estimated VE and 95% Monte Carlo CI",
              plot_name = "seasonality_VE_avg_estV")

#### Absolute bias of VE ----
lollipop_plot3(data = time_var, aes_x ="bias_VE", 
               aes_x_low_ci ="bias_VE_low_CI", aes_x_up_ci = "bias_VE_up_CI",
               mode = "est",
               refline = 0,
               x_break = round(seq(from = -1, to = 0.1, by = 0.1),1),
               x_limits = c(-1, 0.1),
               xlabel = "Bias of VE",
               plot_name = "seasonality_bias_VE_3models")

### Coverage ----
lollipop_plot3(data = time_var, aes_x ="coverage_irr_V", 
               aes_x_low_ci ="coverage_irr_V_low_CI", aes_x_up_ci = "coverage_irr_V_up_CI",
               mode = "est",
               refline = 0.95,
               x_break = round(seq(from = 0.2, to = 1, by = 0.1),1),
               x_limits = c(0.2, 1),
               xlabel = "Coverage of the IRR estimates",
               plot_name = "seasonality_coverage_3models")

### Mean number of events
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
        y = season_id,
        fill = season_id)
  ) +
    
    geom_bar(stat = "identity") +
    
    scale_y_discrete(
      breaks = unique(data$season_id),
      labels = levels(data$season_id)
    ) +
    
    scale_fill_manual(values = c(
        "#2f357c", "#b0799a", "#e69b00", "#355828",
        "#6c5d9e", "#bf3729", "#e48171", "#f5bb50",
        "#9d9cd5", "#17154f", "#f6b3b0", "#ada43b",
        "#1b9e77", "#4d4d4d", "#8c6d31")) +
    
    guides(fill = "none") +
    
    facet_wrap(~ cohort_size_lab, ncol = 2) +
    
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    
    labs(
      y = "Scenario of varying seasonality",
      x = xlabel
    ) +
    coord_flip()
  
  print(p)
  dev.off()
}

mean_events_plot(data = time_var[time_var$methods=="no_calendar",], aes_x ="mean_n_event", 
               x_break = seq(from = 0, to = 600, by = 100),
               x_limits = c(0, 600),
               xlabel = "Mean number of events",
               plot_name = "seasonality_mean_nr_events")