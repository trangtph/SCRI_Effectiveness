---
title: "SCRI Simulation - Results"
author: "Trang Tu"
format: docx
execute: 
  echo: false
  warning: false
  message: false
  output: true
  freeze: auto
editor: visual
---



```{r}
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
```

# 1. Base case

```{r}
Summary_base_case_20251223 <- read_excel("Results/Summary/Summary_base_case_20251223.xlsx")
```


```{r}
#| label: tbl-basecase
#| tbl-cap: "Performance metrics for base case scenarios" 

basecase_table <- Summary_base_case_20251223 %>% 
  mutate("Bias of est_V (MCSE)" = paste0(round(bias_est_V,3), " (", round(bias_est_V_MCSE,3), ")"),
         "Bias of VE (MCSE)" = paste0(round(bias_VE,3), " (", round(bias_VE_MCSE,3), ")"),
         "Empirical SE of est_V (MCSE)" = paste0(round(se_est_V_hat,3), " (", round(se_est_V_hat_MCSE,3), ")"),
         "Coverage of 95% CI of VE (MCSE)" = paste0(round(coverage_irr_V,3), " (", round(coverage_irr_V_MCSE,3), ")")) %>%
  select(sample_size, "Bias of est_V (MCSE)",  "Empirical SE of est_V (MCSE)", "Bias of VE (MCSE)",  "Coverage of 95% CI of VE (MCSE)") %>% rename("Cohort size" = sample_size)
  
kable(basecase_table)
```

