# The SCRI_Effectiveness_git repository

## Version 2025_12_24
This repository contains the code and simulation output of the sample size calculation and the main simulation of the project "Statistical Simulation for Application of Self-Control Risk Interval in Estimating Vaccine Effectiveness" 

It contains the R code, simulation workflow, and results of 1) the preliminary simulation to estimate the required sample size to detect a specified vaccine effectiveness using the SCRI method, and 2) the main simulation.

## Project organization

```
.
├── .gitignore
├── README.md
├── Codes                              <- R scripts for simulation and analysis
|   ├── SCRI_helper_functions.R            <- Utility functions used across scripts
|   ├── SCRI_other_plots.R                 <- Script to produce plots of distributions of infection risk and vaccination date, and sample size calculation
|   ├── SCRI_sim_main_core_functions.R     <- Core functions for data generation, SCRI modelling and summarize the results
|   ├── SCRI_sim_main_execution.R          <- Script to execute the full simulation workflow
|   ├── SCRI_sim_main_results_plots.R      <- Script to produce plots to summarize the simulation results
|   ├── SCRI_sim_main_results_table.QMD    <- Script to produce tables to summarize the simulation results
|   ├── SCRI_sim__main_workflow_functions.R  <- Workflow functions for running large-scale simulations
|   ├── SCRI_sim_samplesize_core_functions.R <- Core functions for data generation, SCRI modelling and summarize the results for sample size calculation
|   ├── SCRI_sim_samplesize_execution.R    <- Script to execute the full simulation workflow for sample size calculation
|   └── SCRI_sim_samplesize_workflow_functions.R   <- Workflow functions for running large-scale simulations for sample size calculation  

├── Results
│   ├── Raw_results_all_scens          <- Output from each simulation run of 45 simulation scenarios
│   ├── Raw_sample_base_Case           <- Output from each simulation run of 3 base case scenarios
│   ├── Raw_sample_size                <- Output from each simulation run of different cohort size (for sample size calculation)
│   └── Summary                        <- Aggregated summary of simulation results
└── Plots                              <- Visualizations for reporting and protocol

```

## Simulation design

### Sample size calculation

Objective: Determine the minimum sample size to achieve 80% power for detecting a 60% vaccine effectiveness, assuming a two-sided α = 0.05.

Setup:

- Sample sizes: 10 cohort sizes from 1,000 to 10,000 individuals (increment = 1,000)
- Number of simulations: 1,000 per sample size
- Two SCRI Models are fitted on vaccinated cases data:
  - Unadjusted SCRI
  - SCRI adjusted for calendar month

- Output: 
  - Empirical power, computed as the proportion of simulations with a vaccine effect p-value < 0.05
  - Absolute bias of the IRR estimates
 
### Main simulation - base case scenarios
Setup:
- Sample sizes: three cohort of size 6000, 10000 and 20000
- Number of simulations: 1,000 per sample size
- Data generating mechanism:
    + Seasonality of infection risk: Gamma distribution, lowest and highest daily risk being 2*10^-4 and 2*10^-3 , respectively
    + Seasonality of vaccination date: none
    + Effectiveness of vaccine: 
        + No effect: day 0-7 and from day 151
        + Peak at 60% on days 16-35
- SCRI model: focal window day 16-35, control window day 3-7, no adjustment for time-varying confounding

### Main simulation - scenarios to quantify bias

The data-generating mechanism and specification of SCRI models for 15 scenarios to quantify bias are described in this document: https://docs.google.com/document/d/1Hr-czYjpjXcnIFoIeJYPfNhpYMxBlLBW5Mz22LbAd-c/edit?usp=sharing. 
Each scenario is run with three cohort sizes  6000, 10000 and 20000.


## Reproducibility
To reproduce the simulation workflow:
1. Clone the repository

```
git clone https://github.com/trangtph/SCRI_Effectiveness.git
```
2. Open `SCRI_sim_main_execution.R`, adjust the paths to the source codes and the desired output location as needed, and run the script in your R environment. 
