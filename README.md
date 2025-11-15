# The SCRI_Effectiveness_git repository

## Version 2025_11_15
This repository contains the code and simulation output of the sample size calculation section for the project "Statistical Simulation for Application of Self-Control Risk Interval in Estimating Vaccine Effectiveness" 

It contains the R code, simulation workflow, and results used to estimate the required sample size to detect a specified vaccine effectiveness using the SCRI method.

## Project organization

```
.
├── .gitignore
├── README.md
├── Codes                        <- R scripts for simulation and analysis
│   ├── SCRI_main_functions.R           <- Core functions for data generation and SCRI modeling
│   ├── SCRI_sim_workflow_functions.R  <- Workflow functions for running large-scale simulations
│   ├── SCRI_helper_functions.R        <- Utility functions used across scripts
│   └── SCRI_sim_execution.R           <- Script to execute the full simulation workflow
├── Results
│   ├── Raw_sample_size          <- Output from each simulation run
│   └── Summary                  <- Aggregated summary of simulation results
└── Plots                        <- Visualizations for reporting and protocol

```

## Simulation design

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

## Reproducibility
To reproduce the simulation workflow:
1. Clone the repository

```
git clone https://github.com/trangtph/SCRI_Effectiveness.git
```
2. Open `SCRI_sim_execution.R`, adjust the paths to the source codes and the desired output location as needed, and run the script in your R environment. 
