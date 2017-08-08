# serology-model
Inference framework for serological data


### Summary

This repository contains code for serological analysis of influenza.

### Guide to files

`main_model.R` Main model fitting and plotting code - calls following source files:

a. `load_data.R` Loads data from Vietnam and China and reshapes to get in format efficient for simulation/inference

b. `sero_functions.R` Functions for model and inference

c. `c_code/c_model.c` Serology model implemented in C

d. `posterior_analysis.R` Reads in MCMC output and plots posteriors and model estimates.