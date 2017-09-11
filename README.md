# serology-model
Inference framework for serological data


### Summary

This repository contains code for serological analysis of influenza. Currently set-up to run simulation study with [new model](http://www.biorxiv.org/content/early/2017/08/31/183111).

### Guide to files

`main_model.R` Main model fitting and plotting code - calls following source files:

a. `sero_functions.R` Functions for model and inference

b. `c_code/c_model.c` Serology model implemented in C

c. `posterior_analysis.R` Reads in MCMC output and plots posteriors and model estimates.