# serosolver
Inference framework for serological data


### Summary

This repository contains code for serological analysis of influenza. Currently set-up to run simulation study with [new model](http://www.biorxiv.org/content/early/2017/08/31/183111).

### Guide to files

`scripts/main_model.R` Simulates some artificial data and attempts to re-estimate the simulation parameters, or attempts to estimate posterior distributions for model parameters given some real data
`scripts/extracting_data.R` extracts the fluscape data into a useable format
`scripts/extract_fluscape_antigenicmap.R` extracts the fluscape antigenic map

### Helpful pointers
All of the functions in `plots.R` are useful for looking at output. They currently only work with a single MCMC chain - priority is to change this to take an MCMC list