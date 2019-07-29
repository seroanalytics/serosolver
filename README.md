
<!-- README.md is generated from README.Rmd. Please edit that file -->
serosolver
==========

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`serosolver` is a modelling and inference package that uses a dynamic model to infer antibody dynamics and infection histories from cross-sectional or longitudinal serological data. The model infers individual-level infection histories, historical attack rates, and patterns of antibody dynamics by accounting for cross-reactive antibody responses and measurement error.

Installation
------------

1.  Install [R](http://cran.r-project.org)

2.  Install the development version of serosolver from [GitHub](https://github.com/seroanalytics/serosolver):

``` r
devtools::install_github("seroanalytics/serosolver")
library(serosolver)
```

Quick start
-----------

Read the [quick start vignette](https://github.com/seroanalytics/serosolver/blob/master/vignettes/serosolver-quick_start_guide.md) to set up and run a simple implementation with a simulation model.

Example
-------

This is a basic example of simulating some serological data and fitting the model using the MCMC framework.

``` r
library(serosolver)
library(plyr)
## Load in example parameter values and antigenic map
data(example_par_tab)
data(example_antigenic_map)

## Get all possible infection times
strain_isolation_times <- unique(example_antigenic_map$inf_years)

## Vector of strains that have titres (note only one representative strain per time)
sampled_viruses <- seq(min(strain_isolation_times), max(strain_isolation_times), by=2)

## Times at which serum samples can be taken
sampling_times <- 2010:2015

## Number of serum samples taken
n_samps <- 2

## Simulate some random attack rates
attack_rates <- runif(length(strain_isolation_times), 0.05, 0.15)

## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=example_par_tab, group=1, n_indiv=50,
                                  strain_isolation_times=strain_isolation_times,
                                  measured_strains=sampled_viruses,
                                  sampling_times=2010:2015, nsamps=n_samps,
                                  antigenic_map=example_antigenic_map,
                                  age_min=10,age_max=75,
                                  attack_rates=attack_rates, repeats=2)

## Pull out the simulated titre data and infection histories
titre_dat <- all_simulated_data$data
ages <- all_simulated_data$ages
example_inf_hist <- all_simulated_data$infection_histories
example_titre_dat <- merge(titre_dat, ages)

## Run the MCMC. We have to remove phi, as running prior version 2.
par_tab <- example_par_tab[example_par_tab$names != "phi",]
res <- run_MCMC(par_tab, example_titre_dat, example_antigenic_map, 
                filename="test", version=2,
                mcmc_pars=c(save_block=10000,thin=10,thin_hist=100,
                            move_sizes=3,swap_propn=0.5,inf_propn=0.5))

## Read in the MCMC chains and plot!
chain <- read.csv(res$chain_file)
inf_chain <- data.table::fread(res$history_file)
plot(coda::as.mcmc(chain[chain$sampno > 10000,c("mu","mu_short","wane")]))
```
