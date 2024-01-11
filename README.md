
<!-- README.md is generated from README.Rmd. Please edit that file -->

# serosolver

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`serosolver` is a modelling and inference package that uses a dynamic
model to infer antibody dynamics and infection histories from
cross-sectional or longitudinal serological data. The model infers
individual-level infection histories, historical attack rates, and
patterns of antibody dynamics by accounting for cross-reactive antibody
responses and measurement error.

## Installation

1.  Install [R](http://cran.r-project.org)

2.  Install the development version of serosolver from
    [GitHub](https://github.com/seroanalytics/serosolver):

``` r
devtools::install_github("seroanalytics/serosolver")
library(serosolver)
```

## Quick start and vignettes

Read the [quick start
vignette](https://seroanalytics.github.io/serosolver/articles/serosolver-quick_start_guide.html)
to set up and run a simple implementation with a simulation model.

There are additional Rmarkdown vignettes for [Case Study
1](https://seroanalytics.github.io/serosolver/articles/cs1_vignette.html)
(longitudinal analysis of influenza A/H1N1p in Hong Kong) and [Case
Study
2](https://seroanalytics.github.io/serosolver/articles/cs2_vignette.html)
(cross-sectional analysis of influenza A/H3N2 in Guangzhou Province,
China), to accompany the analysis in the serosolver paper.

## Example

This is a basic example of simulating some serological data and fitting
the model using the MCMC framework.

``` r
#library(serosolver)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
devtools::load_all("~/Documents/GitHub/serosolver")

## Load in example parameter values and antigenic map
data(example_par_tab)
data(example_antigenic_map)

## Get all possible infection times
possible_exposure_times <- unique(example_antigenic_map$inf_times)

## Vector of antibgens that have biomarker measurements (note only one representative antigen per time)
sampled_antigens <- seq(min(possible_exposure_times), max(possible_exposure_times), by=2)

## Times at which serum samples can be taken
sampling_times <- 2010:2015

## Number of serum samples taken
n_samps <- 2

## Simulate some random attack rates
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)

## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=example_par_tab, group=1, n_indiv=100,
                                  possible_exposure_times=possible_exposure_times,
                                  measured_biomarker_ids = sampled_antigens,
                                  sampling_times=2010:2015, nsamps=n_samps,
                                  antigenic_map=example_antigenic_map,
                                  age_min=10,age_max=75,
                                  attack_rates=attack_rates, repeats=2,
                                  obs_dist=c(1))

## Pull out the simulated titre data and infection histories
antibody_data <- all_simulated_data$antibody_data
example_inf_hist <- all_simulated_data$infection_histories
plot_data(antibody_data,example_inf_hist,possible_exposure_times,5)

## Run the MCMC
# This example uses prior version 2 (i.e. beta prior on phi with parameters alpha, beta)
# We have to remove the explicit specification of phi in the parameter table
res <- run_MCMC(example_par_tab, antibody_data, example_antigenic_map,
                filename="test", prior_version=2,
                mcmc_pars=c(adaptive_period=20000, iterations=80000),
                solve_likelihood=TRUE)

## Read in the MCMC chains and plot posteriors
chain <- read.csv(res$chain_file)
inf_chain <- data.table::fread(res$history_file)
plot(coda::as.mcmc(chain[chain$sampno > 20000,c("boost_long","wane_short","posterior_prob")]))

# Plot model predicted titres for a subset of individuals
plot_model_fits(chain = chain,infection_histories = inf_chain,
                         antibody_data = antibody_data,individuals=c(1:4),
                         antigenic_map=example_antigenic_map,par_tab=par_tab,orientation="cross-sectional")
```
