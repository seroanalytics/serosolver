## Script to run the entire README pipeline, simulating the data, saving the example data structures, and testing the MCMC runs
set.seed(1234)
library(tidyverse)
devtools::document("~/Documents/GitHub/serosolver")
devtools::load_all("~/Documents/GitHub/serosolver")

################################################
## 1. HA NAM DATA
################################################
## Get all possible infection times
antigenic_map <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/antigenic_maps/antigenicMap_vietnam.csv")
par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")

possible_exposure_times <- antigenic_map$inf_times

## Vector of antigens that have biomarker measurements (note only one representative antigen per time)
sampled_antigens <- seq(min(possible_exposure_times), max(possible_exposure_times), by=2)

## Times at which serum samples can be taken
sampling_times <- 2010:2015

## Number of serum samples taken
n_samps <- 5

## Simulate some random attack rates
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)


## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=50,
                                    possible_exposure_times=possible_exposure_times,
                                    measured_biomarker_ids = sampled_antigens,
                                    sampling_times=sampling_times, nsamps=n_samps,
                                    antigenic_map=antigenic_map,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats=2,
                                    data_type=c(1))

## Pull out the simulated titre data and infection histories
antibody_data <- all_simulated_data$antibody_data
true_inf_hist <- all_simulated_data$infection_histories
plot_antibody_data(antibody_data,possible_exposure_times,1:4,infection_histories = true_inf_hist)
example_par_tab <- par_tab
example_antibody_data <- antibody_data
example_inf_hist <- true_inf_hist
example_antigenic_map <- antigenic_map

save(example_antigenic_map,file="~/Documents/GitHub/serosolver/data/example_antigenic_map.RData")
save(example_par_tab,file="~/Documents/GitHub/serosolver/data/example_par_tab.RData")
save(example_antibody_data,file="~/Documents/GitHub/serosolver/data/example_antibody_data.RData")
save(example_inf_hist,file="~/Documents/GitHub/serosolver/data/example_inf_hist.RData")

output <- serosolver(example_par_tab, example_antibody_data, example_antigenic_map,
                                 filename="readme", n_chains=3,parallel=TRUE,
                                 mcmc_pars=c(adaptive_iterations=1000, iterations=5000),verbose=FALSE)

chains <- load_mcmc_chains(location=getwd(),par_tab=example_par_tab,burnin = 1000)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = example_inf_hist,
                individuals=1:5,
                orientation="cross-sectional",
                settings=output$settings)
plot_attack_rates_pointrange( infection_histories = chains$inf_chain,true_ar=all_simulated_data$attack_rates,settings=output$settings)
plot_attack_rates( infection_histories = chains$inf_chain,true_ar=all_simulated_data$attack_rates,settings=output$settings)

################################################
## 2. LONGITUDINAL DATA
################################################
par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")
antigenic_map <- NULL
possible_exposure_times <- 2000:2024

## Vector of antigens that have biomarker measurements (note only one representative antigen per time)
sampled_antigens <- 2024

## Times at which serum samples can be taken
sampling_times <- 2020:2024

## Number of serum samples taken
n_samps <- 5

## Simulate some random attack rates
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.1)

par_tab[par_tab$names %in% c("cr_long","cr_short"),"fixed"] <- 1
par_tab[par_tab$names %in% c("cr_long","cr_short"),"values"] <- 1

par_tab[par_tab$names =="min_measurement","values"] <- 0
par_tab[par_tab$names =="max_measurement","values"] <- 10

par_tab[par_tab$names =="boost_long","values"] <- 2
par_tab[par_tab$names =="boost_short","values"] <- 2

par_tab[par_tab$names =="obs_sd","values"] <- 0.5

## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=100,
                                    possible_exposure_times=possible_exposure_times,
                                    measured_biomarker_ids = sampled_antigens,
                                    sampling_times=sampling_times, nsamps=n_samps,
                                    antigenic_map=antigenic_map,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats=1,
                                    data_type=c(2))

## Pull out the simulated titre data and infection histories
antibody_data <- all_simulated_data$antibody_data
true_inf_hist <- all_simulated_data$infection_histories
plot_antibody_data(antibody_data,possible_exposure_times,1:25,infection_histories = true_inf_hist,study_design="longitudinal") + facet_wrap(~individual)
