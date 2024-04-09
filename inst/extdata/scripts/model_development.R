## Aims
## -- Add age/sex/group-level antibody kinetics
## -- Add additional exposure types
## -- Allow fixing of some exposure entries
## -- Try different samplers
## -- Fix starting titers
set.seed(1)

library(tidyverse)
devtools::document("~/Documents/GitHub/serosolver")
devtools::load_all("~/Documents/GitHub/serosolver")

## Get all possible infection times
antigenic_map <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/antigenic_maps/antigenicMap_vietnam.csv")
possible_exposure_times <- antigenic_map$inf_times


## Vector of antigens that have biomarker measurements (note only one representative antigen per time)
sampled_antigens <- seq(min(possible_exposure_times), max(possible_exposure_times), by=2)

## Times at which serum samples can be taken
sampling_times <- 2010:2015

## Number of serum samples taken
n_samps <- 5

## Simulate some random attack rates
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)

#antigenic_map <- example_antigenic_map[1:25,]
#antigenic_map$inf_times <- possible_exposure_times
par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")

devtools::load_all("~/Documents/GitHub/serosolver")

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
save(example_par_tab,file="~/Documents/GitHub/serosolver/data/example_par_tab.RData")
save(example_antibody_data,file="~/Documents/GitHub/serosolver/data/example_antibody_data.RData")
save(example_inf_hist,file="~/Documents/GitHub/serosolver/data/example_inf_hist.RData")
