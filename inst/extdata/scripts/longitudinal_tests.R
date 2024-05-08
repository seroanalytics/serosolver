## Aims
## -- Add age/sex/group-level antibody kinetics
## -- Add additional exposure types
## -- Allow fixing of some exposure entries
## -- Try different samplers
## -- Fix starting titers
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(tidyverse)
#devtools::document("~/Documents/GitHub/serosolver")
#devtools::load_all("~/Documents/GitHub/serosolver")
install.packages("~/Documents/GitHub/serosolver", repos = NULL, type="source")
library(serosolver)
## Base parameter table
par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")
par_tab[par_tab$names == "cr_long","fixed"] <- 1
par_tab[par_tab$names == "cr_long","values"] <- 1
par_tab[par_tab$names == "cr_short","fixed"] <- 1
par_tab[par_tab$names == "cr_short","values"] <- 1
par_tab[par_tab$names == "wane_short","values"] <- 0.8
## Annual exposures from 2000 onwards
possible_exposure_times <- seq(2000,2024,by=1)

## Vector of antigens that have biomarker measurements (note only one representative antigen per time)
sampled_antigens <- max(possible_exposure_times)

## Times at which serum samples can be taken -- assume randomly distributed through study period
sampling_times <- possible_exposure_times

## Number of serum samples taken
n_samps <- 10

## Simulate some random attack rates
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)

## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, 
                                    n_indiv=100,
                                    possible_exposure_times=possible_exposure_times,
                                    measured_biomarker_ids = sampled_antigens,
                                    sampling_times=sampling_times, nsamps=n_samps,
                                    antigenic_map=NULL,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats=1,
                                    data_type=c(2))

## Pull out the simulated titre data and infection histories
antibody_data <- all_simulated_data$antibody_data
true_inf_hist <- all_simulated_data$infection_histories
plot_antibody_data(antibody_data,possible_exposure_times,1:25,infection_histories = true_inf_hist,study_design="longitudinal") +
  facet_wrap(~individual)

## Run 1: estimate everything from birth
dir.create("~/Documents/local_data/serosolver_testing/dev_long_base")
setwd("~/Documents/local_data/serosolver_testing/dev_long_base")
res <- serosolver(par_tab, antibody_data, NULL,
                  possible_exposure_times=possible_exposure_times,
                  filename="dev_long_base", prior_version=2,n_chains=3,parallel=TRUE,
                  mcmc_pars=c(adaptive_iterations=20000, iterations=50000,proposal_ratio=2),verbose=TRUE,verbose_dev = TRUE,
                  data_type=2)
res$all_diagnostics$p_thetas[[2]] + geom_vline(data=par_tab[par_tab$fixed == 0,] %>% rename(name=names),aes(xintercept=values))
chains <- load_mcmc_chains(location=getwd(),par_tab=par_tab,burnin = 20000,estimated_only = FALSE)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = true_inf_hist,
                individuals=1:25,
                antibody_data=antibody_data,
                orientation="longitudinal",
                subset_biomarker_ids = NULL,
                expand_to_all_times=TRUE,p_ncol=5,
                start_level = "none",
                settings=res$settings) 
plot_attack_rates_pointrange(chains$inf_chain,antibody_data,settings=res$settings,true_ar = all_simulated_data$attack_rates)

## Run 2: truncate when observations were taken and fix starting titer
antibody_data_fixed_start <- antibody_data %>% filter(sample_time > 2010)
antibody_data_fixed_start$birth <- 2011
possible_exposure_times_fixed_start <- possible_exposure_times[possible_exposure_times > 2010]

par_tab[par_tab$names == "wane_long","values"] <- 0
par_tab[par_tab$names == "wane_long","fixed"] <- 1
true_inf_hist_fixed_start <- true_inf_hist[,match(possible_exposure_times_fixed_start,possible_exposure_times)]
true_ar_fixed_start <- data.frame(time=possible_exposure_times_fixed_start,AR=colSums(true_inf_hist_fixed_start)/get_n_alive(antibody_data_fixed_start,possible_exposure_times_fixed_start))
dir.create("~/Documents/local_data/serosolver_testing/dev_long_fixed_start")
setwd("~/Documents/local_data/serosolver_testing/dev_long_fixed_start")
res <- serosolver(par_tab, antibody_data_fixed_start, NULL,
                  possible_exposure_times=possible_exposure_times_fixed_start,
                  filename="dev_long_fixed_start", prior_version=2,n_chains=3,parallel=TRUE,
                  mcmc_pars=c(adaptive_iterations=20000, iterations=50000,proposal_ratio=2),verbose=TRUE,verbose_dev = TRUE,
                  start_level = "min",
                  data_type=2)
res$all_diagnostics$p_thetas[[2]] + geom_vline(data=par_tab[par_tab$fixed == 0,] %>% rename(name=names),aes(xintercept=values))
chains <- load_mcmc_chains(location=getwd(),par_tab=par_tab,burnin = 20000,estimated_only = FALSE)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = true_inf_hist_fixed_start,
                individuals=1:25,
                antibody_data=antibody_data_fixed_start,
                orientation="longitudinal",
                subset_biomarker_ids = NULL,
                expand_to_all_times=TRUE,p_ncol=5,
                start_level = "min",
                settings=res$settings) 
plot_attack_rates_pointrange(chains$inf_chain,antibody_data_fixed_start,settings=res$settings,true_ar = true_ar_fixed_start)
