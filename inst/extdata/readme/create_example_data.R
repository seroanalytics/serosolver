## Script to run the entire README pipeline, simulating the data, saving the example data structures, and testing the MCMC runs
set.seed(1234)
library(serosolver)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
devtools::document("~/Documents/GitHub/serosolver")
devtools::load_all("~/Documents/GitHub/serosolver")
install.packages("~/Documents/GitHub/serosolver",type="source",repos=NULL)

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

par_tab[par_tab$names == "min_measurement","values"] <- 1
par_tab[par_tab$names == "wane_long","values"] <- par_tab[par_tab$names == "wane_short","values"]*0.1 
## Simulate a full serosurvey with these parameters
set.seed(1)
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=50,
                                    possible_exposure_times=possible_exposure_times,
                                    measured_biomarker_ids = sampled_antigens,
                                    sampling_times=sampling_times, nsamps=n_samps,
                                    antigenic_map=antigenic_map,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats=2,
                                    data_type=c(2))

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
                                 filename="chains/readme", n_chains=3,parallel=TRUE,data_type=c(2),
                                 mcmc_pars=c(adaptive_iterations=10000, iterations=20000),verbose=TRUE)
output$all_diagnostics$p_thetas

chains <- load_mcmc_chains(location="chains",par_tab=example_par_tab,burnin = 10000)
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
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)

par_tab[par_tab$names %in% c("cr_long","cr_short"),"fixed"] <- 1
par_tab[par_tab$names %in% c("cr_long","cr_short"),"values"] <- 1

par_tab[par_tab$names =="min_measurement","values"] <- 0
par_tab[par_tab$names =="max_measurement","values"] <- 10

par_tab[par_tab$names =="boost_long","values"] <- 0
par_tab[par_tab$names =="boost_long","fixed"] <- 1
par_tab[par_tab$names =="boost_short","values"] <- 4
par_tab[par_tab$names == "wane_long","values"] <- par_tab[par_tab$names == "wane_short","values"]*0.1 

par_tab[par_tab$names =="obs_sd","values"] <- 1
par_tab[par_tab$names =="boost_delay","values"] <- 1

## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=100,
                                    possible_exposure_times=possible_exposure_times,
                                    measured_biomarker_ids = sampled_antigens,
                                    sampling_times=sampling_times, nsamps=n_samps,
                                    antigenic_map=antigenic_map,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats=2,
                                    data_type=c(2))

## Pull out the simulated titre data and infection histories
antibody_data <- all_simulated_data$antibody_data
true_inf_hist <- all_simulated_data$infection_histories
plot_antibody_data(antibody_data,possible_exposure_times,1:25,infection_histories = true_inf_hist,study_design="longitudinal") + facet_wrap(~individual)

#par_tab[par_tab$names =="boost_long","values"] <- 0
#par_tab[par_tab$names =="boost_long","fixed"] <- 1
#par_tab[par_tab$names =="wane_long","fixed"] <- 0
antibody_data <- antibody_data %>% group_by(individual) %>% mutate(birth=pmax(2020,birth)) %>% ungroup() %>% as.data.frame()
output <- serosolver(par_tab, antibody_data, NULL,possible_exposure_times = 2020:2024,
                     filename="chains/readme", n_chains=3,parallel=TRUE,data_type=c(2),start_level="mean",
                     mcmc_pars=c(adaptive_iterations=10000, iterations=20000),verbose=TRUE)
output$all_diagnostics$p_thetas

chains <- load_mcmc_chains(location="chains",par_tab=example_par_tab,burnin = 10000)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = true_inf_hist[,21:25],
                expand_to_all_times = TRUE,
                start_level="mean",
                settings=output$settings,
                individuals=1:50,
                orientation="longitudinal"
                ) + facet_wrap(~individual)
plot_attack_rates_pointrange( infection_histories = chains$inf_chain,true_ar=all_simulated_data$attack_rates,settings=output$settings)
plot_attack_rates( infection_histories = chains$inf_chain,true_ar=all_simulated_data$attack_rates,settings=output$settings)


################################################
## 3. MULTIPLE BIOMARKER GROUPS
################################################
n_obs_types <- 2
obs_type_dists <- c(1,2) ## Vector of observation models for each observation type -- use 1 for discretized normal (as in previous version) and 2 for continuous, truncated normal
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
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)

par_tab[par_tab$names %in% c("cr_long","cr_short"),"fixed"] <- 1
par_tab[par_tab$names %in% c("cr_long","cr_short"),"values"] <- 1

par_tab[par_tab$names =="min_measurement","values"] <- 0
par_tab[par_tab$names =="max_measurement","values"] <- 10

par_tab[par_tab$names =="boost_long","values"] <- 0
par_tab[par_tab$names =="boost_long","fixed"] <- 1
par_tab[par_tab$names =="boost_short","values"] <- 4
par_tab[par_tab$names == "wane_long","values"] <- par_tab[par_tab$names == "wane_short","values"]*0.1 

par_tab[par_tab$names =="obs_sd","values"] <- 1
par_tab[par_tab$names =="boost_delay","values"] <- 1

## Extend parameter table for each aditional observation type
par_tab$biomarker_group <- 1
par_tab_tmp <- par_tab
if(n_obs_types > 1){
  for(i in 2:n_obs_types){
    par_tab_tmp2 <- par_tab_tmp
    par_tab_tmp2$biomarker_group <- i
    par_tab <- bind_rows(par_tab, par_tab_tmp2 %>% filter(!(names %in% c("alpha","beta"))))
  }
}
par_tab[par_tab$biomarker_group== 1 & par_tab$names %in% c("boost_short"),"values"] <- 2
par_tab[par_tab$biomarker_group== 1 & par_tab$names %in% c("wane_short"),"values"] <- 0.1

all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=100,
                                    possible_exposure_times=possible_exposure_times,
                                    measured_biomarker_ids = sampled_antigens,
                                    sampling_times=sampling_times, nsamps=n_samps,
                                    antigenic_map=antigenic_map,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats=1,
                                    data_type=c(1,2))
antibody_data <- all_simulated_data$antibody_data
true_inf_hist <- all_simulated_data$infection_histories
plot_antibody_data(antibody_data,possible_exposure_times,1:4,infection_histories = true_inf_hist,study_design="longitudinal")
antibody_data <- antibody_data %>% group_by(individual) %>% mutate(birth=pmax(2020,birth)) %>% ungroup() %>% as.data.frame()
output <- serosolver(par_tab, antibody_data, NULL,possible_exposure_times = 2020:2024,
                     filename="chains/readme", n_chains=3,parallel=TRUE,data_type=c(1,2),start_level="mean",
                     mcmc_pars=c(adaptive_iterations=10000, iterations=20000),verbose=TRUE)
output$all_diagnostics$p_thetas
chains <- load_mcmc_chains(location="chains",par_tab=example_par_tab,burnin = 10000)
antibody_data <- antibody_data %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                antibody_data=antibody_data,
                known_infection_history = true_inf_hist[,21:25],
                expand_to_all_times = TRUE,
                nsamp=1000,
                start_level="mean",
                settings=output$settings,
                individuals=1:10,
                orientation="longitudinal",
                subset_biomarker_groups = c(1,2)
)
plot_attack_rates_pointrange( infection_histories = chains$inf_chain,true_ar=all_simulated_data$attack_rates,settings=output$settings)
