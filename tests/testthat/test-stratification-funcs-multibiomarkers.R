## Script to make sure that parameter tables and stratifications are created correctly for different ways of stratifying parameters in serosolver 

## Extend parameter table for each aditional observation type
extend_input_types <- function(par_tab, antigenic_map, n_obs_types){
  antigenic_map_all <- antigenic_map %>% mutate(biomarker_group=1)
  par_tab_all <- par_tab %>% mutate(biomarker_group=1)
  
  if(n_obs_types > 1){
    for(i in 2:n_obs_types){
      par_tab_tmp <- par_tab
      antigenic_map_tmp <- antigenic_map
      par_tab_tmp$biomarker_group <- i
      antigenic_map_tmp$biomarker_group <- i
      par_tab_all <- bind_rows(par_tab_all %>% filter(!(names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"))), par_tab_tmp)
      antigenic_map_all <- bind_rows(antigenic_map_all, antigenic_map_tmp)
    }
  }
  rownames(antigenic_map_all) <- NULL
  rownames(par_tab) <- NULL
  return(list(par_tab=par_tab_all,antigenic_map=antigenic_map_all))
}

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(tidyverse)

devtools::load_all("~/Documents/GitHub/serosolver")
stratification_types <- c("none","single","multiple")
demographic_types <- c("none","fixed","timevarying")
measurement_offsets <- c("none","present","demographics")
biomarker_types <- c("one","multiple")
all_model_combos <- expand_grid(stratification_types,demographic_types)
n_indiv <- 100

############################################################
## 1. No stratification at all

example_par_tab <- read.csv("par_tab.csv")
example_par_tab <- extend_input_types(example_par_tab,example_antigenic_map,2)[[1]]
example_antibody_data <- expand_grid(individual=1:n_indiv,measurement=0,biomarker_id=10,biomarker_group=1,birth=1,repeat_number=1,sample_time=10)
timevarying_demographics <- NULL
demographics <- NULL

par_tab1 <- add_scale_pars(example_par_tab,example_antibody_data) ## Should be unchanged
unique_demographics <- create_demographic_table(example_antibody_data,example_par_tab) ## Should be a table with one row/column labeled all
setup_stratification_table(example_par_tab, unique_demographics) ## Should be a single matrix with one row, 14 columns, all entries 0
get_demographic_groups(example_par_tab,example_antibody_data,NULL,NULL) ## Should match unique_demographics


## 1. Antibody data should just have two extra columns added, demographic_group and population_group, all entries are 1
## 2. the updated demographics/timevarying demographics -- should be NULL here
## 3. Table of unique population groups -- should be NULL here
## 4. list of stratifications for population groups -- should be NULL
## 5. The demographic grouping for each individual
## 6. The demographic grouping for each individual
tmp <- add_stratifying_variables(example_antibody_data, timevarying_demographics, example_par_tab)
unique(tmp$antibody_data$demographic_group)
unique(tmp$antibody_data$population_group)

example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_long","values"] <- 0.75
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_short","values"] <- 1
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "wane_short","values"] <- 0.75

all_simulated_data <- simulate_data(par_tab=example_par_tab, 
                                    n_indiv=n_indiv, ## Number of individuals for the overall simulation
                                    possible_exposure_times=1:10,
                                    measured_biomarker_ids = 10, 
                                    sampling_times=1:10, 
                                    age_min=10,
                                    age_max=10,
                                    nsamps=5,
                                    antigenic_map=NULL, 
                                    attack_rates=simulate_attack_rates(1:10), 
                                    repeats=1, 
                                    missing_data = 0,
                                    data_type=c(1,2), 
                                    demographics=NULL,
                                    measurement_bias = NULL,
                                    verbose=TRUE)
prior_func <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  par_names <- par_tab$names
  f <- function(pars){
    prior_p <- 0
    ## ==================================================
    ## IMPORTANT -- prior function for the model coefficients. This places a standard normal prior on each of the coefficients. 
    ## This is quite strong but shrinks their effects to 0
    prior_p <- sum(dnorm(pars[coef_pars],0,1,log=TRUE))
    prior_p
  }
}
res <- serosolver(example_par_tab, 
                  as.data.frame(all_simulated_data$antibody_data), 
                  demographics=NULL,
                  antigenic_map=NULL,
                  possible_exposure_times = 1:10,
                  prior_func=prior_func,
                  filename="tmp",
                  n_chains=3, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=c(adaptive_iterations=100000, iterations=100000,proposal_ratio=1,thin_inf_hist=10), 
                  verbose=TRUE,
                  data_type=c(1,2),
                  measurement_bias= NULL) 
res$all_diagnostics$p_thetas[[1]]
chains <- load_mcmc_chains(getwd(),example_par_tab,burnin=100000,estimated_only = TRUE)
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_line(aes(x=samp_no,y=value,col=as.factor(chain_no))) +
  facet_wrap(~name,scales="free")
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_density(aes(x=value)) +
  facet_wrap(~name,scales="free")
############################################################
break

############################################################
## 2. Fixed stratification using demographics
example_par_tab <- read.csv("par_tab.csv")
example_par_tab <- extend_input_types(example_par_tab,example_antigenic_map,2)[[1]]

example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_long","values"] <- 0.75
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_short","values"] <- 1
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "wane_short","values"] <- 0.75

demographics <- data.frame(individual=1:n_indiv,urban=sample(c(0,1),n_indiv,replace=TRUE))
example_par_tab[example_par_tab$names == "boost_long","stratification"] <- "urban"
example_antibody_data <- expand_grid(individual=1:n_indiv,measurement=0,biomarker_id=10,biomarker_group=1,birth=1,repeat_number=1,sample_time=10)
example_antibody_data <- align_antibody_demographic_dat(example_antibody_data,demographics)
par_tab1 <- add_scale_pars(example_par_tab,example_antibody_data) ## Should have one extra row with the boost_long urban_1 coef
unique_demographics <- create_demographic_table(example_antibody_data,example_par_tab) ## Should be a table with one column labeled urban, tow entries
setup_stratification_table(example_par_tab, unique_demographics) ## Should be a single matrix with two rows, 14 columns, all entries other than one should be 1
tmp1 <- get_demographic_groups(example_par_tab,example_antibody_data,NULL,NULL) ## Should match unique_demographics


## 1. Antibody data should just have two extra columns added, demographic_group and population_group, all entries are 1
## 2. the updated demographics/timevarying demographics -- should be NULL here
## 3. Table of unique population groups -- should be NULL here
## 4. list of stratifications for population groups -- should be NULL
## 5. The demographic grouping for each individual
## 6. The demographic grouping for each individual
tmp <- add_stratifying_variables(example_antibody_data, timevarying_demographics, example_par_tab,use_demographic_groups=tmp1$use_demographic_groups)
unique(tmp$antibody_data$demographic_group)
unique(tmp$antibody_data$population_group)
par_tab1[par_tab1$names == "boost_long_biomarker_1_coef_urban_1","values"] <- 1
par_tab1[par_tab1$names == "boost_long_biomarker_2_coef_urban_1","values"] <- -1

all_simulated_data <- simulate_data(par_tab=par_tab1, 
                                    n_indiv=n_indiv, ## Number of individuals for the overall simulation
                                    possible_exposure_times=1:10,
                                    measured_biomarker_ids = 10, 
                                    sampling_times=1:10, 
                                    age_min=10,
                                    age_max=10,
                                    nsamps=5,
                                    antigenic_map=NULL, 
                                    attack_rates=simulate_attack_rates(1:10), 
                                    repeats=1, 
                                    missing_data = 0,
                                    data_type=c(1,2), 
                                    demographics=demographics,
                                    measurement_bias = NULL,
                                    verbose=TRUE)

res <- serosolver(example_par_tab, 
                  as.data.frame(all_simulated_data$antibody_data), 
                  demographics=NULL,
                  antigenic_map=NULL,
                  possible_exposure_times = 1:10,
                  prior_func=prior_func,
                  filename="tmp",
                  n_chains=3, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=c(adaptive_iterations=10000, iterations=10000,proposal_ratio=1,thin_inf_hist=10), 
                  verbose=TRUE,
                  data_type=c(1,2),
                  measurement_bias= NULL) 
res$all_diagnostics$p_thetas[[1]]
res$all_diagnostics$p_thetas[[7]]

chains <- load_mcmc_chains(getwd(),par_tab,burnin=10000,estimated_only = FALSE)
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_line(aes(x=samp_no,y=value,col=as.factor(chain_no))) +
  facet_wrap(~name,scales="free")
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_density(aes(x=value)) +
  facet_wrap(~name,scales="free")
############################################################

############################################################
## 3. Timevarying stratification
example_par_tab <- read.csv("par_tab.csv")
example_par_tab <- extend_input_types(example_par_tab,example_antigenic_map,2)[[1]]

example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_long","values"] <- 0.75
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_short","values"] <- 1
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "wane_short","values"] <- 0.75

example_antibody_data <- expand_grid(individual=1:n_indiv,measurement=0,biomarker_id=10,biomarker_group=1,birth=1,repeat_number=1,sample_time=10)
demographics <- data.frame(individual=1:n_indiv,urban=sample(c(0,1),n_indiv,replace=TRUE))
timevarying_demographics <- data.frame(individual=1:n_indiv,urban=sample(c(0,1),n_indiv,replace=TRUE),birth=floor(runif(n_indiv,1,11))) %>% 
  expand_grid(time=1:10) %>% arrange(individual,time) %>% mutate(age_group=if_else(birth >= 5, 0, 1)) 

example_par_tab[example_par_tab$names == "boost_long","stratification"] <- "age_group"

par_tab1 <- add_scale_pars(example_par_tab,example_antibody_data,timevarying_demographics) ## Should have one extra row with the boost_long age_group 1 coef
unique_demographics <- create_demographic_table(timevarying_demographics,example_par_tab) ## Should be a table with one column labeled age group, tow entries
setup_stratification_table(example_par_tab, unique_demographics) ## Should be a single matrix with two rows, 14 columns, all entries other than one should be 1
tmp1 <- get_demographic_groups(example_par_tab,example_antibody_data,timevarying_demographics,NULL) ## Should match unique_demographics and have true for timevarying demographic
tmp <- add_stratifying_variables(example_antibody_data, timevarying_demographics, example_par_tab,use_demographic_groups=tmp1$use_demographic_groups)
unique(tmp$antibody_data$demographic_group)
unique(tmp$antibody_data$population_group)

par_tab1[par_tab1$names == "boost_long_biomarker_1_coef_age_group_1","values"] <- 1
par_tab1[par_tab1$names == "boost_long_biomarker_2_coef_age_group_1","values"] <- -1


all_simulated_data <- simulate_data(par_tab=par_tab1, 
                                    n_indiv=n_indiv, ## Number of individuals for the overall simulation
                                    possible_exposure_times=1:10,
                                    measured_biomarker_ids = 10, 
                                    sampling_times=1:10, 
                                    age_min=10,
                                    age_max=10,
                                    nsamps=5,
                                    antigenic_map=NULL, 
                                    attack_rates=simulate_attack_rates(1:10), 
                                    repeats=1, 
                                    missing_data = 0,
                                    data_type=c(1,2), 
                                    demographics=timevarying_demographics,
                                    measurement_bias = NULL,
                                    verbose=TRUE)

res <- serosolver(example_par_tab, 
                  as.data.frame(all_simulated_data$antibody_data), 
                  demographics=timevarying_demographics,
                  antigenic_map=NULL,
                  possible_exposure_times = 1:10,
                  prior_func=prior_func,
                  filename="tmp",
                  n_chains=3, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=c(adaptive_iterations=10000, iterations=10000,proposal_ratio=1,thin_inf_hist=10), 
                  verbose=TRUE,
                  data_type=c(1,2),
                  measurement_bias= NULL) 
res$all_diagnostics$p_thetas[[1]]
res$all_diagnostics$p_thetas[[7]]
chains <- load_mcmc_chains(getwd(),par_tab,burnin=10000,estimated_only = FALSE)
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_line(aes(x=samp_no,y=value,col=as.factor(chain_no))) +
  facet_wrap(~name,scales="free")
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_density(aes(x=value)) +
  facet_wrap(~name,scales="free")

############################################################

############################################################
## 4. Multiple stratifications, fixed demographics
example_par_tab <- read.csv("par_tab.csv")
example_par_tab <- extend_input_types(example_par_tab,example_antigenic_map,2)[[1]]

example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_long","values"] <- 0.75
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_short","values"] <- 1
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "wane_short","values"] <- 0.75

example_antibody_data <- expand_grid(individual=1:n_indiv,measurement=0,biomarker_id=10,biomarker_group=1,birth=1,repeat_number=1,sample_time=10)

demographics <- data.frame(individual=1:n_indiv,urban=sample(c(0,1),n_indiv,replace=TRUE), other=sample(c(0,1,2),n_indiv,replace=TRUE))
example_par_tab[example_par_tab$names == "boost_long","stratification"] <- "urban, other"
example_par_tab[example_par_tab$names == "boost_short","stratification"] <- "urban, other"

example_antibody_data <- align_antibody_demographic_dat(example_antibody_data,demographics)
par_tab1 <- add_scale_pars(example_par_tab,example_antibody_data) 
unique_demographics <- create_demographic_table(example_antibody_data,example_par_tab) 
setup_stratification_table(example_par_tab, unique_demographics) 
tmp <- get_demographic_groups(example_par_tab,example_antibody_data,NULL,NULL) 
tmp <- add_stratifying_variables(example_antibody_data, NULL, example_par_tab,use_demographic_groups=tmp$use_demographic_groups)
unique(tmp$antibody_data$demographic_group)
unique(tmp$antibody_data$population_group)

all_simulated_data <- simulate_data(par_tab=par_tab1, 
                                    n_indiv=n_indiv, ## Number of individuals for the overall simulation
                                    possible_exposure_times=1:10,
                                    measured_biomarker_ids = 10, 
                                    sampling_times=1:10, 
                                    age_min=10,
                                    age_max=10,
                                    nsamps=5,
                                    antigenic_map=NULL, 
                                    attack_rates=simulate_attack_rates(1:10), 
                                    repeats=1, 
                                    missing_data = 0,
                                    data_type=c(1,2), 
                                    demographics=demographics,
                                    measurement_bias = NULL,
                                    verbose=TRUE)

res <- serosolver(example_par_tab, 
                  as.data.frame(all_simulated_data$antibody_data), 
                  demographics=NULL,
                  antigenic_map=NULL,
                  possible_exposure_times = 1:10,
                  prior_func=prior_func,
                  filename="tmp",
                  n_chains=3, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=c(adaptive_iterations=10000, iterations=10000,proposal_ratio=1,thin_inf_hist=10), 
                  verbose=TRUE,
                  data_type=c(1,2),
                  measurement_bias= NULL) 
res$all_diagnostics$p_thetas[[1]]
res$all_diagnostics$p_thetas[[7]]
chains <- load_mcmc_chains(getwd(),example_par_tab,burnin=10000,estimated_only = TRUE)
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_line(aes(x=samp_no,y=value,col=as.factor(chain_no))) +
  facet_wrap(~name,scales="free")
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_density(aes(x=value)) +
  facet_wrap(~name,scales="free")

############################################################

############################################################
## 5. Multiple stratifications, timevarying demographics
example_par_tab <- read.csv("par_tab.csv")
example_par_tab <- extend_input_types(par_tab,example_antigenic_map,2)[[1]]

example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_long","values"] <- 0.75
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "boost_short","values"] <- 1
example_par_tab[example_par_tab$biomarker_group == 2 & example_par_tab$names == "wane_short","values"] <- 0.75

example_antibody_data <- expand_grid(individual=1:n_indiv,measurement=0,biomarker_id=10,biomarker_group=1,birth=1,repeat_number=1,sample_time=10)

timevarying_demographics <- data.frame(individual=1:n_indiv,urban=sample(c(0,1),n_indiv,replace=TRUE),
                                       other=sample(c(0,1,2),n_indiv,replace=TRUE),
                                       birth=floor(runif(n_indiv,1,11))) %>% 
  expand_grid(time=1:10) %>% arrange(individual,time) %>% mutate(age_group=if_else(birth >= 5, 0, 1)) 

example_par_tab[example_par_tab$names == "boost_long","stratification"] <- "urban, other, age_group"
example_par_tab[example_par_tab$names == "boost_short","stratification"] <- "age_group"
example_par_tab[example_par_tab$names == "infection_model_prior_shape1","stratification"] <- "age_group, urban"
example_par_tab[example_par_tab$names == "infection_model_prior_shape2","stratification"] <- "age_group, urban"

#example_antibody_data <- align_antibody_demographic_dat(example_antibody_data,timevarying_demographics)
par_tab1 <- add_scale_pars(example_par_tab,timevarying_demographics) 
unique_demographics <- create_demographic_table(timevarying_demographics,example_par_tab) 
setup_stratification_table(example_par_tab, unique_demographics) 
tmp <- get_demographic_groups(example_par_tab,example_antibody_data,timevarying_demographics,NULL) 
tmp <- add_stratifying_variables(example_antibody_data, timevarying_demographics, example_par_tab,use_demographic_groups=tmp$use_demographic_groups)
unique(tmp$antibody_data$demographic_group)
unique(tmp$antibody_data$population_group)

par_tab1[par_tab1$names == "boost_long_biomarker_1_coef_urban_1","values"] <- 1
par_tab1[par_tab1$names == "boost_long_biomarker_1_coef_other_1","values"] <- 0.25
par_tab1[par_tab1$names == "boost_long_biomarker_1_coef_other_2","values"] <- 0.75
par_tab1[par_tab1$names == "boost_long_biomarker_1_coef_age_group_1","values"] <- -1
par_tab1[par_tab1$names == "boost_short_biomarker_1_coef_age_group_1","values"] <- -0.5

all_simulated_data <- simulate_data(par_tab=par_tab1, 
                                    n_indiv=n_indiv, ## Number of individuals for the overall simulation
                                    possible_exposure_times=1:10,
                                    measured_biomarker_ids = 10, 
                                    sampling_times=1:10, 
                                    age_min=10,
                                    age_max=10,
                                    nsamps=5,
                                    antigenic_map=NULL, 
                                    attack_rates=simulate_attack_rates(1:10), 
                                    repeats=1, 
                                    missing_data = 0,
                                    data_type=c(1,2), 
                                    demographics=timevarying_demographics,
                                    measurement_bias = NULL,
                                    verbose=TRUE)

res <- serosolver(example_par_tab, 
                  as.data.frame(all_simulated_data$antibody_data), 
                  demographics=timevarying_demographics,
                  antigenic_map=NULL,
                  possible_exposure_times = 1:10,
                  prior_func=prior_func,
                  filename="tmp",
                  n_chains=3, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=c(adaptive_iterations=10000, iterations=10000,proposal_ratio=1,thin_inf_hist=10), 
                  verbose=TRUE,
                  data_type=c(1,2),
                  measurement_bias= NULL) 
res$all_diagnostics$p_thetas[[1]]
res$all_diagnostics$p_thetas[[7]]
chains <- load_mcmc_chains(getwd(),par_tab,burnin=10000,estimated_only = FALSE)
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_line(aes(x=samp_no,y=value,col=as.factor(chain_no))) +
  facet_wrap(~name,scales="free")
ggplot(chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no))) + geom_density(aes(x=value)) +
  facet_wrap(~name,scales="free")

############################################################