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
devtools::document("~/Documents/GitHub/serosolver")
devtools::load_all("~/Documents/GitHub/serosolver")

## Get all possible infection times
antigenic_map <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/antigenic_maps/antigenicMap_vietnam.csv")
par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")

possible_exposure_times <- seq(2000,2024,by=1)

possible_exposure_times <- c(seq(2000,2020,by=4),seq(2021,2024,by=1))

## Vector of antigens that have biomarker measurements (note only one representative antigen per time)
sampled_antigens <- max(possible_exposure_times)

## Times at which serum samples can be taken
sampling_times <- 2020:2024

## Number of serum samples taken
n_samps <- 5

## Simulate some random attack rates
attack_rates <- runif(length(possible_exposure_times), 0.05, 0.15)


devtools::load_all("~/Documents/GitHub/serosolver")
par_tab[par_tab$names == "wane_long","fixed"] <- 1
par_tab[par_tab$names == "wane_long","values"] <- 0
## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=50,
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

setwd("~/Documents/GitHub/serosolver")
dir.create("test_long")

#setwd("~/Documents/GitHub/serosolver")
#install.packages("~/Documents/GitHub/serosolver",type="source",repos=NULL)

#devtools::document("~/Documents/GitHub/serosolver")
#Rcpp::compileAttributes()
#devtools::load_all("~/Documents/GitHub/serosolver")
par_tab[par_tab$names == "cr_long","fixed"] <- 1
par_tab[par_tab$names == "cr_long","values"] <- 1
par_tab[par_tab$names == "cr_short","fixed"] <- 1
par_tab[par_tab$names == "cr_short","values"] <- 1
par_tab[par_tab$names == "boost_long","values"] <- 0
par_tab[par_tab$names == "boost_long","fixed"] <- 1
antibody_data1 <- antibody_data %>% group_by(individual) %>% mutate(birth=min(sample_time)) %>% as.data.frame()
res <- serosolver(par_tab, antibody_data1, NULL,
                  possible_exposure_times=possible_exposure_times,#2010:max(possible_exposure_times),
                  filename="test_long/readme", prior_version=2,n_chains=1,parallel=FALSE,
                  mcmc_pars=c(adaptive_iterations=20000, iterations=50000,proposal_ratio=2),verbose=TRUE,verbose_dev = TRUE,
                  data_type=2,
                  start_level="mean")
chains <- load_mcmc_chains(location="test_long",par_tab=par_tab,burnin = 20000,unfixed=FALSE)

plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = true_inf_hist[,match(2010:max(possible_exposure_times),possible_exposure_times)],
                individuals=indiv_plots[i]:indiv_plots_bot[i],
                antibody_data=antibody_data1,
                orientation="longitudinal",
                subset_biomarker_ids = NULL,
                expand_to_all_times=FALSE,p_ncol=1,
                settings=res$settings) 

indiv_plots <- seq(1,max(antibody_data$individual)-5,by=4)
indiv_plots_bot <- indiv_plots[2:length(indiv_plots)]
indiv_plots <- indiv_plots[1:(length(indiv_plots)-1)]
pdf("tmp.pdf")
for(i in seq_along(indiv_plots)){
  p1 <- plot_model_fits(chain = chains$theta_chain,
                  infection_histories = chains$inf_chain,
                  known_infection_history = true_inf_hist[,match(2010:max(possible_exposure_times),possible_exposure_times)],
                  antibody_data = antibody_data1,individuals=indiv_plots[i]:indiv_plots_bot[i],
                  antigenic_map=NULL,
                  possible_exposure_times = 2010:max(possible_exposure_times),
                  nsamp=100,
                  par_tab=par_tab,
                  orientation="longitudinal",
                  subset_biomarker_ids = NULL,
                  expand_to_all_times=FALSE,
                  start_level="mean",p_ncol=1) 
  
  p2 <- plot_cumulative_infection_histories(chains$inf_chain,0,indiv_plots[i]:indiv_plots_bot[i],
                                            true_inf_hist[,match(2010:max(possible_exposure_times),possible_exposure_times)],
                                      possible_exposure_times=2010:max(possible_exposure_times)
                                        )[[1]]
  
  print(p1 | p2)
}
dev.off()
