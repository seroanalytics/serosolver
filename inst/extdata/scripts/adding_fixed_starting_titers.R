## Test antibody kinetics models
setwd("~/Documents/GitHub/serosolver")
devtools::document("~/Documents/GitHub/serosolver")
Rcpp::compileAttributes()

devtools::load_all("~/Documents/GitHub/serosolver")

library(tidyverse)
#example_antigenic_map[1:10,c("x_coord","y_coord")] <- 1
antigenic_map <- melt_antigenic_coords(example_antigenic_map[1:25,c("x_coord","y_coord")])

antigenic_map_long <- matrix(create_cross_reactivity_vector(antigenic_map, 0.1),ncol=1)
antigenic_map_short <- matrix(create_cross_reactivity_vector(antigenic_map, 0.03),ncol=1)
exposure_times <- seq(1,25,by=1)
birth <- 3
infection_times <- c(10,15)
infection_indices <- infection_times-1
sample_times <- seq(birth,25,by=1)
biomarker_ids <- exposure_times - 1

start_titres <- runif(length(biomarker_ids), 0,2)
biomarker_ids <- seq(0,24,by=4)

print(start_titres)
y <- antibody_model_individual_wrapper(3,2,1,0.5,0.02,0.05,3,
                                       start_titres,
                                       length(exposure_times),
                                  infection_times,
                                  infection_indices,
                                  biomarker_ids,
                                    sample_times,
                                  antigenic_map_long,
                                  antigenic_map_short)
data.frame(sample_times = rep(sample_times,each=length(biomarker_ids)),biomarker_ids = rep(biomarker_ids,length(sample_times)),titre=y) %>%
  ggplot() + geom_line(aes(x=sample_times,y=titre,col=as.factor(biomarker_ids))) + scale_y_continuous(limits=c(0,8)) + scale_x_continuous(limits=c(0,max(sample_times)))

library(data.table)
library(plyr)

par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")

devtools::load_all("~/Documents/GitHub/serosolver")

example_inf_hist1 <- matrix(0, nrow=50,ncol=48)

antibody_data <- example_antibody_data #%>% filter(repeat_number == 1)
antibody_data <- antibody_data %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
f <- create_posterior_func(par_tab,antibody_data,example_antigenic_map,function_type=4,
                           start_level="max")

antibody_data$y <- f(par_tab$values, example_inf_hist1)

ggplot(antibody_data %>% filter(individual %in% 1:10, biomarker_id %in% seq(1968,2014,by=8))) + 
  geom_point(aes(x=sample_time,y=measurement)) +
  geom_line(aes(x=sample_time,y=y),col="red") +
  facet_grid(individual~biomarker_id)

ggplot(antibody_data %>% filter(individual %in% 1:10)) + 
  geom_point(aes(x=biomarker_id,y=measurement)) +
  geom_line(aes(x=biomarker_id,y=y),col="red") +
  facet_grid(individual~sample_time)

#library(serosolver)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)

serosolver(par_tab, example_antibody_data, example_antigenic_map,
                       filename="readme", prior_version=2,n_chains=1,parallel=FALSE,
                       mcmc_pars=c(adaptive_iterations=10000, iterations=50000),verbose=TRUE,
           start_level="min")

chains <- load_mcmc_chains(location=getwd(),par_tab=example_par_tab,burnin = 100000,unfixed=TRUE)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = example_inf_hist,
                antibody_data = example_antibody_data,individuals=1:4,
                antigenic_map=example_antigenic_map,
                par_tab=example_par_tab,
                orientation="cross-sectional")
