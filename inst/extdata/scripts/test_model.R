## Test antibody kinetics models
setwd("~/Documents/GitHub/serosolver")
devtools::load_all("~/Documents/GitHub/serosolver")
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)

par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")

example_inf_hist <- matrix(0, nrow=5,ncol=25)
possible_exposure_times <- seq(1,25,by=1)

example_inf_hist <- simulate_infection_histories(runif(length(possible_exposure_times),0,0.15),possible_exposure_times,possible_exposure_times,rep(1,5))[[1]]

antibody_data <- as.data.frame(expand_grid(individual=1:5, biomarker_group=1,sample_time=possible_exposure_times,biomarker_id=25,repeat_number=1,birth=1,measurement=0))
f <- create_posterior_func(par_tab,antibody_data,NULL,function_type=4,
                           possible_exposure_times = possible_exposure_times,
                           start_level="none")

antibody_data$measurement <- f(par_tab$values, example_inf_hist)
ggplot(antibody_data) + 
  geom_line(aes(x=sample_time,y=measurement)) +
  facet_wrap(~individual)
