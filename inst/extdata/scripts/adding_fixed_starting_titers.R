## Test antibody kinetics models
setwd("~/Documents/GitHub/serosolver")
devtools::document("~/Documents/GitHub/serosolver")
Rcpp::compileAttributes()
devtools::load_all("~/Documents/GitHub/serosolver")

library(data.table)
library(plyr)

par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")
par_tab[par_tab$names == "wane_long","values"] <- 0.1
devtools::load_all("~/Documents/GitHub/serosolver")
example_inf_hist1 <- matrix(0, nrow=50,ncol=48)

antibody_data <- example_antibody_data #%>% filter(repeat_number == 1)
antibody_data <- antibody_data %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
f <- create_posterior_func(par_tab,antibody_data %>% mutate(birth = 2010),example_antigenic_map,function_type=4,
                           start_level="mean")
antibody_data %>% filter(individual == 1,biomarker_id == 1992) %>% dplyr::summarize(min_y=min(measurement))
antibody_data %>% filter(individual == 1,biomarker_id == 1992) %>% dplyr::summarize(min_y=min(y))

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
           start_level="none")

chains <- load_mcmc_chains(location=getwd(),par_tab=par_tab,burnin = 10000,unfixed=TRUE)
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = NULL,
                antibody_data = example_antibody_data,individuals=1:5,
                antigenic_map=example_antigenic_map,
                nsamp=100,
                par_tab=par_tab,
                orientation="cross-sectional")
