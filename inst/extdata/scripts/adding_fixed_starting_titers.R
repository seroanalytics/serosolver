## Test antibody kinetics models
setwd("~/Documents/GitHub/serosolver")
install.packages("~/Documents/GitHub/serosolver",type="source",repos=NULL)

#devtools::document("~/Documents/GitHub/serosolver")
#Rcpp::compileAttributes()
devtools::load_all("~/Documents/GitHub/serosolver")
#library(serosolver)
rerun<-FALSE
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)

par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")
par_tab[par_tab$names == "wane_long","values"] <- 0.25
par_tab[par_tab$names=="boost_delay","values"] <- 1
devtools::load_all("~/Documents/GitHub/serosolver")
example_inf_hist1 <- matrix(0, nrow=50,ncol=48)

antibody_data <- example_antibody_data %>% #mutate(birth = 2010) %>% 
  arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
f <- create_posterior_func(par_tab,antibody_data,example_antigenic_map,function_type=4,
                           possible_exposure_times = 2010:2015,
                           start_level="none")

antibody_data$y <- f(par_tab$values, example_inf_hist[,43:48])
ggplot(antibody_data %>% filter(individual %in% 1:10, biomarker_id %in% seq(1968,2014,by=8))) + 
  geom_point(aes(x=sample_time,y=measurement)) +
  geom_line(aes(x=sample_time,y=y),col="red") +
  facet_grid(individual~biomarker_id)

ggplot(antibody_data %>% filter(individual %in% 1:10)) + 
  geom_point(aes(x=biomarker_id,y=measurement)) +
  geom_line(aes(x=biomarker_id,y=y),col="red") +
  facet_grid(individual~sample_time)
antibody_data <- antibody_data %>% select(-y)

par_tab[par_tab$names == "wane_long","values"] <- 0
par_tab[par_tab$names == "wane_long","fixed"] <- 1
if(rerun){
#  devtools::load_all("~/Documents/GitHub/serosolver")
  
res <- serosolver(par_tab, antibody_data %>% filter(individual == 1), example_antigenic_map,
                  possible_exposure_times=NULL,
                       filename="chains/readme", prior_version=2,n_chains=3,parallel=FALSE,
                       mcmc_pars=c(adaptive_iterations=10000, iterations=5000,proposal_ratio=2),verbose=TRUE,
           start_level="none")
}
chains <- load_mcmc_chains(location="chains",par_tab=par_tab,burnin = 5000,unfixed=FALSE)

break
devtools::load_all("~/Documents/GitHub/serosolver")
plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = NULL,
                antibody_data = antibody_data,individuals=20:25,
                antigenic_map=example_antigenic_map,
                possible_exposure_times = NULL,
                nsamp=100,
                par_tab=par_tab,
                orientation="longitudinal",
                subset_biomarker_ids = seq(1968,2015,by=8),
                start_level="none") 

plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = example_inf_hist,
                antibody_data = antibody_data,individuals=1:5,
                antigenic_map=example_antigenic_map,
                possible_exposure_times = NULL,
                nsamp=100,
                par_tab=par_tab,
                orientation="cross-sectional",
                subset_biomarker_ids = 1968:2015,
                start_level="none") 
