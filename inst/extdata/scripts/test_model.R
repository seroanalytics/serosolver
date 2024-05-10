## Test antibody kinetics models
setwd("~/Documents/GitHub/serosolver")
dir.create("~/Documents/local_data/serosolver_testing/test_stratification")
setwd("~/Documents/local_data/serosolver_testing/test_stratification")

devtools::load_all("~/Documents/GitHub/serosolver")
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)

set.seed(1)
par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")
par_tab[par_tab$names %in% c("cr_long","cr_short"),"fixed"] <- 1
par_tab[par_tab$names == "wane_long","values"] <- 0
par_tab[par_tab$names == "wane_long","fixed"] <- 1
par_tab[par_tab$names == "antigenic_seniority","fixed"] <- 1
par_tab[par_tab$names == 'infection_model_prior_shape1',"values"] <- 1/3
par_tab[par_tab$names == 'infection_model_prior_shape2',"values"] <- 1/3
n_indiv <- 25
possible_exposure_times <- seq(1,25,by=1)

example_inf_hist <- simulate_infection_histories(runif(length(possible_exposure_times),0,0.15),possible_exposure_times,possible_exposure_times,rep(1,5))[[1]]

example_inf_hist <- matrix(sample(c(0,1),size=n_indiv*length(possible_exposure_times),prob = c(0.95,0.05),replace=TRUE), nrow=n_indiv,ncol=25)
#example_inf_hist[,3] <- 1
#example_inf_hist[,15] <- 1

antibody_data <- as.data.frame(expand_grid(individual=1:n_indiv, biomarker_group=1,sample_time=seq(1,25,by=1),biomarker_id=25,repeat_number=1,birth=1,measurement=0))
#antibody_data <- antibody_data %>% mutate(demographic_group = if_else(individual %in% 1:3, 1, 2))
antibody_data <- antibody_data %>% mutate(age_group = if_else(individual %in% 1:16, if_else(individual %in% 1:8, 1,2), 3))
antibody_data <- antibody_data %>% mutate(sex = if_else(individual %in% 1:12, 1,2))

#antibody_data$demographic_group <- 1
#antibody_data <- antibody_data %>%mutate(demographic_group = if_else(individual %in% 1:6, if_else(individual %in% 1:3, 1,2), 3))
par_tab$stratification <- NA
par_tab$biomarker_group <- 1
#par_tab[par_tab$names %in% c("boost_long"),"stratification"] <- c("sex")
#par_tab[par_tab$names %in% c("wane_short"),"stratification"] <- c("age_group")
par_tab1 <- add_scale_pars(par_tab,antibody_data)
#par_tab1[16,"values"] <- -1
#par_tab1[18,"values"] <- -1
#par_tab1[16,"values"] <- -1
#par_tab1[17,"values"] <- 1
#par_tab1[19,"values"] <- -1
#par_tab1[20,"values"] <- 1
f <- create_posterior_func(par_tab1,antibody_data,NULL,function_type=3,
                           possible_exposure_times = possible_exposure_times,
                           start_level="none")

antibody_data$measurement <- f(par_tab1$values, example_inf_hist)
ggplot(antibody_data) + 
  geom_line(aes(x=sample_time,y=measurement)) +
  facet_wrap(~individual)
antibody_data <- bind_rows(antibody_data %>% mutate(repeat_number=1),antibody_data %>% sample_frac(0.2) %>% mutate(repeat_number=2))

f <- create_posterior_func(par_tab1,antibody_data,NULL,function_type=1,
                           possible_exposure_times = possible_exposure_times,
                           start_level="none")
f(par_tab1$values, example_inf_hist)

antibody_data <- antibody_data %>% arrange(individual, sample_time, repeat_number)

#par_tab$biomarker_group <- 1
#par_tab$stratification <- NA
#antibody_data$measurement <- antibody_data$measurement + rnorm(nrow(antibody_data),0,0.1)
antibody_data$measurement <- pmin(antibody_data$measurement, 8)
antibody_data$measurement <- pmax(antibody_data$measurement, 0)
#antibody_data$measurement <- floor(antibody_data$measurement)
#Rprof(tmp<-tempfile())
plot_antibody_data(antibody_data,possible_exposure_times,n_indivs=5,study_design = "longitudinal") + facet_wrap(~individual,ncol=1)
res <- serosolver(par_tab, antibody_data , NULL,start_inf_hist = NULL,
                  possible_exposure_times=possible_exposure_times,
                  filename="dev_long_fixed_start", prior_version=2,n_chains=1,parallel=FALSE,
                  mcmc_pars=c(adaptive_iterations=20000, iterations=50000,proposal_ratio=1,thin_inf_hist=100),verbose=TRUE,verbose_dev = TRUE,
                  start_level = "none",
                  data_type=2)
#Rprof(NULL)
#rprofsum <- summaryRprof(tmp)
#rprofsum$by.self %>% arrange(-self.pct) %>%head
#res$settings$par_tab <- res$settings$par_tab %>% filter(par_type != 4)
chains <- load_mcmc_chains(getwd(),par_tab,burnin=000)
plot_mcmc_diagnostics(getwd(),par_tab,0)$p_thetas[[1]]
#plot(coda::as.mcmc(chains$theta_chain))
true_ar <- data.frame(time=possible_exposure_times,AR=colSums(example_inf_hist)/get_n_alive(antibody_data,possible_exposure_times))
antibody_data$population_group <- 1

plot_attack_rates_pointrange(chains$inf_chain,possible_exposure_times=possible_exposure_times,antibody_data=antibody_data,settings = res$settings,prior_pars = list(prior_version="2",infection_model_prior_shape1=2,infection_model_prior_shape2=5),true_ar=true_ar)
plot_model_fits(chains$theta_chain,chains$inf_chain,antibody_data=antibody_data,individuals = 1:25,settings=res$settings,orientation="longitudinal",known_infection_history = example_inf_hist,data_type = 1,
                expand_to_all_times = TRUE) + facet_wrap(~individual,ncol=5)
plot_cumulative_infection_histories(chains$inf_chain,indivs=1:9,possible_exposure_times = possible_exposure_times)[[1]]
res$all_diagnostics$p_thetas[[2]]
