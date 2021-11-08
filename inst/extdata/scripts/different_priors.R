library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(doMC)
library(foreach)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
setwd("~/Drive/Influenza/serosolver/infection_history_prior_methods_results/")
## How many individuals to simulate?
n_indiv <- 3
buckets <- 1 

par_tab <- read.csv("~/Drive/Influenza/serosolver/age_dependent_boosting/data/parTab_base.csv",stringsAsFactors=FALSE)
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(2,10)
par_tab[par_tab$names == "wane","values"] <- 0.8
par_tab[par_tab$names == "wane","values"] <- par_tab[par_tab$names == "wane","values"]/buckets
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(2,10)
par_tab_lambda <- par_tab
par_tab <- par_tab[par_tab$names != "phi",]
## Add rows for each phi value to be inferred


alpha <- par_tab[par_tab$names == "alpha","values"]
beta <- par_tab[par_tab$names == "beta","values"]

sampling_times <- seq(2009*buckets, 2009*buckets, by=1)
year_min <- 1968*buckets
year_max <- 2009*buckets
age_min <- 75*buckets
age_max <- 75*buckets

## Read in and generate the antigenic map to read strain relationships from
antigenic_map <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, buckets)
fit_dat <- fit_dat[fit_dat$inf_years >= year_min & fit_dat$inf_years <= year_max,]
strain_isolation_times <- unique(fit_dat$inf_years)

tmp <- par_tab_lambda[par_tab_lambda$names == "phi",]
for(i in 1:(length(strain_isolation_times)-1)){
  par_tab_lambda <- rbind(par_tab_lambda, tmp)
}


simInfPars=c("mean"=0.15/buckets,"sd"=0.5,"bigMean"=0.5/buckets/2,"logSD"=1)
attack_rates <- simulate_attack_rates(strain_isolation_times, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])

## Simulate data
dat <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indiv, buckets=buckets,strain_isolation_times=strain_isolation_times,
                     sampling_times=sampling_times, nsamps=2, antigenic_map=fit_dat, titre_sensoring=0, age_min=age_min*buckets,
                     age_max=age_max*buckets,
                     attack_rates=attack_rates)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets

titre_dat <- dat[[1]]
titre_dat <- titre_dat[titre_dat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]
titre_dat <- merge(titre_dat, ages)

## Starting infection histories based on data
start_inf <- setup_infection_histories_new(titre_dat, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
age_mask <- create_age_mask(ages[,2], strain_isolation_times)

## Generate starting locations for MCMC
start_tab <- par_tab
for(i in 1:nrow(start_tab)){
  if(start_tab[i,"fixed"] == 0){
    start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"], 
                                   start_tab[i,"upper_start"])
  }
}

DOBs <- unique(titre_dat[,c("individual","DOB")])[,2]
age_mask <- create_age_mask(DOBs, strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat,strain_isolation_times)
masks <- data.frame(cbind(age_mask, strain_mask))
n_alive_tmp <- sapply(seq(1,length(strain_isolation_times)), function(x)
  nrow(masks[masks$age_mask <=x & masks$strain_mask >= x,]))

## Specify paramters controlling the MCMC procedure
mcmc_pars <- c("save_block"=1000,"thin"=10,"thin_hist"=100,"iterations"=50000,"adaptive_period"=10000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0,"year_swap_propn"=0.5,"swap_propn"=0.5,
               "inf_propn"=1,"hist_sample_prob"=1,"move_size"=2, "hist_opt"=0,"popt"=0.44,"opt_freq"=2000)

f <- create_posterior_func(par_tab,titre_dat,fit_dat,2)
y <- f(par_tab$values, start_inf)



res_gibbs <- run_MCMC(start_tab, titre_dat=titre_dat, 
         antigenic_map = fit_dat,mcmc_pars=mcmc_pars,
         mvr_pars=NULL, filename="chain_no_lik/gibbs/gibbs",PRIOR_FUNC=NULL,
         version=2,  
         start_inf_hist=start_inf,mu_indices=NULL,measurement_random_effects=FALSE,
         measurement_indices=NULL, solve_likelihood=FALSE,
         temp=1)

res_lambda <- run_MCMC(par_tab_lambda, titre_dat=titre_dat, 
                      antigenic_map = fit_dat,mcmc_pars=mcmc_pars,
                      mvr_pars=NULL, filename="chain_no_lik/lambda/lambda",
                      create_posterior_func_fast, PRIOR_FUNC=NULL,
                      version=1,  
                      start_inf_hist=start_inf,mu_indices=NULL,measurement_random_effects=FALSE,
                      measurement_indices=NULL,
                      fast_version=TRUE, solve_likelihood=FALSE,
                      temp=1)


res_beta_binomial <- run_MCMC(start_tab, titre_dat=titre_dat, 
                      antigenic_map = fit_dat,mcmc_pars=mcmc_pars,
                      mvr_pars=NULL, filename="chain_no_lik/beta_binomial/beta_binomial",
                      PRIOR_FUNC=NULL,
                      version=3,  
                      start_inf_hist=start_inf,mu_indices=NULL,measurement_random_effects=FALSE,
                      measurement_indices=NULL,
                      fast_version=TRUE, solve_likelihood=FALSE,
                      temp=1)


res_overall <- run_MCMC(start_tab, titre_dat=titre_dat, 
                              antigenic_map = fit_dat,mcmc_pars=mcmc_pars,
                              mvr_pars=NULL, filename="chain_no_lik/overall/overall",
                              create_posterior_func, PRIOR_FUNC=NULL,
                              version=4,  
                              start_inf_hist=start_inf,mu_indices=NULL,measurement_random_effects=FALSE,
                              measurement_indices=NULL,
                              fast_version=FALSE, solve_likelihood=FALSE,
                              temp=1)


all_chain_gibbs <- load_mcmc_chains(location="chain_no_lik/gibbs",thin=1,burnin=10000,
                                    par_tab=par_tab,unfixed=TRUE,convert_mcmc=TRUE)
inf_chain_gibbs <- all_chain_gibbs$inf_chain
y_gibbs <- generate_cumulative_inf_plots(inf_chain_gibbs,burnin = 0,1,nsamp=100,
                                         strain_isolation_times = strain_isolation_times,
                                         pad_chain=FALSE,number_col = 2)
all_chain_lambda <- load_mcmc_chains(location="chain_no_lik/lambda",thin=1,burnin=10000,
                                     par_tab=par_tab,unfixed=TRUE,convert_mcmc=TRUE)
inf_chain_lambda <- all_chain_lambda$inf_chain
y_lambda <- generate_cumulative_inf_plots(inf_chain_lambda,burnin = 0,1,nsamp=100,
                                          strain_isolation_times = strain_isolation_times,
                                          pad_chain=FALSE,number_col = 2)
all_chains_bb <- load_mcmc_chains(location="chain_no_lik/beta_binomial",thin=1,burnin=10000,
                                  par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)
inf_chain_bb <- all_chains_bb$inf_chain
y_beta_binomial <- generate_cumulative_inf_plots(inf_chain_bb,burnin = 0,1,nsamp=100,
                                                 strain_isolation_times = strain_isolation_times,
                                                 pad_chain=FALSE,number_col = 2)
all_chain_overall <- load_mcmc_chains(location="chain_no_lik/overall",thin=1,burnin=10000,
                                     par_tab=par_tab,unfixed=TRUE,convert_mcmc=TRUE)
inf_chain_overall <- all_chain_overall$inf_chain

y_overall <- generate_cumulative_inf_plots(inf_chain_overall,burnin = 0,1,nsamp=100,
                                          strain_isolation_times = strain_isolation_times,
                                          pad_chain=FALSE,number_col = 2)

p1 <- y_gibbs[[1]] + 
  theme(legend.position="none",strip.background = element_blank(), strip.text.x = element_blank()) + 
  ggtitle("Beta prior on time, integrated out ϕ ") + scale_y_continuous(limits=c(0,40))+ 
  scale_x_continuous(limits=c(1980,(max(titre_dat[titre_dat$individual == 1,"samples"]))))
p2 <- y_lambda[[1]]+  
  theme(legend.position="none",strip.background = element_blank(), strip.text.x = element_blank()) + 
  ggtitle("Beta prior on time, infer ϕ \n(uniform prior on ϕ)") + scale_y_continuous(limits=c(0,35))+ 
  scale_x_continuous(limits=c(1980,(max(titre_dat[titre_dat$individual == 1,"samples"]))))
p3 <- y_beta_binomial[[1]]+ 
  theme(legend.position="none",strip.background = element_blank(), strip.text.x = element_blank()) +  
  ggtitle("Beta prior on lifetime") + scale_y_continuous(limits=c(0,40)) + 
  scale_x_continuous(limits=c(1980,(max(titre_dat[titre_dat$individual == 1,"samples"]))))
p4 <- y_overall[[1]]+ 
  theme(legend.position="none",strip.background = element_blank(), strip.text.x = element_blank()) +  
  ggtitle("Beta prior on overall") + scale_y_continuous(limits=c(0,40))+ 
  scale_x_continuous(limits=c(1980,(max(titre_dat[titre_dat$individual == 1,"samples"]))))

library(cowplot)

p <- plot_grid(p1,p2,p3,p4,ncol=2, align="hv")

svg("cumu_histories_2.svg",width=8,height=6,family="Arial")
print(p)
dev.off()

png("cumu_histories_2.png",width=8,height=6,family="Arial",units="in",res=300)
print(p)
dev.off()

