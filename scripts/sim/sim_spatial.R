######################
## JAMES HAY 13.08.2018 - jameshay218@gmail.com
## This script fits the serosolver antibody kinetics model to simulated HI titre data
## This particular script uses a gibbs proposal step to resample infection histories
## which integrates out the annual force of infection parameters.

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

filename <- "chains/sim_spatial"

## How many individuals to simulate?
n_indiv_per_group <- 50
n_groups <- 5
buckets <- 4 ## Set to 1 for annual model. Greater than 1 gives subannual (eg. buckets = 2 is infection period every half year)

## Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
par_tab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_base.csv",stringsAsFactors=FALSE)
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1,1)
par_tab[par_tab$names == "wane","values"] <- 1
par_tab[par_tab$names == "wane","values"] <- par_tab[par_tab$names == "wane","values"]/buckets
par_tab <- par_tab[par_tab$names != "phi",]


sampling_times <- seq(2010*buckets, 2015*buckets, by=1)
year_min <- 1968*buckets
year_max <- 2015*buckets
age_min <- 6*buckets
age_max <- 75*buckets

## Read in and generate the antigenic map to read strain relationships from
antigenic_map <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, buckets)
fit_dat <- fit_dat[fit_dat$inf_years >= year_min & fit_dat$inf_years <= year_max,]
strain_isolation_times <- unique(fit_dat$inf_years)

## Simulate data
## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets
sim_ar_means <- runif(n_groups,0.1,0.3)

max_indiv <- 0
spatial_attack_rates <- matrix(nrow=n_groups, ncol=length(strain_isolation_times))
titre_dat <- NULL
inf_hist <- NULL
for(i in 1:n_groups){
  ## Group 1
  sim_inf_pars=c("mean"=sim_ar_means[i]/buckets,"sd"=0.5,"bigMean"=0.5/buckets/2,"logSD"=1)
  attack_rates <- simulate_attack_rates(strain_isolation_times, sim_inf_pars["mean"],sim_inf_pars["sd"],TRUE,sim_inf_pars["bigMean"])
  dat <- simulate_data(par_tab, i, n_indiv_per_group, buckets,strain_isolation_times,
                       sampling_times, 2, antigenic_map=fit_dat, 0, age_min*buckets,age_max*buckets,
                      attack_rates=attack_rates)

  titre_dat_tmp <- dat[[1]]
  titre_dat_tmp <- titre_dat_tmp[titre_dat_tmp$virus %in% viruses,]
  inf_hist_tmp <- dat[[2]]
  ages <- dat[[3]]
  AR <- dat[[4]]
  spatial_attack_rates[i,] <- AR[,2]
  titre_dat_tmp <- merge(titre_dat_tmp, ages)
  titre_dat_tmp$individual <- titre_dat_tmp$individual + max_indiv
  titre_dat_tmp$group <- i
  titre_dat <- rbind(titre_dat, titre_dat_tmp)
  inf_hist <- rbind(inf_hist, inf_hist_tmp)
  max_indiv <- max(titre_dat$individual)
}

## If we want to use or save pre-simulated data
#write.table(titre_dat,"~/net/home/serosolver/data/sim_1000_data.csv",sep=",",row.names=FALSE)
#write.table(inf_hist,"~/net/home/serosolver/data/sim_1000_infHist.csv",sep=",",row.names=FALSE)
#write.table(ages,"~/net/home/serosolver/data/sim_1000_ages.csv",sep=",",row.names=FALSE)
#write.table(AR,"~/net/home/serosolver/data/sim_1000_AR.csv",sep=",",row.names=FALSE)

## Starting infection histories based on data
start_inf <- setup_infection_histories_new(titre_dat, unique(fit_dat$inf_years), space=5,titre_cutoff=2, sample_prob=0.1)
age_mask <- create_age_mask(ages[,2], strain_isolation_times)

## Generate starting locations for MCMC
start_tab <- par_tab
for(i in 1:nrow(start_tab)){
  if(start_tab[i,"fixed"] == 0){
    start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"], 
                                   start_tab[i,"upper_start"])
  }
}

## Specify paramters controlling the MCMC procedure

devtools::load_all()
f <- create_posterior_func_fast(par_tab,titre_dat, fit_dat, 1, TRUE)
probs <- f(par_tab$values, inf_hist)
sampled_indivs <- 1:100
n_infs <- rep(1, length(sampled_indivs))
alpha <- 1
beta <- 1
swap_propn <- 0.5
swap_dist <- 3
temp <- 1
proposal <- create_posterior_func_fast(par_tab,titre_dat, fit_dat, function_type=2, TRUE)
proposal(par_tab$values, inf_hist, probs, sampled_indivs, alpha, beta, n_infs, swap_propn, swap_dist, temp)

y <- f(par_tab$values, start_inf)
index <- 1
while(!is.finite(sum(y)) & index < 100){
  print(index)
  start_tab <- par_tab
  for(i in 1:nrow(start_tab)){
    if(start_tab[i,"fixed"] == 0){
      start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"], 
                                     start_tab[i,"upper_start"])
    }
  }
  
  start_inf <- setup_infection_histories_new_2(titre_dat,strain_isolation_times,5,2,sample_prob=0.1)
  y <- f(start_tab$values, start_inf)
  index <- index + 1
}
## Run the MCMC using the inputs generated above
#startTab$values <- par_tab$values
#startInf <- infHist

devtools::load_all()
mcmc_pars <- c("iterations" = 50000, "popt" = 0.44, "popt_hist" = 0.44, "opt_freq" = 2000, "thin" = 1,
               "adaptive_period" = 10000,
               "save_block" = 1000, "thin_hist" = 10, "hist_sample_prob" = 0.5, "switch_sample" = 2, "burnin" = 0,
               "inf_propn" = 0.5, "move_size" = 5, "hist_opt" = 1, "swap_propn" = 0.5,
               "hist_switch_prob" = 0, "year_swap_propn" = 1
)


group_ids_vec <- unique(titre_dat[,c("individual","group")])[,"group"]-1
n_infections <- sum_infections_by_group(inf_hist,group_ids_vec,n_groups)
n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
print(inf_mat_prior_group_cpp(n_infections, n_alive, 1,1))
print(inf_mat_prior_cpp(inf_hist[1:n_indiv_per_group,],n_alive[1,],1,1))


res <- run_MCMC(start_tab, titre_dat=titre_dat, 
                antigenic_map=fit_dat,mcmc_pars=mcmc_pars,
                mvrPars=NULL, filename=filename,
               CREATE_POSTERIOR_FUNC=create_posterior_func_fast, PRIOR_FUNC=NULL,
                version=2,  n_alive=NULL,
                start_inf_hist=start_inf,fast_version = TRUE,
                temp=1)
beepr::beep(4)


#########################
## Processing outputs
#########################
## Density/trace plots
chain <- read.csv(res$chain_file)
chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
#plot(coda::as.mcmc(chain))

## Plot inferred attack rates against true simulated attack rates
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]


titre_dat1 <- titre_dat
titre_dat1$group <- 1
filename1 <- "chains/sim_spatial_2"
res <- run_MCMC(par_tab, titre_dat=titre_dat1, 
                antigenic_map=fit_dat,mcmc_pars=mcmc_pars,
                mvrPars=NULL, filename=filename1,
                CREATE_POSTERIOR_FUNC=create_posterior_func_fast, PRIOR_FUNC=NULL,
                version=2,  n_alive=NULL,
                start_inf_hist=start_inf,fast_version = TRUE,
                temp=1)
beepr::beep(4)

infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
#n_alive <- sapply(1:length(strain_isolation_times), function(x) length(age_mask[age_mask<=x]))
n_alive <- get_n_alive(titre_dat, strain_isolation_times)
inf_prop <- colSums(inf_hist)/n_alive
inf_prop <- data.frame(AR=inf_prop,year=strain_isolation_times)
AR_p <- plot_attack_rates(infChain, titre_dat, strain_isolation_times, pad_chain=FALSE,
                          n_alive,prior_pars=list(prior_version=2,alpha=1,beta=1),plot_den=TRUE)
AR_p <- AR_p + 
  geom_point(data=inf_prop,aes(x=year,y=AR),col="purple") 

