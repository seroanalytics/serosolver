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

filename <- "chains/sim_temp_test1"

## How many individuals to simulate?
n_indiv <- 100
buckets <- 1 ## Set to 1 for annual model. Greater than 1 gives subannual (eg. buckets = 2 is infection period every half year)

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

simInfPars=c("mean"=0.15/buckets,"sd"=0.5,"bigMean"=0.5/buckets/2,"logSD"=1)
##useSIR
#attackRates  <- simulate_ars_buckets(strain_isolation_times, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
#attackRates  <- simulate_attack_rates(strain_isolation_times, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
#attackRates  <- simulate_ars_spline(strain_isolation_times, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"], knots,theta)
attackRates <- simulate_attack_rates(strain_isolation_times, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])

## Simulate data
dat <- simulate_data(par_tab, 1, n_indiv, buckets,strain_isolation_times,
                     sampling_times, 3, antigenic_map=fit_dat, 0, age_min*buckets,age_max*buckets,
                    attack_rates=attackRates)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets

titre_dat <- dat[[1]]
#titre_dat <- titre_dat[titre_dat$virus %in% viruses,]
inf_hist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]
titre_dat <- merge(titre_dat, ages)
titre_dat$group <- 1

## If we want to use or save pre-simulated data
#write.table(titre_dat,"~/net/home/serosolver/data/sim_1000_data.csv",sep=",",row.names=FALSE)
#write.table(inf_hist,"~/net/home/serosolver/data/sim_1000_infHist.csv",sep=",",row.names=FALSE)
#write.table(ages,"~/net/home/serosolver/data/sim_1000_ages.csv",sep=",",row.names=FALSE)
#write.table(AR,"~/net/home/serosolver/data/sim_1000_AR.csv",sep=",",row.names=FALSE)

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
  start_inf <- setup_infection_histories_new_2(titre_dat,strain_isolation_times,5,2)
  y <- f(par_tab$values, start_inf)
  index <- index + 1
}
## Run the MCMC using the inputs generated above
#startTab$values <- par_tab$values
#startInf <- infHist

devtools::load_all()
mcmc_pars <- c("iterations" = 50000, "popt" = 0.44, "popt_hist" = 0.44, "opt_freq" = 2000, "thin" = 1,
               "adaptive_period" = 10000,
               "save_block" = 100, "thin_hist" = 10, "hist_sample_prob" = 0.5, "switch_sample" = 2, "burnin" = 0,
               "inf_propn" = 0.5, "move_size" = 5, "hist_opt" = 1, "swap_propn" = 0.5,
               "hist_switch_prob" = 0, "year_swap_propn" = 1
)


group_ids_vec <- unique(titre_dat[,c("individual","group")])[,"group"]-1
n_groups <- length(unique(group_ids_vec))
n_infections <- sum_infections_by_group(inf_hist,group_ids_vec,n_groups)
n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
print(inf_mat_prior_group_cpp(n_infections, n_alive, 1,1))
print(inf_mat_prior_cpp(inf_hist,n_alive[1,],1,1))

n_groups <- 3
indiv_groups <- NULL
for(i in 1:length(unique(titre_dat$individual))){
  indiv_groups[i] <- sample(1:n_groups, 1)
}

group_table <- data.frame("individual"=1:length(indiv_groups),"group"=indiv_groups)
titre_dat <- merge(titre_dat[,c("individual","titre","samples","virus","run","DOB")],group_table)


group_ids_vec <- unique(titre_dat[,c("individual","group")])[,"group"]-1
n_groups <- length(unique(group_ids_vec))
n_infections <- sum_infections_by_group(inf_hist,group_ids_vec,n_groups)
n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
n_alive_overall <- get_n_alive(titre_dat,strain_isolation_times)
print(inf_mat_prior_group_cpp(n_infections, n_alive, 1,1))
print(inf_mat_prior_cpp(inf_hist,n_alive_overall,1,1))

res <- run_MCMC(par_tab, titre_dat=titre_dat, 
                antigenic_map=fit_dat,mcmc_pars=mcmc_pars,
                mvrPars=NULL, filename="chains/sim_temp_test1",
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
plot(coda::as.mcmc(chain))

titre_dat1 <- titre_dat
titre_dat1$group <- 1
res <- run_MCMC(par_tab, titre_dat=titre_dat1, 
                antigenic_map=fit_dat,mcmc_pars=mcmc_pars,
                mvrPars=NULL, filename="chains/sim_temp_test2",
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
plot(coda::as.mcmc(chain))

## Plot inferred attack rates against true simulated attack rates
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
n_alive <- sapply(1:length(strain_isolation_times), function(x) length(age_mask[age_mask<=x]))
n_alive <- get_n_alive(titre_dat, strain_isolation_times)
inf_prop <- colSums(inf_hist)/n_alive
inf_prop <- data.frame(AR=inf_prop,year=strain_isolation_times)
AR_p <- plot_attack_rates(infChain, titre_dat, strain_isolation_times, n_alive)
AR_p <- AR_p +  geom_point(data=AR,aes(x=year,y=AR)) + 
  geom_point(data=inf_prop,aes(x=year,y=AR),col="purple") + scale_y_continuous(expand=c(0,0),limits=c(0,1))

## Density/trace plots on total number of infections
indivs <- sample(n_indiv, 10)
sampd <- sample(n_indiv,20)
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_strain <- max(infChain$j)
data.table::setkey(infChain, "j","sampno")
n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]

inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~j)


n_indiv <- max(infChain$i)
data.table::setkey(infChain, "i","sampno")
i_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]

inf_chain_p_i <- ggplot(i_inf_chain) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~i)

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
ps <- generate_cumulative_inf_plots(res$history_file, 0, 
                                    sampd, infHist, startInf,strainIsolationTimes, 100,ages,numberCol=4)


## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
plot_infection_histories(chain, infChain, titreDat, sample(1:n_indiv, 10), fit_dat, par_tab,100)

f <- create_post_func(par_tab,titreDat, fit_dat, NULL, 99, ageMask, mu_indices=NULL,measurement_indices=NULL,temp=1000000)
tmpInf <- f(par_tab$values, startInf,1,1,1,1,1,1)
samps <- 50000
tmpInf <- startInf
omg <- matrix(ncol=ncol(startInf)+1,nrow=nrow(startInf)*(samps))
omg <- matrix(ncol=nrow(startInf)+1,nrow=samps)
index <- 1
for(i in 1:samps){
  #print(i)
  tmp <- f(par_tab$values, tmpInf, 1,1,1,1,3)
  tmpInf <- tmp
  tmp <- rowSums(tmp)
  tmp <- c(tmp,i)
  omg[i,] <- tmp
  #omg[index:(index+nrow(startInf)-1),] <- tmp
  #index <- index + nrow(startInf)
}
colnames(omg) <- c(1:nrow(startInf),"sampno")
omg <- as.data.frame(omg)
wow <- reshape2::melt(omg,id.vars="sampno")
ggplot(wow) + geom_histogram(aes(x=value),binwidth=2) + facet_wrap(~variable)

