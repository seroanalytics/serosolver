## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo=TRUE,
  comment = "#>",
  eval=TRUE
)

## ----message=FALSE, warning=FALSE, r,eval=TRUE---------------------------
# Required to run serosolver
#devtools::install_github("seroanalytics/serosolver")
library(serosolver)
library(plyr)
library(data.table)

## Required for this analysis
library(reshape2)
library(foreach)
library(doParallel)
library(bayesplot)
library(coda)
library(ggplot2)
library(viridis)

# set up cluster
set.seed(1234)
cl <- makeCluster(5)

## Note that this vignette was generated on a Windows machine,
## and the setup for parallelisation is different on a Linux or Mac:

if(Sys.info()[["sysname"]]=="Darwin" | Sys.info()[["sysname"]]=="Linux"){
  library(doMC)
  library(doRNG)
  registerDoMC(cores=5)
}else{
  registerDoParallel(cl)
}


## ---- eval=TRUE----------------------------------------------------------
filename <- "case_study_1"
resolution <- 4 ## set to 4 for quarterly resolution
sample_years <- 2009:2012

serosolver::describe_proposals()
prior_version <- 2

## ------------------------------------------------------------------------
## Read in titre data
# unvaccinated
input_dat_path <- system.file("extdata", "HKdata_h1n1_unvac.csv", package = "serosolver")
input_dat <- read.csv(file = input_dat_path, header = TRUE)
# vaccinated
# input_dat_path2 <- system.file("extdata", "HKdata_h1n1_vac.csv", package = "serosolver")
# input_dat_vac <- read.csv(file = input_dat_path2, header = TRUE)

indivs <- unique(input_dat$individual) #all individuals

# Subset data for indivs
titre_dat <- input_dat[input_dat$individual %in% indivs,c("individual","virus","titre","samples","DOB")]
titre_dat$individual <- match(titre_dat$individual, indivs)

titre_dat <- unique(titre_dat)
titre_dat <- plyr::ddply(titre_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x),"group"=1))
print(head(titre_dat))

## ------------------------------------------------------------------------
## Read in raw coordinates
antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)

## Convert to form expected by serosolver
antigenic_map <- generate_antigenic_map(antigenic_coords, resolution)

## Restrict entries to years of interest. Entries in antigenic_map determine
## the times that individual can be infected ie. the dimensions of the infection
## history matrix.
antigenic_map <- antigenic_map[antigenic_map$inf_years>=(sample_years[1]*resolution+1) & antigenic_map$inf_years<=sample_years[4]*resolution,]

## Change all coordinates to the same as the first pair of coordinates (because we assume the virus is the same at every time point)
antigenic_map[2:nrow(antigenic_map),c('x_coord','y_coord')] <- antigenic_map[1,c('x_coord','y_coord')]
print(head(antigenic_map))

strain_isolation_times <- unique(antigenic_map$inf_years)

## ------------------------------------------------------------------------
par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(par_tab_path, stringsAsFactors=FALSE)

## Set parameters for beta and alpha to 1
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1,1)
## Maximum recordable log titre in these data is 9
par_tab[par_tab$names == "MAX_TITRE","values"] <- 9

## Remove phi parameters, as these are integrated out under prior version 2
par_tab <- par_tab[par_tab$names != "phi",]

## Fix all long term parameters to 0
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 #mu, tau, sigma1, and sigma2 are fixed
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 # set these values to 0

## ----Run mcmc 1, eval=TRUE, include=TRUE---------------------------------
## Distinct filename for each chain
no_chains <- 5
filenames <- paste0(filename, "_",1:no_chains)
chain_path <- sub("par_tab_base.csv","",par_tab_path)
chain_path_real <- paste0(chain_path, "cs1_real/")
chain_path_sim <- paste0(chain_path, "cs1_sim/")

## Create the posterior solving function that will be used in the MCMC framework 
par_tab[par_tab$names == "mu_short","lower_bound"] <- 1
model_func <- create_posterior_func(par_tab=par_tab,
                            titre_dat=titre_dat,
                            antigenic_map=antigenic_map,
                            version=prior_version) # function in posteriors.R
  
## Generate results in parallel
res <- foreach(x = filenames, .packages = c('serosolver','data.table','plyr')) %dopar% {
  ## Not all random starting conditions return finite likelihood, so for each chain generate random
  ## conditions until we get one with a finite likelihood
  start_prob <- -Inf
  while(!is.finite(start_prob)){
    ## Generating starting antibody kinetics parameters
    start_tab <- generate_start_tab(par_tab)
    
    ## Generate starting infection history
    start_inf <- setup_infection_histories_new_2(titre_dat, strain_isolation_times, space=3,titre_cutoff=4)
    start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
  }
  
  res <- run_MCMC(par_tab = start_tab, 
                  titre_dat = titre_dat,
                  antigenic_map = antigenic_map,
                  start_inf_hist = start_inf, 
                  mcmc_pars = c("iterations"=500000,"adaptive_period"=100000,
                                "thin"=500,"thin_hist"=5000,
                                "save_block"=1000,"inf_propn"=1, 
                                "hist_sample_prob"=1,"hist_switch_prob"=0.8,
                                "year_swap_propn"=1),
                  filename = paste0(chain_path_real,x), 
                  CREATE_POSTERIOR_FUNC = create_posterior_func, 
                  version = prior_version)
}

## ----message=FALSE, warning=FALSE----------------------------------------
## Read in the MCMC chains
# Note that `thin` here is in addition to any thinning done during the fitting
# Chain length values in load function need to be consistent with MCMC run
all_chains <- load_mcmc_chains(location=chain_path_real,thin=1,burnin=100000,
                             par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)
## Alternative, load the included MCMC chains rather than re-running
## load(cs1_chains_real)
## all_chains <- cs1_chains_real

print(summary(all_chains))

## ----fig.height=8, fig.width=11------------------------------------------
## Get the MCMC chains as a list
list_chains <- all_chains$theta_list_chains
## Look at diagnostics for the free parameters
list_chains1 <- lapply(list_chains, function(x) x[,c("mu_short","wane","error",
                                                     "total_infections",
                                                     "lnlike","prior_prob")])

## Gelman-Rubin diagnostics to assess between-chain convergence for each parameter
print(gelman.diag(as.mcmc.list(list_chains1)))
gelman.plot(as.mcmc.list(list_chains1))

## Effective sample size for each parameter
print(effectiveSize(as.mcmc.list(list_chains1)))

## Posterior estimates for each parameter
print(summary(as.mcmc.list(list_chains1)))


## Plot the MCMC trace using the `bayesplot` package
color_scheme_set("viridis")
p_theta_trace <- mcmc_trace(list_chains1)
print(p_theta_trace)

## ----fig.height=8, fig.width=11------------------------------------------
## Extract infection history chain
inf_chain <- all_chains$inf_chain

## Look at inferred attack rates
p_ar <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times, pad_chain=FALSE,
                          plot_den = TRUE,prior_pars=list(prior_version=prior_version, 
                                                          alpha=par_tab[par_tab$names=="alpha","values"],
                                                          beta=par_tab[par_tab$names=="beta","values"])) 
print(p_ar)

## Calculate convergence diagnostics and summary statistics on infection histories
## Important to scale all infection estimates by number alive from titre_dat
n_alive <- get_n_alive_group(titre_dat, strain_isolation_times,melt=TRUE)

## This function generates a number of MCMC outputs
ps_infhist <- plot_posteriors_infhist(inf_chain=inf_chain, 
                                      years=strain_isolation_times, 
                                      samples = 100,  # Needs to be smaller than length of sampled chain 
                                      n_alive=n_alive)


## Posterior mean, median, 95% credible intervals and effective sample size
## on per time attack rates
print(head(ps_infhist[["estimates"]]$by_year))

## Posterior mean, median, 95% credible intervals and effective sample size
## on per individual total number of infections
print(head(ps_infhist[["estimates"]]$by_indiv))

## Check for agreement between inferred cumulative infection histories 
## for some individuals
p_indiv_inf_hists <- generate_cumulative_inf_plots(inf_chain, indivs=1:9, pad_chain=FALSE,
                                                   nsamp = 100, # Needs to be smaller than length of sampled chain 
                                                   strain_isolation_times = strain_isolation_times,
                                                   number_col=3)
print(p_indiv_inf_hists[[1]])

## Posterior probability that infections occured at given times per individual
print(p_indiv_inf_hists[[2]])

## ----fig.height=8, fig.width=11------------------------------------------
## get_titre_predictions expects only a single MCMC chain, so
## subset for only one chain
chain <- as.data.frame(all_chains$theta_chain)
chain1 <- chain[chain$chain_no == 1,]
inf_chain1 <- inf_chain[inf_chain$chain_no == 1,]

titre_preds <- get_titre_predictions(chain = chain1, 
                                     infection_histories = inf_chain1, 
                                     titre_dat = titre_dat, 
                                     individuals = unique(titre_dat$individual),
                                     nsamp = 100, # Needs to be smaller than length of sampled chain 
                                     antigenic_map = antigenic_map, 
                                     par_tab = par_tab,expand_titredat=FALSE)
to_use <- titre_preds$predictions
print(head(to_use))

## Using ggplot
titre_pred_p <- ggplot(to_use[to_use$individual %in% c(1:8,13),])+
  geom_ribbon(aes(x=samples,ymin=lower, ymax=upper),fill="gray90")+
  geom_ribbon(aes(x=samples,ymin=lower_50, ymax=upper_50),fill="gray70")+
  geom_line(aes(x=samples, y=median))+
  geom_point(aes(x=samples, y=titre))+
  coord_cartesian(ylim=c(0,8))+
  ylab("log titre") +
  xlab("Time of virus circulation") +
  theme_classic() +
  facet_wrap(~individual)
titre_pred_p

## ----message=FALSE, warning=FALSE----------------------------------------
## Read in MCMC chains from fitting
all_chains <- load_mcmc_chains(location=chain_path_real,thin=1,burnin=100000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=FALSE)

## Alternative, load the included MCMC chains rather than re-running
## load(cs1_chains_real_b)
## all_chains <- cs1_chains_real_b

## Find samples that were in both theta and inf hist chains
chain <- all_chains$theta_chain
inf_chain <- all_chains$inf_chain
intersect_samps <- intersect(unique(inf_chain$sampno), unique(chain$sampno))
chain <- chain[chain$sampno %in% intersect_samps,]

## Find the parameter values that gave the highest posterior probability
which_mle <- chain[which.max(chain$lnlike),c("sampno","chain_no")]
mle_theta_pars <- chain[chain$sampno == which_mle$sampno & chain$chain_no == which_mle$chain_no,]

## Store total infections to compare later
mle_total_infs <- mle_theta_pars[,"total_infections"]
mle_theta_pars <- mle_theta_pars[,par_tab$names]
mle_inf_hist <- inf_chain[inf_chain$sampno == which_mle$sampno & inf_chain$chain_no == which_mle$chain_no,]

## Generate full infection history matrix using provided function
mle_inf_hist <- expand_summary_inf_chain(mle_inf_hist[,c("sampno","j","i","x")])
## Find number of infections per year from this infection history
no_infs <- colSums(mle_inf_hist[,3:ncol(mle_inf_hist)])

## If missing time points in simulated attack rates
if(length(no_infs) < length(strain_isolation_times)){
  diff_lengths <- length(strain_isolation_times) - length(no_infs)
  no_infs <- c(no_infs, rep(0, diff_lengths))
}

## Find attack rate per year
n_alive <- get_n_alive(titre_dat, strain_isolation_times)
attack_rates <- no_infs/n_alive

## ----fig.height=8, fig.width=11------------------------------------------
set.seed(1234)

sim_par_tab <- par_tab
sim_par_tab$values <- as.numeric(mle_theta_pars)
sim_par_tab[sim_par_tab$names %in% c("alpha","beta"),"values"] <- c(1/3,1/3)
sim_par_tab[sim_par_tab$names == "wane","values"] <- 1
sim_par_tab[sim_par_tab$names == "wane","values"] <- par_tab[par_tab$names == "wane","values"]/resolution
sim_par_tab[sim_par_tab$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 
sim_par_tab[sim_par_tab$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 
sim_par_tab[sim_par_tab$names == "MAX_TITRE","values"] <- 9 

sampling_times <- seq(2009*resolution + 1, 2012*resolution, by=1)

age_min <- 6*resolution
age_max <- 6*resolution
n_indiv <- length(unique(titre_dat$individual))

dat <- simulate_data(par_tab = sim_par_tab, 
                     n_indiv = n_indiv,
                     buckets = resolution,
                     strain_isolation_times = strain_isolation_times,
                     sampling_times = sampling_times, 
                     nsamps = 4, 
                     antigenic_map = antigenic_map, 
                     age_min = age_min,
                     age_max = age_max,
                     attack_rates=attack_rates,
                     repeats = 1)



## Inspect simulated antibody titre data and infection histories
sim_titre_dat <- dat[["data"]]
sim_infection_histories <- dat[["infection_histories"]]

## Store total infections to compare later
actual_total_infections <- sum(sim_infection_histories)

plot_data(sim_titre_dat, sim_infection_histories, strain_isolation_times,n_samps = 5)

## Use titres only against same viruses tested in real data
viruses <- unique(titre_dat$virus)
sim_titre_dat <- sim_titre_dat[sim_titre_dat$virus %in% viruses, ]
sim_ages <- dat[["ages"]]
sim_titre_dat <- merge(sim_titre_dat, sim_ages)
sim_ar <- dat[["attack_rates"]]

## ----Run mcmc sim, eval=FALSE, include=TRUE------------------------------
#  filename <- "case_study_1_sim"
#  
#  ## Distinct filename for each chain
#  no_chains <- 5
#  filenames <- paste0(filename, "_",1:no_chains)
#  
#  ## Create the posterior solving function that will be used in the MCMC framework
#  model_func <- create_posterior_func(par_tab=sim_par_tab,
#                                  titre_dat=sim_titre_dat,
#                                  antigenic_map=antigenic_map,
#                                  version=2) # function in posteriors.R
#  
#  ## Generate results in parallel
#  res <- foreach(x = filenames, .packages = c('serosolver','data.table','plyr')) %dopar% {
#    ## Not all random starting conditions return finite likelihood, so for each chain generate random
#    ## conditions until we get one with a finite likelihood
#    start_prob <- -Inf
#    while(!is.finite(start_prob)){
#      ## Generate starting values for theta
#      start_tab <- generate_start_tab(par_tab)
#      ## Generate starting infection history
#      start_inf <- setup_infection_histories_new_2(sim_titre_dat, strain_isolation_times, space=3,titre_cutoff=4)
#      start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
#    }
#  
#    res <- run_MCMC(par_tab = start_tab,
#                    titre_dat = sim_titre_dat,
#                    antigenic_map = antigenic_map,
#                    start_inf_hist = start_inf,
#                    mcmc_pars = c("iterations"=500000,"adaptive_period"=100000,"thin"=500,
#                                  "thin_hist"=5000,"save_block"=1000,
#                                  "hist_switch_prob"=0.8, "hist_sample_prob"=1),
#                    filename = paste0(chain_path_sim,x),
#                    CREATE_POSTERIOR_FUNC = create_posterior_func,
#                    version = prior_version)
#  }

