## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo=TRUE,
  comment = "#>",
  eval=TRUE
)

## ----message=FALSE, warning=FALSE, r,eval=TRUE---------------------------
# Required to run serosolver
devtools::install_github("seroanalytics/serosolver")
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
registerDoParallel(cl)

## Note that this vignette was generated on a Windows machine,
## and the setup for parallelisation is different on a Linux machine
## for Linux machine:
#  library(doMC)
#  library(doRNG)
#  registerDoMC(cores=5)

## ---- eval=TRUE----------------------------------------------------------
filename <- "case_study_2"
resolution <- 1 ## eg. this would be set to 12 for monthly resolution
sample_year <- 2009

serosolver::describe_proposals()
prior_version <- 2

## ------------------------------------------------------------------------
## Read in data
raw_dat_path <- system.file("extdata", "Fluscape_HI_data.csv", package = "serosolver")
raw_dat <- read.csv(file = raw_dat_path, stringsAsFactors = FALSE)
print(head(raw_dat))

## Add indexing column for each individual
raw_dat$individual <- 1:nrow(raw_dat)

## Convert data to long format
melted_dat <- reshape2::melt(raw_dat, id.vars=c("individual","Age"),stringsAsFactors=FALSE)

## Modify column names to meet serosolver's expectations
colnames(melted_dat) <- c("individual","DOB","virus","titre")
melted_dat$virus <- as.character(melted_dat$virus)

## Extract circulation years for each virus code, which will be used 
## by serosolver as the circulation time
melted_dat$virus <- as.numeric(sapply(melted_dat$virus, function(x) strsplit(x,split = "HI.H3N2.")[[1]][2]))

## Clean and log transform the data
melted_dat <- melted_dat[complete.cases(melted_dat),]
melted_dat[melted_dat$titre == 0,"titre"] <- 5
melted_dat$titre <- log2(melted_dat$titre/5)

## Convert ages to DOB
melted_dat$DOB <- sample_year - melted_dat$DOB

## All samples taken at the same time
melted_dat$samples <- sample_year

## Add column for titre repeats, enumerating for each measurement for the same virus/sample/individual
melted_dat <- plyr::ddply(melted_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x),"group"=1))

## Rename to data expected by serosolver
titre_dat <- melted_dat
print(head(titre_dat))

## ------------------------------------------------------------------------
## Read in raw coordinates
antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
antigenic_coords <- read.csv(file = antigenic_coords_path, stringsAsFactors=FALSE)
print(head(antigenic_coords))

## Convert to form expected by serosolver
antigenic_map <- generate_antigenic_map(antigenic_coords, resolution)
print(head(antigenic_map))

## More flexible version of the above function
virus_key <- c(
    "HK68" = 1968, "EN72" = 1972, "VI75" = 1975, "TX77" = 1977, "BK79" = 1979, "SI87" = 1987, "BE89" = 1989, "BJ89" = 1989,
    "BE92" = 1992, "WU95" = 1995, "SY97" = 1997, "FU02" = 2002, "CA04" = 2004, "WI05" = 2005, "PE06" = 2006
  )
antigenic_coords$Strain <- virus_key[antigenic_coords$Strain]
antigenic_map <- generate_antigenic_map_flexible(antigenic_coords)

## Restrict entries to years of interest. Entries in antigenic_map determine
## the times that individual can be infected ie. the dimensions of the infection
## history matrix.
antigenic_map <- antigenic_map[antigenic_map$inf_years >= 1968 & antigenic_map$inf_years <= sample_year,]
strain_isolation_times <- unique(antigenic_map$inf_years)

## ------------------------------------------------------------------------
par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(file = par_tab_path, stringsAsFactors=FALSE)

## Set parameters for Beta prior on infection histories
beta_pars <- find_beta_prior_mode(0.15,4)
par_tab[par_tab$names == "alpha","values"] <- beta_pars$alpha
par_tab[par_tab$names == "beta","values"] <- beta_pars$beta
## Maximum recordable log titre in these data is 8
par_tab[par_tab$names == "MAX_TITRE","values"] <- 8

## Remove phi parameters, as these are integrated out under prior version 2
par_tab <- par_tab[par_tab$names != "phi",]

## Fix all short term parameters to 0
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"fixed"] <- 1 # mu_short, waning and sigma2 are fixed
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"values"] <- 0 # set these values to 0

## ----Run mcmc 1, eval=TRUE, include=TRUE---------------------------------
## Distinct filename for each chain
no_chains <- 5
filenames <- paste0(filename, "_",1:no_chains)
chain_path <- sub("par_tab_base.csv","",par_tab_path)
chain_path_real <- paste0(chain_path, "cs2_real/")
chain_path_sim <- paste0(chain_path, "cs2_sim/")

## Create the posterior solving function that will be used in the MCMC framework 
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
                  mcmc_pars = c("iterations"=500000,"adaptive_period"=100000,"thin"=100,
                                "thin_hist"=500,"save_block"=1000,
                                "inf_propn"=1, "hist_sample_prob"=1,
                                "hist_switch_prob"=0.8, "year_swap_propn"=1),
                  filename = paste0(chain_path_real,x), 
                  CREATE_POSTERIOR_FUNC = create_posterior_func, 
                  version = prior_version)
}

## ----message=FALSE, warning=FALSE----------------------------------------
## Read in the MCMC chains
## Note that `thin` here is in addition to any thinning done during the fitting
all_chains <- load_mcmc_chains(location=chain_path_real,thin=1,burnin=100000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)

## Alternative, load the included MCMC chains rather than re-running
## load(cs2_chains_real)
## all_chains <- cs2_chains_real

print(summary(all_chains))

