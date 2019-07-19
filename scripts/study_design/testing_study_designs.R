code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

## Specify paramters controlling the MCMC procedure
mcmc_pars <- c("iterations"=500000,"adaptive_period"=500000, "burnin"=50000,
               "hist_sample_prob"=0.5,"thin"=200,"thin_hist"=200,
               "swap_propn"=0.5,"hist_switch_prob"=0.2,"year_swap_propn"=0.5)


## Antigenic map for cross reactivity parameters
antigenic_map <- read.csv(file.path(code.dir,"data/fonville_map_approx.csv"),
                         stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, buckets=1)

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968,]
strainIsolationTimes <- unique(fit_dat$inf_years)

## Read in parameter table
par_tab <- read.csv(file.path(code.dir,"inputs/par_tab_study_design.csv"))


library(doParallel)
library(foreach)

#create cluster with desired number of cores
cl<-makeCluster(4)

#register cluster
registerDoParallel(cl)


###mcmc chains
for(f in 1:3){
  foreach(i=1:3)%do%{
    for(type in c("L","CS")){
      filename <- paste0(i,"study",f,"_", type, sep="")
      
      titre_dat <- read.csv(paste("data/study_design/",filename,"dat.csv",sep = ""))
      ages <- read.csv(paste("data/study_design/",filename,"ages.csv",sep = ""))
      
      res <- run_MCMC(par_tab = par_tab, titre_dat = merge(titre_dat,ages, by = "individual"), 
                      antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                      mvr_pars = NULL, start_inf_hist = NULL, filename=paste("chains/", filename,sep=""),
                      CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = NULL,
                      version = 2,  
                      fast_version = TRUE)
    }
  }
}


stopCluster(cl)