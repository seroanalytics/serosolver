code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

library(plyr)
library(plotrix)

## Specify paramters controlling the MCMC procedure
resolution <- 12

#use the real SA data to inform the sampling times 
#this file is created in formatting_data
titre_dat_real <- read.csv("data/south_africatitre_dat_full.csv")

## Antigenic map for cross reactivity parameters
antigenic_map <- read.csv(file.path(code.dir,"data/fonville_map_approx.csv"),stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, resolution)
resolution <- 12
n_years <- 3
## All possible circulation times
fit_dat <- fit_dat[1 : (n_years*12),]
fit_dat[2:nrow(fit_dat),c('x_coord','y_coord')] <- fit_dat[1,c('x_coord','y_coord')]
dim(fit_dat)[1] / resolution
max_sample_time <- max(titre_dat_real$samples)
# propose infection times 4 years prior to the sampling time
fit_dat$inf_years <- seq((max_sample_time - n_years*resolution) + 1, max_sample_time, by = 1)
strain_isolation_times <- unique(fit_dat$inf_years)


#####unfixed and with multivariate normal distribution
###no antigenic senority
par_tab <- read.csv("inputs/par_tab_inf_data.csv")
par_tab$values[which(par_tab$names == "wane")] <- par_tab$values[which(par_tab$names == "wane")] / resolution

## Fix all long term parameters to 0
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 #mu, tau, sigma1, and sigma2 are fixed
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 # set these values to 0

tmp <- par_tab[par_tab$names == "phi",]
for(i in 1:(length(strain_isolation_times)-1)){
  par_tab <- rbind(par_tab, tmp)
}

x <- seq(0, 0.2, l = 6)
attack_rates <- rep(c(x, rev(x)), n_years)

sampling_times <- unique(titre_dat_real$samples)
n_indiv <- length(unique(titre_dat_real$individual))
set.seed(54)
all_simulated_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indiv,
              strain_isolation_times=strain_isolation_times,
              sampling_times=sampling_times, nsamps=2,
              antigenic_map=fit_dat,
              age_min=10,age_max=75,
              attack_rates=attack_rates, repeats = 2)

titre_dat <- all_simulated_data$data
titre_dat <- titre_dat[which(titre_dat$virus == unique(titre_dat_real$virus)), ]
ages <- all_simulated_data$ages
inf_hist <- all_simulated_data$infection_histories
titre_dat <- merge(titre_dat, ages)
head(titre_dat)

n_alive <- get_n_alive(titre_dat, strain_isolation_times)
any(colSums(inf_hist)/n_alive >1)


start_inf_hist <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indiv,
                                    strain_isolation_times=strain_isolation_times,
                                    sampling_times=sampling_times, nsamps=2,
                                    antigenic_map=fit_dat,
                                    age_min=10,age_max=75,
                                    attack_rates=attack_rates, repeats = 2)$infection_histories


#simulate PCR data
delta <- par_tab$values[which(par_tab$names == "delta")]
inf_data <- as.matrix(inf_hist)
inf_data <- sapply(1:dim(inf_hist)[2],function(x) {
  inf_vec <- inf_data[,x]
  inf_vec[which(inf_vec == 1)]<-rbinom(length(which(inf_vec == 1)), 1, delta)
  inf_data[,x] <- inf_vec
})

mcmc_pars <- c("iterations"=200000,"adaptive_period"=100000, "burnin"=50000,
               "hist_sample_prob"=0.5,"thin"=5,"thin_hist"=20,"swap_propn"=0.5,
               "hist_switch_prob"=0.2,"year_swap_propn"=0.5,
               "switch_sample" = 2)

run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
                antigenic_map = fit_dat, mcmc_pars =mcmc_pars,
                mvr_pars = NULL, inf_dat = inf_data, start_inf_hist = NULL, filename="chains/sim_fix_delta_full_inf_data",
                CREATE_POSTERIOR_FUNC=create_posterior_func_fast, CREATE_PRIOR_FUNC=NULL,
                version=1,
                fast_version=TRUE)


par_tab$fixed[par_tab$names == "delta"] <- 1
inf_data_partial <- inf_data
inf_data_partial[ , c(1:14, 26:36)] <- NA
run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
         antigenic_map = fit_dat, mcmc_pars =mcmc_pars,
         mvr_pars = NULL, inf_dat = inf_data_partial, start_inf_hist = NULL, filename="chains/sim_fix_delta_partial_inf_data",
         CREATE_POSTERIOR_FUNC=create_posterior_func_fast, CREATE_PRIOR_FUNC=NULL,
         version=1,
         fast_version=TRUE)


filenames <- c("sim_fix_delta_full_inf_data",
               "sim_fix_delta_partial_inf_data")


for(filename in filenames){
  chain <- read.csv(paste("chains/",filename,"_chain.csv", sep = ""))
  chain<- chain[16929: 26929 ,]
  chain<- chain[chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]

  inf_chain  <- data.table::fread(paste("chains/",filename,"_infection_histories.csv", sep = ""))
  inf_chain<- inf_chain[16929: 26929 ,]
  inf_chain  <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]

  plot_attack_rates(inf_chain, titre_dat, strain_isolation_times)
  
  par(mfrow =c (1, 2))
  plot(chain$likelihood, type = "l")
  plot(chain$delta, type ="l")
}