code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

library(data.table)
library(ggplot2)

# read in data files
titre_dat <- read.csv("data/south_africa/titre_dat_full.csv")

infection_data_2016 <- read.csv("data/south_africa/infection_data_2016.csv",
                                check.names = F)

#obtain infection data for the indiividuals no's
#which are in titre_dat
ids <- unique(titre_dat$individual)
infection_data_2016 <- infection_data_2016[which(infection_data_2016$ind_id %in% ids),]
infection_data_2016$ind_id <- factor(infection_data_2016$ind_id, levels = levels(titre_dat$individual))

#then add NA's to the sampling times without data
PCR_dates <- as.numeric(colnames(infection_data_2016)[-1])

## Antigenic map for cross reactivity parameters
antigenic_map <- read.csv(file.path(code.dir,"data/fonville_map_approx.csv"),stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, resolution)
resolution <- 12
## All possible circulation times
fit_dat <- fit_dat[1 : (3*12),]
fit_dat[2:nrow(fit_dat),c('x_coord','y_coord')] <- fit_dat[1,c('x_coord','y_coord')]
dim(fit_dat)[1] / resolution
max_sample_time <- max(titre_dat$samples)
# propose infection times 4 years prior to the sampling time
fit_dat$inf_years <- seq((max_sample_time - 3*resolution) + 1, max_sample_time, by = 1)
strain_isolation_times <- unique(fit_dat$inf_years)

# set up full size infection history matrix then fill in the PCR data
# this is done here because it depends on strain isolation times
inf_data <- matrix(NA, ncol = length(strain_isolation_times), 
                   nrow = length(unique(titre_dat$individual)))
colnames(inf_data) <- c(strain_isolation_times)
inf_data[ , which(strain_isolation_times %in%  PCR_dates) + 1 ] <- as.matrix(infection_data_2016[ , -1])

#make ids numeric in titre_dat
titre_dat$individual <- as.numeric(titre_dat$individual)

####MCMC
par_tab <- read.csv("inputs/par_tab_inf_data.csv")
par_tab$values[which(par_tab$names == "wane")] <- par_tab$values[which(par_tab$names == "wane")] / resolution

## Fix all long term parameters to 0
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 #mu, tau, sigma1, and sigma2 are fixed
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 # set these values to 0

tmp <- par_tab[par_tab$names == "phi",]
for(i in 1:(length(strain_isolation_times)-1)){
  par_tab <- rbind(par_tab, tmp)
}

#delta will need to be fixed to make some sense out of results 
#but fix higher at 0.9
par_tab$fixed[par_tab$names == "delta"] <- 1
par_tab$values[par_tab$names == "delta"] <- 0.9


mcmc_pars <- c("iterations"=400000,"adaptive_period"=200000, "burnin"=50000,
               "hist_sample_prob"=0.5,"thin"=5,"thin_hist"=20,"swap_propn"=0.5,
               "hist_switch_prob"=0.2,"year_swap_propn"=0.5,
               "switch_sample" = 2)


res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
                antigenic_map = fit_dat, mcmc_pars =mcmc_pars,
                mvr_pars = NULL, inf_dat = inf_data, start_inf_hist = NULL, filename="chains/A_H1N1",
                CREATE_POSTERIOR_FUNC=create_posterior_func_fast, CREATE_PRIOR_FUNC=NULL,
                version=1,
                fast_version=TRUE)


res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
                antigenic_map = fit_dat, mcmc_pars =mcmc_pars,
                mvr_pars = NULL, inf_dat = NULL, start_inf_hist = NULL, filename="chains/A_H1N1_no_inf_data",
                CREATE_POSTERIOR_FUNC=create_posterior_func_fast, CREATE_PRIOR_FUNC=NULL,
                version=1,
                fast_version=TRUE)

which_phi <- which(par_tab$names == "phi")
high_risk_index <- c(30, 31, 32, 33, 42, 43, 44, 45, 54, 55, 56, 57)
low_risk_index <-   which_phi[-which(which_phi %in% high_risk_index) ]

par_tab$values[low_risk_index] <- 0.05 

create_prior_func <- function(par_tab){
  par_names <- par_tab$names
  
  prior <- function(pars){
    
    which_phi <- which(par_tab$names == "phi")
    high_risk_index <- c(30, 31, 32, 33, 42, 43, 44, 45, 54, 55, 56, 57)
    low_risk_index <-   which_phi[-which(which_phi %in% high_risk_index) ]
    
    delta_index<-which(par_names== "delta") 
    
    c <- dunif(pars[delta_index],0.7, 1,log=TRUE)

    a <- sapply(high_risk_index, function(x) dbeta(pars[x], shape1=3.1, shape2=2.2,log=TRUE) ) 
    b <- sapply(low_risk_index, function(x) dbeta(pars[x], shape1 = 1, shape2 = 50,log=TRUE) )

    return(sum(a,b, c))
  }
}


par_tab$fixed[par_tab$names == "delta"] <- 0
par_tab$fixed[par_tab$names == "error"] <- 1
res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
                antigenic_map = fit_dat, mcmc_pars =mcmc_pars,
                mvr_pars = NULL, inf_dat = inf_data, start_inf_hist = NULL, filename="chains/A_H1N1_prior",
                CREATE_POSTERIOR_FUNC=create_posterior_func_fast, CREATE_PRIOR_FUNC=create_prior_func,
                version=1,
                fast_version=TRUE)

par_tab$fixed[par_tab$names == "delta"] <- 1
res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
                antigenic_map = fit_dat, mcmc_pars =mcmc_pars,
                mvr_pars = NULL, inf_dat = NULL, start_inf_hist = NULL, filename="chains/A_H1N1_no_inf_data_prior",
                CREATE_POSTERIOR_FUNC=create_posterior_func_fast, CREATE_PRIOR_FUNC=create_prior_func,
                version=1,
                fast_version=TRUE)

#filename <- "A_H1N1_prior"
filename <- "A_H1N1_no_inf_data_prior"
chain <- read.csv(paste("~/Documents/GitHub/serosolver/chains/",filename,"_chain.csv", sep = ""))
chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
#chain <- chain[chain$sampno >= (max(unique(chain$sampno)) - 100000),] # if the chain is still too big to handle

head(chain)

inf_chain  <- data.table::fread(paste("~/Documents/GitHub/serosolver/chains/",filename,"_infection_histories.csv", sep = ""))
inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
#inf_chain <- inf_chain[inf_chain$sampno >= (max(unique(chain$sampno)) - 100000),] # if the chain is still too big to handle

n_alive <- get_n_alive(titre_dat, strain_isolation_times)
quantiles_res <- find_quantiles(inf_chain, titre_dat, 
                                strain_isolation_times)


plot_attack_rates(inf_chain, titre_dat, strain_isolation_times)


png(paste(filename, "likelihood.png", sep = ""), res = 200, units = "px", width = 1500, height = 1000)
par(mfrow =c (1, 1))
plot(chain$likelihood, type = "l")
dev.off()

png(paste(filename, "pars.png", sep = ""), res = 200, units = "px", width = 1500, height = 750)
par(mfrow = c(1, 3))
hist(chain$delta, xlab = "Sensitivity", main = "", col = "lightgrey", breaks = 15)
hist(chain$wane, xlab = "Waning", main = "",  col = "lightgrey",breaks = 15)
hist(chain$mu_short, xlab = "Short term boost", main = "",  col = "lightgrey", breaks = 15)
dev.off()

