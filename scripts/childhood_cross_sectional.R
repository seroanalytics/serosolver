# Set working directory and load code
code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()


## Antigenic map for cross reactivity parameters
fit_dat <- read.csv("data/map_UK.csv")
strain_isolation_times <- fit_dat$inf_years

## Read in parameter table
par_tab <- read.csv("~/Google Drive/Influenza childhood UK/par_tab_UK.csv")
tmp <- par_tab[par_tab$names == "phi",]
for(i in 1:(length(strain_isolation_times)-1)){
  par_tab <- rbind(par_tab, tmp)
}

## How many individuals
n_indiv <- 200

# What ages
min_age <- 5
max_age <- 11

# Constant attck rate
set.seed(54)
attack_rates <- rnorm(length(strain_isolation_times), 0.2, 0.05)
plot(attack_rates, type = "l", ylim = c(0, 1))
attack_rates_df <- data.frame( AR = attack_rates, inf_years = strain_isolation_times)


set.seed(54)
dat <- simulate_data(par_tab = par_tab, group = 1, n_indiv = n_indiv, buckets = 1,
                     strain_isolation_times, sampling_times = 2014, nsamps = 1,
                     antigenic_map = fit_dat,
                     age_min = min_age, age_max = max_age,
                     attack_rates = attack_rates)


# Recent panel only
titre_dat <- dat[[1]]
titre_dat <- titre_dat[titre_dat$virus %in% c(2014, 2013, 2012),]
dim(titre_dat) # 300 assays 
write.table(titre_dat,paste("data/recent_titre_dat.csv",sep="_"),row.names=FALSE,sep=",")

# Full panel (recent and some historic)
# Test 100 indivividuals against 1 strain (each year)
titre_dat <- dat[[1]]
titre_dat <- titre_dat[titre_dat$virus %in% c(2014, 2013, 2012, 2007, 2005, 2002),]
dim(titre_dat) #600 assays
write.table(titre_dat,paste("data/full_titre_dat.csv",sep="_"),row.names=FALSE,sep=",")

# ages
ages <- dat$ages
write.table(ages,paste("data/UK_ages.csv",sep="_"),row.names=FALSE,sep=",")

# infection history
inf_hist <- dat$infection_histories
write.table(inf_hist,paste("data/UK_inf_hist.csv",sep="_"),row.names=FALSE,sep=",")


######### mcmc ############

# Priors on antibody process parameters
create_prior_func <- function(par_tab){
  par_names <- par_tab$names
  
  prior <- function(pars){
    
    wane_index <- which(par_names == "wane")
    mu_index <- which(par_names == "mu")
    a_index <- which(par_names=="a") 
    sigma1_index <- which(par_names=='sigma1') #mean 0.1
    b_index<-which(par_names=='b') #mean is 0.2
    
    #a <- dgamma(pars[mu_index],shape=9,rate=5,log=TRUE)
    a <- dgamma(pars[mu_index], shape=32691.41,rate=18161.62,log=TRUE) #tight prior on long term boost
    b <- dbeta(pars[wane_index],shape1=6,shape2=1,log=TRUE)
    c <- dbeta(pars[a_index],shape1 = 60,shape2=60,log=TRUE)
    d <- dbeta(pars[sigma1_index],shape1 = 90,shape2=813,log=TRUE)
    e <- dbeta(pars[b_index],shape1 = 5,shape2=5,log=TRUE)
    
    return(sum(a,b,c,d,e))
  }
}

mcmc_pars <- c("iterations"=200000,"adaptive_period"= 200000, "burnin"= 10000,"hist_sample_prob"=0.5,"thin"=500,"thin_hist"=200,"swap_propn"=0.5,"hist_switch_prob"=0.2,"year_swap_propn"=0.5)

#Try fix waning 
par_tab$fixed[which(par_tab$names == "wane")] <- 1

# Recent panel
titre_dat <- read.csv("data/recent_titre_dat.csv")
titre_dat <- merge(titre_dat, ages)
res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat, 
                antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                mvr_pars = NULL, start_inf_hist = NULL, 
                filename = "chains/recent",
                CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = create_prior_func,
                version = 1,  
                fast_version = TRUE)



# Full panel
titre_dat <- read.csv("data/full_titre_dat.csv")
titre_dat <- merge(titre_dat, ages)
res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat, 
                antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                mvr_pars = NULL, start_inf_hist = NULL, 
                filename = "chains/full",
                CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = create_prior_func,
                version = 1,  
                fast_version = TRUE)




######### plots ############

make_plots(filename = "recent")
make_plots(filename = "full")




