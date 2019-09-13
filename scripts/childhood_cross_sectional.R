# Set working directory and load code
code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

library(plyr)
library(data.table)

## Antigenic map for cross reactivity parameters
fit_dat <- read.csv("data/map_UK.csv")
strain_isolation_times <- fit_dat$inf_years

## Read in parameter table
par_tab <- read.csv("inputs/parTab_wane_UK_children.csv")
tmp <- par_tab[par_tab$names == "phi",]
for(i in 1:(length(strain_isolation_times)-1)){
  par_tab <- rbind(par_tab, tmp)
}

## How many individuals
n_indiv <- 100

# What ages
min_age <- 5
max_age <- 11

# Constant attck rate
set.seed(50)
attack_rates <- rlnorm(length(strain_isolation_times), log(0.18), 0.5) # Mean = exp( log(0.18)+0.5^2/2 )
par(mar=c(3,3,1,1),mgp=c(2,0.7,0))
plot(strain_isolation_times,attack_rates, type = "l", ylim = c(0, 0.5),ylab="simulated attack rate",xlab="year")

dev.copy(pdf,paste("outputs/attack_rates.pdf",sep=""),width=6,height=4)
dev.off()


attack_rates_df <- data.frame( AR = attack_rates, inf_years = strain_isolation_times)


set.seed(50)
dat <- simulate_data(par_tab = par_tab, group = 1, n_indiv = n_indiv, buckets = 1,
                     strain_isolation_times, sampling_times = c(2015), nsamps = 1,
                     antigenic_map = fit_dat,
                     age_min = min_age, age_max = max_age,
                     attack_rates = attack_rates)


# Recent panel only
titre_dat <- dat[[1]]
titre_dat <- titre_dat[titre_dat$virus %in% c(2015, 2012, 2007),]
dim(titre_dat) # 300 assays 
write.table(titre_dat,paste("data/recent_titre_dat.csv",sep="_"),row.names=FALSE,sep=",")

# Full panel (recent and some historic)
# Test 100 indivividuals against 1 strain (each year)
titre_dat <- dat[[1]]
titre_dat <- titre_dat[titre_dat$virus %in% c(2015, 2014, 2012, 2009, 2007, 2005, 2002),]
#titre_dat <- titre_dat[titre_dat$virus %in% c(2002:2014),]
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
    
    a <- dgamma(pars[mu_index],shape=1000,rate=500,log=TRUE)
    #a <- dgamma(pars[mu_index], shape=32691.41,rate=18161.62,log=TRUE) #tight prior on long term boost
    b <- dbeta(pars[wane_index],shape1=6,shape2=1,log=TRUE)
    c <- dbeta(pars[a_index],shape1 = 60,shape2=60,log=TRUE)
    d <- dbeta(pars[sigma1_index],shape1 = 90,shape2=510,log=TRUE)
    e <- dbeta(pars[b_index],shape1 = 5,shape2=20,log=TRUE)
    
    return(sum(a,b,c,d,e))
  }
}

mcmc_pars <- c("iterations"=1e5,"adaptive_period"= 1e5, "burnin"= 1e4,"hist_sample_prob"=0.5,"thin"=200,"thin_hist"=100,"swap_propn"=0.5,"hist_switch_prob"=0.2,"year_swap_propn"=0.5)

# plot(seq(0,3,0.1),dgamma(seq(0,3,0.1), shape=100,rate=50,log=F))
# plot(seq(0,0.3,0.01),dbeta(seq(0,0.3,0.01), shape1=90, shape2=510,log=F))

#Try fix waning 
par_tab$fixed[which(par_tab$names == "wane")] <- 1

# Recent panel
titre_dat <- read.csv("data/recent_titre_dat.csv")
titre_dat <- merge(titre_dat, ages)
resR <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat,
                antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                mvr_pars = NULL, start_inf_hist = NULL,
                filename = "chains/recent",
                CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = create_prior_func,
                version = 1,
                fast_version = TRUE)



# Full panel
titre_dat <- read.csv("data/full_titre_dat.csv")
titre_dat <- merge(titre_dat, ages)
resF <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat, 
                antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                mvr_pars = NULL, start_inf_hist = NULL, 
                filename = "chains/full",
                CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = create_prior_func,
                version = 1,  
                fast_version = TRUE)




######### plots ############



# make_plots(filename = "recent")
# make_plots(filename = "full")

res <- resF

## Read in the MCMC chains and plot posteriors
chain <- read.csv(res$chain_file)
inf_chain <- data.table::fread(res$history_file)
# plot(coda::as.mcmc(chain[chain$sampno > 20000,c("mu","sigma1","a")]))

# Plot model predicted titres for a subset of individuals
plot_infection_histories(chain = chain,infection_histories = inf_chain,
                         titre_dat = titre_dat,individuals=c(1:4),
                         antigenic_map=fit_dat,par_tab=par_tab)


# Plot inference against true histories

# p_indiv_inf_hists <- generate_cumulative_inf_plots(inf_chain,indivs=1:9,
#                                                    pad_chain=FALSE,
#                                                    real_inf_hist=inf_hist,
#                                                    strain_isolation_times = strain_isolation_times,
#                                                    number_col=3)
# 
# print(p_indiv_inf_hists[[2]])

# Plot attack rates
# p1 <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times, n_alive = NULL,
#                               pointsize = 1, fatten = 1, pad_chain=TRUE, prior_pars=NULL,plot_den=FALSE) +
#       geom_line(data=attack_rates_df,aes(x=inf_years, y=AR))
# 
# p1

# Plot infection history comparison

histogram_comparison("full")


histogram_comparison <- function(file_name="") {
  
  indivs=1:n_indiv; pad_chain=FALSE; real_inf_hist=inf_hist; strain_isolation_times = strain_isolation_times; number_col=3; nsamp = 100; subset_years=NULL

  if(is.null(inf_chain$chain_no)){
    inf_chain$chain_no <- 1
  }  

  samps <- sample(unique(inf_chain$sampno), nsamp)
  inf_chain <- inf_chain[inf_chain$sampno %in% samps, ]
  
  ## Get number of probability that infected in a given time point per individual and year
  inf_chain1 <- inf_chain[inf_chain$i %in% indivs, ]
  if(!is.null(subset_years)) inf_chain1 <- inf_chain1[inf_chain1$j %in% subset_years,]
  data.table::setkey(inf_chain1, "i", "j","chain_no")
  #max_sampno <- length(unique(inf_chain1$sampno))
  max_sampno <- nsamp
  
  ## Number of samples with a 1 divided by total samples
  densities <- inf_chain1[, list(V1 = sum(x) / max_sampno), by = key(inf_chain1)]
  
  ## Convert to real time points
  densities$j <- as.numeric(strain_isolation_times[densities$j])
  densities$i <- as.numeric(densities$i)
  
  ## If someone wasn't infected in a given year at all, then need a 0
  strain_isolation_times1 <- strain_isolation_times
  if(!is.null(subset_years)) strain_isolation_times1 <- strain_isolation_times[subset_years]
  all_combos <- data.table(expand.grid(i = indivs, j = strain_isolation_times1, chain_no=unique(inf_chain$chain_no)))
  all_combos$j <- as.numeric(all_combos$j)
  all_combos$i <- as.numeric(all_combos$i)
  all_combos <- data.table::fsetdiff(all_combos[, c("i", "j","chain_no")], densities[, c("i", "j","chain_no")])
  all_combos$V1 <- 0
  densities <- rbind(all_combos, densities)
  
  infection_history1 <- NULL
  if (!is.null(real_inf_hist)) {
    infection_history1 <- as.data.frame(real_inf_hist)
    infection_history1 <- infection_history1[indivs, ]
    infection_history1$individual <- indivs
    
    colnames(infection_history1) <- c(strain_isolation_times, "i")
    infection_history1 <- reshape2::melt(infection_history1, id.vars = "i")
    infection_history1$variable <- as.numeric(as.character(infection_history1$variable))
    #infection_history1 <- infection_history1[infection_history1$value == 1, ]
  }
  
  # Sum over individuals 
  
  collect_hist <- NULL
  
  for(ii in 1:n_indiv){
    
    true_hist <- infection_history1 %>% filter(i==ii) %>% select(value) %>% sum()
    
    estimated_hist <- densities %>% filter(i==ii) %>% select(V1) %>% sum()
    
    collect_hist <- rbind(collect_hist,c(ii,true_hist,estimated_hist) )
    
  }
  
  collect_hist <- as_tibble(collect_hist); names(collect_hist) <- c("id","true_n","estimated_n")
  
  par(mfrow=c(1,1),mgp=c(2,0.7,0),mar=c(3,3,1,1))
  
  hist(collect_hist$estimated_n,breaks=seq(-0.5,8.5,1),xlab="number of infections",main="",col="grey",yaxs="i",ylim=c(0,40))
  text(x=4,y=15,labels = "Estimated",col="dark grey",adj=0)
  text(x=4,y=20,labels = "'True'",col=rgb(0,0.8,0),adj=0)
  hist(collect_hist$true_n,breaks=seq(-8.5,8.5,1),add=T,col=rgb(0,0.8,0,0.2))
  
  dev.copy(png,paste("outputs/histogram.png",sep=""),res=200,units="cm",width=14,height=10)
  dev.off()

  hist(collect_hist$estimated_n - collect_hist$true_n,breaks=seq(-8.5,8.5,1),xlab="estimated - true number of infections",main="",col="light blue",yaxs="i")
  
  output_table <- as.data.frame(table(abs(round(collect_hist$estimated_n - collect_hist$true_n)))/n_indiv )
  
  write_csv(output_table,paste0("outputs/accuracy_",file_name,"_n",n_indiv,"_age",max_age,".csv"))

}





