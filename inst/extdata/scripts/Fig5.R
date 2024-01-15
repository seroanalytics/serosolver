######################
## Author: 15.11.2019 - James Hay jameshay218@gmail.com
## Description: 
##  This script fits the serosolver antibody kinetics model to simulated data to
##  match the the Hong Kong HI titre data. Use this to create Figure 5 in the main text

library(serosolver)
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(tidyr)
library(doParallel)
library(foreach)
library(ggpubr)
library(bayesplot)
library(viridis)
library(ggthemes)
library(cowplot)
library(grid)
library(gridExtra)


# Set working directory
setwd("E:/James/Google Drive/Influenza/serosolver/methods_paper/PLOS Comp Biol/Results/case_study_1/")

# Setup -------------------------------------------------------------------
serosolver <- FALSE

## We'll be parallelising a few chains
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

# set seed
set.seed(0)
# filenames
filename <- "case_study_1_sim"
filenames <- paste0(filename, "_",1:5)

# Simulate data and setup inputs -------------------------------------------------------------------
# How many individuals to simulate?
n_indiv <- 311
n_samp <- 4 ## Number of blood samples per individual
buckets <- 4 ## Set to 1 for annual model. Greater than 1 gives subannual (eg. buckets = 2 is infection period every half year)
repeats <- 1 ## Number of repeats measurements per titre
# Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(par_tab_path, stringsAsFactors=FALSE)

par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1/3,1/3)
par_tab[par_tab$names == "wane","values"] <- 1
par_tab[par_tab$names == "wane","values"] <- par_tab[par_tab$names == "wane","values"]/buckets
par_tab[par_tab$names %in% c("tau","sigma1","sigma2"),"fixed"] <- 1 # mu, tau, and sigma are fixed
par_tab[par_tab$names %in% c("tau","sigma1","sigma2"),"values"] <- 0 # set value of mu and tau to 0
par_tab[par_tab$names == "MAX_TITRE","values"] <- 9 # set max titre to 9

par_tab <- par_tab[par_tab$names != "phi",] # remove lambda row of parameter table

# set par values to MLE from real data
#par_tab[par_tab$names %in% c("mu_short","wane","error"),"values"] <- c(4.02,0.03,0.79)
par_tab[par_tab$names %in% c("mu","mu_short","wane","error"),"values"] <- c(3.327, 2.841,  0.443, 0.629)

samplingTimes <- seq(2009*buckets + 1, 2012*buckets, by=1)
yearMin <- 2009*buckets + 1
yearMax <- 2012*buckets
age_min <- 6*buckets
age_max <- 6*buckets

# Read in and generate the antigenic map to read strain relationships from
strain_isolation_times <- seq(yearMin, yearMax, by=1)
simInfPars=c("mean"=0.2/buckets,"sd"=0.5,"bigMean"=0.5/buckets/2,"logSD"=1)
# use SIR
#attack_rates <- simulate_attack_rates(strain_isolation_times, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
attack_rates <- c(0.0643, 0.0225, 0.0129, 0.0547, 0.0514, 0.0386, 0.0643, 0.0129, 
                  0.00322, 0.0772, 0.0675, 0.011) # estimated from real data

# Simulate data
dat <- simulate_data(par_tab, group=1, n_indiv=n_indiv, buckets=buckets,
                     strain_isolation_times=strain_isolation_times,
                     sampling_times=samplingTimes, nsamps = n_samp, 
                     antigenic_map=NULL, 
                     titre_sensoring=0, 
                     age_min=age_min*buckets,age_max=age_max*buckets,
                     attack_rates=attack_rates,repeats = repeats)

mydat <- dat$data[dat$data$virus == 8037,]
titre_dat <- merge(mydat,dat$ages,by='individual') # add DOB column
indivs <- unique(titre_dat$individual) #all individuals
true_inf_hist <- dat$infection_histories

# Run MCMC framework ------------------------------------------------------
#setwd("Q:/serosolver_manuscript/case_study_1") # PC path
x <- filenames
# default mcmcPars (note: don't need to pass to serosolver unless you want non-default options)
mcmc_pars <- c("iterations"=1000000,"target_acceptance_rate_theta"=0.44,"target_acceptance_rate_inf_hist"=0.44,"adaptive_frequency"=1000,"thin"=100,"adaptive_iterations"=200000, 
              "save_block"=1000, "thin_inf_hist"=1000,"proposal_inf_hist_indiv_prop"=0.5,"proposal_ratio"=2, "burnin"=0, "proposal_inf_hist_time_prop"=1, 
              "proposal_inf_hist_distance"=3,"proposal_inf_hist_adaptive"=1,"proposal_inf_hist_indiv_swap_ratio"=0.5,"proposal_inf_hist_group_swap_ratio"=0.5,"proposal_inf_hist_group_swap_prop"=1)

if(serosolver) {
  #res <- for(j in 1:length(filenames)){
  res <- foreach(x = filenames) %dopar% {
  # starting parameter values
    start_prob <- -Inf
    while(!is.finite(start_prob)){
      start_tab <- par_tab
      for(i in 1:nrow(start_tab)){
        if(start_tab[i,"fixed"] == 0){
          start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"], 
                                         start_tab[i,"upper_start"])
        }
      }
      start_inf <- setup_infection_histories_titre(titre_dat, strain_isolation_times, space=3,titre_cutoff=4)
      f <- create_posterior_func(start_tab,titre_dat,fit_dat,version=2) # function in posteriors.R
      start_prob <- sum(f(start_tab$values, start_inf)[[1]])
    }
    
    res <- serosolver(par_tab = start_tab, titre_dat = titre_dat,antigenic_map = fit_dat,start_inf_hist=start_inf, 
                    mcmc_pars = mcmc_pars,
                    filename = paste0("chains_sim/",x), CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2)
  }
}
beepr::beep(4)

######################
# For Debugging #
######################
# start_tab[start_tab$names %in% c("mu","mu_short","tau","wane","sigma1","sigma2"),"values"] <- c(1.80, 2.70, 0.05, 0.25, 0.10, 0.03) # set value of mu and tau to 0
# f <- create_posterior_func(start_tab,titre_dat,fit_dat,version=2) # function in posteriors.R
# sum(f(start_tab$values, start_inf))
# save(start_tab,titre_dat,fit_dat,start_inf,file='forJames3.rData')

######################
# Processing outputs #
######################
## Read in the MCMC chains automatically
all_chains <- load_mcmc_chains(location="chains_sim",thin=10,burnin=200000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)

## Get the MCMC chain list
list_chains <- all_chains$theta_list_chains
## Look at diagnostics for the free parameters
list_chains1 <- lapply(list_chains, function(x) x[,c("mu","mu_short","wane","error",
                                                     "total_infections","lnlike","prior_prob")])
param_estimates <- round(apply(list_chains1[[1]],2,mean),2)
gelman.diag(as.mcmc.list(list_chains1))
gelman.plot(as.mcmc.list(list_chains1))
effectiveSize(as.mcmc.list(list_chains1))

color_scheme_set("viridis")
#pdf(paste0("plots/",filename,"_parameter_trace.pdf"))
mcmc_trace(list_chains1)
#dev.off()

inf_chain <- all_chains$inf_chain
#plot_infection_history_chains_time(inf_chain,pad_chain=FALSE)[[1]]
#plot_infection_history_chains_indiv(inf_chain,indivs = 201:320,pad_chain=FALSE)[[1]]

inf_chain1 <- inf_chain[inf_chain$chain_no == 1,]

# Plot inferred antibody titres
chain <- as.data.frame(all_chains$theta_chain)
chain1 <- chain[chain$chain_no == 1,]
#chain <- read.csv("chains/case_study_1_1_chain.csv")
par_tab <- par_tab[par_tab$names != "phi",]
rand_indivs <- sample(unique(titre_dat$individual),10*5)
rand_indivs <- rand_indivs[order(rand_indivs)]
rand_indivs <- c(19,43,130,298)
titre_preds <- get_titre_predictions(chain = chain1, infection_histories = inf_chain1, titre_dat = titre_dat, 
                                     individuals = rand_indivs,nsamp = 1000,
                                     antigenic_map = NULL, strain_isolation_times = strain_isolation_times,
                                     par_tab = par_tab,expand_titredat = TRUE)
to_use <- titre_preds$predicted_observations
model_preds <- titre_preds$predictions
to_use$individual <- rand_indivs[to_use$individual]

labels2 <- c("2009-Q1","2009-Q2","2009-Q3","2009-Q4",
             "2010-Q1","2010-Q2","2010-Q3","2010-Q4",
             "2011-Q1","2011-Q2","2011-Q3","2011-Q4")

inf_hist_densities <- titre_preds$histories
inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
inf_hist_densities$xmax <- inf_hist_densities$variable+0.5

true_inf_hist_melt <- melt(true_inf_hist)
colnames(true_inf_hist_melt) <- c("individual","variable","inf")
true_inf_hist_melt <- true_inf_hist_melt[true_inf_hist_melt$inf == 1,]
true_inf_hist_melt$variable <- strain_isolation_times[true_inf_hist_melt$variable]
true_inf_hist_melt <- true_inf_hist_melt[true_inf_hist_melt$individual %in% rand_indivs,]


inf_hist_densities_show <- inf_hist_densities[inf_hist_densities$value > 0.001,]

titre_pred_p <- ggplot(to_use) +
  #geom_vline(xintercept=strain_isolation_times,col="grey70",size=0.02,alpha=0.5) +
  #geom_hline(yintercept=0:10,col="grey70",size=0.02,alpha=0.5) +
  geom_rect(data=inf_hist_densities_show,
            aes(xmin=xmin,xmax=xmax,fill=value),ymin=-1,ymax=11)+
  geom_ribbon(aes(x=samples,ymin=lower, ymax=upper),alpha=0.4, fill="#56B4E9",size=0.2)+
  geom_ribbon(data=model_preds[model_preds$individual %in% rand_indivs,], 
              aes(x=samples,ymin=lower,ymax=upper),alpha=0.7,fill="#56B4E9",size=0.2) + 
  geom_line(data=model_preds, aes(x=samples, y=median),linetype="dotted",color="grey10")+
  geom_rect(ymin=9,ymax=11,xmin=0,xmax=9000,fill="grey70")+
  geom_rect(ymin=-2,ymax=0,xmin=0,xmax=9000,fill="grey70")+
  scale_x_continuous(expand=c(0,0),labels=labels2[seq(1,12,by=2)],breaks=strain_isolation_times[seq(1,12,by=2)]) +
  scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
  guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                              barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
  geom_vline(data=true_inf_hist_melt,aes(xintercept=variable)) +
  geom_point(data=titre_dat[titre_dat$individual %in% rand_indivs,], aes(x=samples, y=titre),fill="#009E73",size=1,shape=23,col="black")+
  ylab("log titre") +
  xlab("Time of virus circulation") +
  theme_pubr()+
  theme(strip.background = element_blank(),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        legend.margin = margin(-1,-1,-3,-1),
        #legend.position="none",
        strip.text=element_blank(),
        #panel.grid.major=element_line(colour="grey70",size=0.1),
        #panel.grid.minor=element_line(colour="grey70",size=0.1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        plot.margin=margin(r=5,t=5,l=5))+
  coord_cartesian(ylim=c(0,10),xlim=c(8036.5,8048.5)) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  facet_wrap(~individual,ncol=2)
titre_pred_p


###
## Plot attack rates
###
inf_chain <- all_chains$inf_chain
inf_chain <- pad_inf_chain(inf_chain)
n_alive <- get_n_alive(titre_dat, strain_isolation_times)
data.table::setkey(inf_chain, "samp_no", "j","chain_no")
tmp <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.1, 0.5, 0.9,  0.975)))
colnames(quantiles) <- c("j", "lower", "lower_50","median","upper_50","upper")
quantiles[c("lower", "lower_50","median","upper_50","upper")] <- quantiles[c("lower", "lower_50","median","upper_50","upper")]/n_alive[quantiles$j]
quantiles$year <- strain_isolation_times[quantiles$j]
quantiles$taken <- quantiles$year %in% unique(titre_dat$samples)
quantiles$vac_status <- c(rep('Unvaccinated',dim(quantiles)[1]))

quantiles_all <- quantiles

## Colour depending on vac_status
colour_fills_unvac <- c("#009E73","#56B4E9")
colour_fills_age <- c("#CC79A7","#0072B2","#D55E00")

ymax = 0.15

ar <- data.frame("year"=strain_isolation_times,"ar"=colSums(dat$infection_histories)/n_alive)
quantiles_all$true <- ar$ar
tmp_quantiles_median <- melt(quantiles_all[,c("median","true","year")],id.vars=c("year"))


p <- ggplot(quantiles_all) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper), alpha = 0.25,fill="#56B4E9") +
  geom_ribbon(aes(x = year, ymin = lower_50, ymax = upper_50), alpha = 0.5,fill="#56B4E9") +
  geom_line(data=tmp_quantiles_median, aes(x = year, y = value,col=variable,group=variable),size=0.75) +
  geom_point(data=tmp_quantiles_median, aes(x = year, y = value,col=variable), size = 0.75,col="#56B4E9") +
  #geom_line(data=ar,aes(x=year,y=ar),col="#009E73",size=1,linetype="dashed") +
  geom_point(data=ar,aes(x=year,y=ar),fill="#009E73",size=2,shape=23,col="black") +
  scale_y_continuous(limits = c(-0.005, ymax), expand = c(0, 0),breaks=seq(0,ymax,by=0.05)) +
  scale_x_continuous(expand = c(0, 0.1), breaks = strain_isolation_times, labels = labels2) +
  theme_pubr() +
  scale_color_manual(name = '', values=c("#56B4E9", "#009E73"),labels = c('Median estimate','True AR')) +
  ylab("Estimated per capita\n incidence (per quarter)") +
  xlab("Time of virus circulation") +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        #legend.key=element_rect(color=NA),
        #legend.title = element_blank(),
        legend.text=element_text(size=8,family="sans"),
        legend.position = c(0.7,0.99),
        legend.direction = "horizontal",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin=margin(l=10,r=5,t=5))
p

## Theta plot
theta_chain <- as.data.frame(all_chains$theta_chain)
theta_chain <- melt(theta_chain,id.vars=c("samp_no","chain_no"))
theta_chain$variable <- as.character(theta_chain$variable)

#ranges <- data.frame("variable"=c("mu","tau","sigma1"),"lower"=c(1.5,0,0.06),"upper"=c(3,0.1,0.13))

par_key <- c("mu"="mu[l]","mu_short"="mu[s]","wane"="omega","error"="epsilon","total_infections"="sum(Z[i])")
#par_key <- c("mu_short","wane","error")

theta_chain$variable <- par_key[theta_chain$variable]
theta_chain <- theta_chain[theta_chain$variable %in% par_key,]
quantiles <- ddply(theta_chain, ~variable, function(x) quantile(x$value, c(0.025,0.5,0.975)))
quantiles_melted <- melt(quantiles)
colnames(quantiles_melted) <- c("variable","quantile","val")
#quantiles_melted$variable <- as.character(quantiles_melted$variable)
quantiles_melted2 <- quantiles_melted[quantiles_melted$variable%in%par_key,]
quantiles_melted2$variable <- factor(quantiles_melted2$variable,levels=par_key)
theta_chain$variable <- factor(theta_chain$variable, levels=par_key)

use_vars <- par_key

correct_vals <- par_tab[par_tab$names %in% names(par_key),c("names","values")]
correct_vals$names <- par_key[correct_vals$names]
colnames(correct_vals) <- c("variable","values")
correct_vals <- rbind(correct_vals, data.frame("variable"="sum(Z[i])","values"=sum(true_inf_hist)))
correct_vals$variable <- factor(correct_vals$variable, levels=par_key)

densities <- ddply(theta_chain,~variable, function(x) {
  tmp <- density(x$value)
  ret <- data.frame(tmp$x,tmp$y)
})
densities <- densities[densities$variable != "sum(Z[i])",]
sum_inf_table <- table(theta_chain[theta_chain$variable == "sum(Z[i])","value"])
densities_sum_inf <- data.frame("variable"="sum(Z[i])","tmp.x"=as.numeric(names(sum_inf_table)),"tmp.y"=as.numeric(sum_inf_table)/sum(sum_inf_table))
densities <- rbind(densities, densities_sum_inf)
probs <- c(0,0.025,0.975,1)
quantiles <- ddply(theta_chain, ~variable, function(x) quantile(x$value, probs))
factorised <- c()
median_bits_x <- NULL
median_bits_y <- NULL

for(par in par_key){
  tmp_den <- densities[densities$variable == par,]
  factorised <- c(factorised, factor(findInterval(tmp_den$tmp.x,as.numeric(quantiles[quantiles$variable == par,2:4]))))
  tmp_den_func <- approxfun(density(theta_chain[theta_chain$variable == par,"value"]))
  tmp_median <- quantile(theta_chain[theta_chain$variable == par,"value"], 0.5)
  y_tmp <- tmp_den_func(tmp_median)
  median_bits_x[par] <- tmp_median
  median_bits_y[par] <- y_tmp
}
densities$quant <- factorised
densities[densities$variable == "sum(Z[i])","quant"] <- 1
median_segments <- data.frame(variable=par_key,x=median_bits_x,y=median_bits_y)

hacked_plots <- NULL
for(par in par_key){
  hacked_plots[[par]] <-  ggplot() + 
    #geom_blank(data=ranges,aes(xmin=lower,xmax=upper)) +
    geom_ribbon(data=densities[densities$variable == par,], aes(ymin=0,ymax=tmp.y,x=tmp.x),fill="grey80",col="black") +
    geom_ribbon(data=densities[densities$quant %in% c(3) & densities$variable == par,], 
                aes(ymin=0,ymax=tmp.y,x=tmp.x),fill="grey60",col="black") +
    geom_vline(data=correct_vals[correct_vals$variable == par,], aes(xintercept=values),size=1,col="#009E73") +
    geom_linerange(data=median_segments[median_segments$variable == par,],aes(x=x,ymin=0,ymax=y)) +
    scale_y_continuous(expand=c(0,0,0.05,0)) +
    ylab("Posterior density") +
    xlab("Value") +
    theme_pubr() +
    facet_wrap(~variable, labeller=label_parsed,scales="free",ncol=2) +
    theme(
      strip.background = element_blank(),
      strip.text=element_text(hjust=1),
      strip.text.x = element_text(size=12,family="sans",hjust=0.5),
      axis.text.x=element_text(size=8),
      axis.text.y=element_text(size=8),
      axis.title=element_blank(),
      plot.margin=margin(l=0,t=0,b=0,r=-10)
    ) 
}

par <- "sum(Z[i])"
hacked_plots[["mu[s]"]] <- hacked_plots[["mu[s]"]] + theme(plot.margin=margin(l=-5))
hacked_plots[["sum(Z[i])"]] <-  ggplot() + 
  #geom_blank(data=ranges,aes(xmin=lower,xmax=upper)) +
  geom_bar(data=densities[densities$variable == par,], aes(y=tmp.y,x=tmp.x),fill="grey80",col="black",stat="identity") +
  geom_vline(data=correct_vals[correct_vals$variable == par,], aes(xintercept=values),size=1,col="#009E73") +
  geom_vline(data=median_segments[median_segments$variable == par,],aes(xintercept=x)) +
  scale_y_continuous(expand=c(0,0,0.05,0)) +
  ylab("Posterior density") +
  xlab("Value") +
  theme_pubr() +
  facet_wrap(~variable, labeller=label_parsed,scales="free",ncol=2) +
  theme(
    strip.background = element_blank(),
    strip.text=element_text(hjust=1),
    strip.text.x = element_text(size=12,family="sans",hjust=0.5),
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=8),
    axis.title=element_blank(),
    plot.margin=margin(l=-5,t=0,b=0,r=5)
  ) 

p3 <- plot_grid(hacked_plots[[1]],hacked_plots[[2]],hacked_plots[[3]],hacked_plots[[5]],align="hv",axis="tlbr")
p3
y.grob <- textGrob("Posterior density", 
                   gp=gpar(fontsize=10), rot=90)

x.grob <- textGrob("Value", 
                   gp=gpar(fontsize=10))
p3 <- grid.arrange(arrangeGrob(p3, left = y.grob, bottom = x.grob))

sim_p <- plot_grid(p, plot_grid(p3, titre_pred_p,align="v",axis="tb",labels=c("B","C")),ncol=1,labels=c("A",""))
sim_p

svg("Fig5.svg",width=6.5,height=5.5)
plot(sim_p)
dev.off()
cairo_pdf("Fig5.pdf",width=6.5,height=5.5)
plot(sim_p)
dev.off()


