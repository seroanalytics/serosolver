######################
## KYLIE AINSLIE 25.10.2018 - ainslie.kylie@gmail.com
## Updated 15.11.2019 - jameshay218@gmail.com
##      Changed case study to estimate long and short term antibody boosting
## This script fits the serosolver antibody kinetics model to the Hong Kong HI titre data
## First, we fit to unvaccinated people and then to vaccinated to see differences in antibody
## kinetics parameter estimates (if any)
## This particular script uses a gibbs proposal step to resample infection histories
## which integrates out the annual force of infection parameters.
library(serosolver)
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(tidyr)
library(doMC)
library(doParallel)
library(foreach)
library(ggpubr)
library(bayesplot)
library(viridis)
library(ggthemes)
library(cowplot)
#library(doRNG)

serosolver <- TRUE

## Filename prefix for all saved outputs
filename <- "case_study_1_test"   # The general output filename
filenames <- paste0(filename, "_",1:5)

## We'll be parallelising a few chains
 cl <- makeCluster(detectCores(), type='PSOCK')
#registerDoParallel(cl)
registerDoMC(cores=5)

## Set working directory and load code from serosolver package
#  devtools::load_all("Q:/serosolver")
 #devtools::load_all("~/Documents/serosolver")

# read in data
# setwd("/Volumes/kainslie/serosolver_manuscript/case_study_1") # Mac path
# setwd("Q:/serosolver_manuscript/case_study_1") # PC path
setwd("~/Drive/Influenza/serosolver/methods_paper/PLOS Comp Biol/Results/case_study_1/")

## doRNG handles seeding in foreach loop
set.seed(0)

buckets <- 4
n_indiv <- 311   # How many individuals to fit to?

## Read in titre data (unvac)
inputDat <- read.csv(file = "../../../case_study_1/data/HKdata_h1n1_unvac.csv",header = TRUE)

indivs <- unique(inputDat$individual) #all individuals
# Subset data for indivs
titre_dat <- inputDat[inputDat$individual %in% indivs,c("individual","virus","titre","samples","DOB")]
titre_dat$individual <- match(titre_dat$individual, indivs)
titre_dat$group <- c(rep(1,nrow(titre_dat)))

titre_dat <- unique(titre_dat)
titre_dat <- plyr::ddply(titre_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))

# Read in titre data (vac)
inputDat2 <- read.csv(file = "../../../case_study_1/data/HKdata_h1n1_vac.csv",header = TRUE)

indivs2 <- unique(inputDat2$individual) #all individuals
# Subset data for indivs
titre_dat2 <- inputDat2[inputDat2$individual %in% indivs2,c("individual","virus","titre","samples","DOB")]
titre_dat2$individual <- match(titre_dat2$individual, indivs2)
titre_dat2$group <- c(rep(1,nrow(titre_dat2)))

titre_dat2 <- unique(titre_dat2)
titre_dat2 <- plyr::ddply(titre_dat2,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))

# Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
par_tab <- read.csv("../../../case_study_1/data/parTab_base.csv",stringsAsFactors=FALSE)
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1/3,1/3)

# Read in and generate the antigenic map to read strain relationships from
antigenicMap <- read.csv("../../../case_study_1/data/fonville_map_approx.csv",stringsAsFactors=FALSE)

fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- fit_dat[fit_dat$inf_years>=(2009*buckets+1) & fit_dat$inf_years<=2012*buckets,]
fit_dat[2:nrow(fit_dat),c('x_coord','y_coord')] <- fit_dat[1,c('x_coord','y_coord')] # make all coordinates the same

strain_isolation_times <- unique(fit_dat$inf_years)

# Multivariate proposals or univariate? Use univariate for now
covMat <- diag(nrow(par_tab))
scale <- 0.00005
w <- 1
mvrPars <- list(covMat, scale, w)
mvrPars <- NULL
# unique rows for each individual

unique_indiv <- titre_dat[!duplicated(titre_dat$individual),]
ageMask <- create_age_mask(unique_indiv$DOB, strain_isolation_times)
mcmcPars <- c("iterations"=500000,"target_acceptance_rate_theta"=0.44,"target_acceptance_rate_inf_hist"=0.44,"adaptive_frequency"=1000,"thin"=10,"adaptive_iterations"=100000, 
              "save_block"=1000, "thin2"=100,"histSampleProb"=1,"proposal_ratio"=2, "burnin"=0, "proposal_inf_hist_time_prop"=0.5, 
              "moveSize"=3,"histOpt"=1,"swapPropn"=0.5,"proposal_inf_hist_group_swap_ratio"=0.5,"proposal_inf_hist_group_swap_prop"=1)

################################
## Unvaccinated participants
################################
# starting parameter values 
#res <- for(j in 1:length(filenames)){

if(serosolver) {
  res <- foreach(x = filenames) %dopar% {
    ## Not all random starting conditions return finite likelihood, so for each chain generate random
    ## conditions until we get one with a finite likelihood
    start_prob <- -Inf
    while(!is.finite(start_prob)){
      if (is.na(start_prob)){
        cat('start_prob is NA')
        break
      }
      startTab <- par_tab
      for(i in 1:nrow(startTab)){
        if(startTab[i,"fixed"] == 0){
          startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                        startTab[i,"upper_start"])
        }
      }
      startInf <- setup_infection_histories_new_2(titre_dat, strain_isolation_times, space=3,titre_cutoff=4)
      startTab1 <- startTab[startTab$names != "phi",] # remove phi row of parameter table
      startTab1[startTab1$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 # mu, tau, and sigma are fixed
      startTab1[startTab1$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 # set value of mu and tau to 0
      startTab1[startTab1$names == "MAX_TITRE","values"] <- max(titre_dat$titre) # set max titre to 9
      f <- create_posterior_func(startTab1,titre_dat,fit_dat,version=2) # function in posteriors.R
      start_prob <- sum(f(startTab1$values, startInf)[[1]])
    }
    res_unvac <- serosolver(par_tab = startTab1, titre_dat = titre_dat,antigenic_map = fit_dat,start_inf_hist=startInf, 
                    mcmc_pars = mcmcPars,
                    filename = paste0("chains_unvac_simple/",x), CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2, temp=1)
  }
  beepr::beep(4)

################################
## Vaccinated participants
################################
# starting parameter values 
  #res <- for(j in 1:length(filenames)){
  res <- foreach(x = filenames) %dopar% {
    ## Not all random starting conditions return finite likelihood, so for each chain generate random
    ## conditions until we get one with a finite likelihood
    start_prob <- -Inf
    while(!is.finite(start_prob)){
      if (is.na(start_prob)){
        cat('start_prob is NA')
        break
      }
      startTab <- par_tab
      for(i in 1:nrow(startTab)){
        if(startTab[i,"fixed"] == 0){
          startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                        startTab[i,"upper_start"])
        }
      }
      startInf <- setup_infection_histories_new_2(titre_dat2, strain_isolation_times, space=3,titre_cutoff=4)
      startTab1 <- startTab[startTab$names != "phi",] # remove phi row of parameter table
      startTab1[startTab1$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 # mu, tau, and sigma are fixed
      startTab1[startTab1$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 # set value of mu and tau to 0
      startTab1[startTab1$names == "MAX_TITRE","values"] <- max(titre_dat2$titre) # set max titre to 9
      f <- create_posterior_func(startTab1,titre_dat2,fit_dat,version=2) # function in posteriors.R
      start_prob <- sum(f(startTab1$values, startInf)[[1]])
    }
    res_vac <- serosolver(par_tab = startTab1, titre_dat = titre_dat2,antigenic_map = fit_dat,start_inf_hist=startInf, 
                        mcmc_pars = mcmcPars,
                        filename = paste0("chains_vac_simple/",x), CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2, temp=1)
  }
  beepr::beep(4)
}

## Read in the MCMC chains automatically
all_chains <- load_mcmc_chains(location="chains_unvac",thin=10,burnin=100000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)

all_chains_vac <- load_mcmc_chains(location="chains_vac",thin=10,burnin=100000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)

## Get the MCMC chain list
list_chains <- all_chains$theta_list_chains
## Look at diagnostics for the free parameters
list_chains1 <- lapply(list_chains, function(x) x[,c("mu_short","wane","error","total_infections","lnlike")])
gelman.diag(as.mcmc.list(list_chains1))
gelman.plot(as.mcmc.list(list_chains1))
effectiveSize(as.mcmc.list(list_chains1))
summary(as.mcmc.list(list_chains1))

## Get the MCMC chain list
list_chains_vac <- all_chains_vac$theta_list_chains
## Look at diagnostics for the free parameters
list_chains_vac1 <- lapply(list_chains_vac, function(x) x[,c("mu_short","wane","error","total_infections","lnlike")])
gelman.diag(as.mcmc.list(list_chains_vac1))
gelman.plot(as.mcmc.list(list_chains_vac1))
effectiveSize(as.mcmc.list(list_chains_vac1))
summary(as.mcmc.list(list_chains_vac1))


## Plot the MCMC trace
color_scheme_set("viridis")

#pdf(paste0("plots/",filename,"_parameter_trace.pdf"))
mcmc_trace(list_chains1)
#dev.off()


chain1 <- all_chains$theta_chain

######################
# Summary statistics #
######################
# parameter estimates
myresults <- matrix(c(rep(0,4*7)),nrow=4)
rownames(myresults) <- c("mu_short","wane","error")
colnames(myresults) <- c("mean","sd","2.5%","25%","50%","75%","97.5%")

myresults[,"mean"] <- round(apply(chain1[,c("mu","mu_short","wane","error")],2,mean),3)
myresults[,"sd"] <- round(apply(chain1[,c("mu","mu_short","wane","error")],2,sd),3)  
myresults[,3:7] <- t(round(apply(chain1[,c("mu","mu_short","wane","error")],2,quantile,probs=c(0.025,0.1,0.5,0.9,0.975)),3))  
myresults      
###############
## Combine inferred infection histories with meta data
###############
inf_chain <- all_chains$inf_chain
is <- unique(titre_dat$individual)
js <- unique(inf_chain$j)
samp_nos <- unique(inf_chain$samp_no)
chain_nos <- unique(inf_chain$chain_no)
expanded_values <- data.table::CJ(
  i = is,
  j = js,
  samp_no = samp_nos,
  chain_no = chain_nos
)
diff_infs <- fsetdiff(expanded_values, inf_chain[, c("i", "j", "samp_no","chain_no")])
diff_infs$x <- 0
inf_chain <- rbind(inf_chain, diff_infs)

data.table::setkey(inf_chain, "i", "samp_no","chain_no")
n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
setkey(n_inf_chain_i, "i")
n_inf_chain <- n_inf_chain_i[,list(median_infs=median(V1)), 
                             by=key(n_inf_chain_i)]
colnames(n_inf_chain)[1] <- "individual"
setkey(n_inf_chain, "individual")
inputDat <- data.table(inputDat)
setkey(inputDat, "individual")
wow <- merge(n_inf_chain, inputDat)
wow$Age_1 <- factor(wow$Age_1, levels=c("<19","19-64",">64"))
age_ids <- unique(wow[,c("individual","Age_1")])

age_dist <- ggplot(wow) + 
  geom_boxplot(aes(group=Age_1,y=median_infs,x=Age_1)) +
  theme_bw() +
  ylab("Distribution of median number of inferred infections") +
  xlab("Age group")


########################
## PLOTS
#########################
## Plot MCMC trace on attack rates
n_alive <- get_n_alive_group(titre_dat, strain_isolation_times,melt_dat=TRUE)
n_alive
inf_chain <- all_chains$inf_chain
indivs <-unique(titre_dat[titre_dat$samples %in% strain_isolation_times[length(strain_isolation_times)],"individual"])
#pdf(paste0("plots/",filename,"_attack_rate_trace.pdf"))
plot_infection_history_chains_time(inf_chain,0,NULL,n_alive=n_alive,pad_chain=FALSE)[[1]]
#dev.off()

## MCMC trace plot on individual total infections
#pdf(paste0("plots/",filename,"_indiv_hist_trace.pdf"))
plot_infection_history_chains_indiv(all_chains$inf_chain,burnin = 0,1:25,pad_chain=FALSE)[[1]]
#dev.off()

## Generate posterior distributions for infection histories
y <- generate_cumulative_inf_plots(all_chains$inf_chain,burnin = 0,1:10,nsamp=100,
                                   strain_isolation_times = strain_isolation_times,
                                   pad_chain=FALSE,number_col = 2,subset_years = NULL)

#pdf(paste0("plots/",filename,"_cumu_infhist.pdf"))
plot(y[[1]])
#dev.off()

#pdf(paste0("plots/",filename,"_infhist.pdf"))
plot(y[[2]])
#dev.off()


# Plot inferred antibody titres
chain <- as.data.frame(all_chains$theta_chain)
chain1 <- chain[chain$chain_no == 1,]
#chain <- read.csv("chains/case_study_1_1_chain.csv")
par_tab <- par_tab[par_tab$names != "phi",]
inf_chain1 <- inf_chain[inf_chain$chain_no == 1,]
#rand_indivs <- sample(unique(titre_dat$individual),6)
#rand_indivs <- rand_indivs[order(rand_indivs)]
#rand_indivs <- c(2,10,21,36,95,195)
rand_indivs <- c(2,21,36,195)
titre_preds <- get_titre_predictions(chain = chain1, infection_histories = inf_chain1, titre_dat = titre_dat, 
                                     individuals = rand_indivs,nsamp = 1000,
                                    antigenic_map = fit_dat, par_tab = par_tab,expand_titredat = TRUE)
to_use <- titre_preds$predicted_observations
model_preds <- titre_preds$predictions
to_use$individual <- rand_indivs[to_use$individual]

labels2 <- c("2009-Q1","2009-Q2","2009-Q3","2009-Q4",
             "2010-Q1","2010-Q2","2010-Q3","2010-Q4",
             "2011-Q1","2011-Q2","2011-Q3","2011-Q4")

inf_hist_densities <- titre_preds$histories
inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
inf_hist_densities$xmax <- inf_hist_densities$variable+0.5

#inf_hist_densities_show <- inf_hist_densities[inf_hist_densities$value > 0.001,]

titre_pred_p <- ggplot(to_use) +
  geom_rect(data=inf_hist_densities,
            aes(xmin=xmin,xmax=xmax,fill=value),ymin=-1,ymax=11)+
  #geom_vline(xintercept=strain_isolation_times,col="grey70",size=0.02,alpha=0.5) +
  #geom_hline(yintercept=0:10,col="grey70",size=0.02,alpha=0.5) +
  geom_ribbon(aes(x=samples,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",size=0.2)+
  geom_ribbon(data=model_preds[model_preds$individual %in% rand_indivs,], 
              aes(x=samples,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",size=0.2) + 
  geom_line(data=model_preds, aes(x=samples, y=median),linetype="dotted",color="grey10")+
  geom_rect(ymin=9,ymax=11,xmin=0,xmax=9000,fill="grey70")+
  geom_rect(ymin=-2,ymax=0,xmin=0,xmax=9000,fill="grey70")+
  scale_x_continuous(expand=c(0,0),labels=labels2[seq(1,12,by=2)],breaks=fit_dat$inf_years[seq(1,12,by=2)]) +
  scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
  guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                              barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
  geom_point(data=titre_dat[titre_dat$individual %in% rand_indivs,], aes(x=samples, y=titre),shape=23, 
             col="black",size=1,fill=viridis(1)[1])+
  ylab("log titre") +
  xlab("Time of virus circulation") +
  theme_pubr()+
  theme(strip.background = element_blank(),
        legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.margin = margin(-1,-1,-3,-1),
        strip.text=element_blank(),
        #panel.grid.major=element_line(color="grey70",size=0.1),
        #panel.grid.major=element_line(colour="grey70",size=0.1),
        #panel.grid.minor=element_line(colour="grey70",size=0.1),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        plot.margin=margin(r=15,t=5,l=5))+
  coord_cartesian(ylim=c(0,10),xlim=c(8036.5,8048.5)) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  facet_wrap(~individual,ncol=2)
titre_pred_p



# plot unvac and vac on same plot
# unvaccinated
inf_chain <- all_chains$inf_chain
inf_chain <- pad_inf_chain(inf_chain)
n_alive <- get_n_alive(titre_dat, strain_isolation_times)
data.table::setkey(inf_chain, "samp_no", "j","chain_no")
tmp <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.25, 0.5, 0.75,  0.975)))
colnames(quantiles) <- c("j", "lower", "lower_50","median","upper_50","upper")
quantiles[c("lower", "lower_50","median","upper_50","upper")] <- quantiles[c("lower", "lower_50","median","upper_50","upper")] / n_alive[quantiles$j]
quantiles$year <- strain_isolation_times[quantiles$j]
quantiles$taken <- quantiles$year %in% unique(titre_dat$samples)
quantiles$vac_status <- c(rep('Unvaccinated',dim(quantiles)[1]))

# vaccinated
inf_chain2 <- all_chains_vac$inf_chain
inf_chain2 <- inf_chain2[inf_chain2$chain_no == 1,]
inf_chain2 <- pad_inf_chain(inf_chain2)
n_alive2 <- get_n_alive(titre_dat2, strain_isolation_times)
data.table::setkey(inf_chain2, "samp_no", "j","chain_no")
tmp <- inf_chain2[, list(V1 = sum(x)), by = key(inf_chain2)]

quantiles2 <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.1, 0.5, 0.9,  0.975)))
colnames(quantiles2) <- c("j", "lower", "lower_50","median","upper_50","upper")
quantiles2[c("lower", "lower_50","median","upper_50","upper")] <- quantiles2[c("lower", "lower_50","median","upper_50","upper")] / n_alive2[quantiles2$j]
quantiles2$year <- strain_isolation_times[quantiles2$j]
quantiles2$taken <- quantiles2$year %in% unique(titre_dat2$samples)
quantiles2$vac_status <- c(rep('Vaccinated',dim(quantiles2)[1]))

quantiles_all <- rbind(quantiles,quantiles2)
## Colour depending on vac_status
colour_fills_unvac <- c("#E69F00","#0072B2")
colour_fills_age <- c("#CC79A7","#009E73","#56B4E9")

strain_isolation_times1 <- strain_isolation_times 

ymax <- 0.5
quantiles_all$vac_status <- factor(quantiles_all$vac_status, levels=c("Vaccinated","Unvaccinated"))

p <- ggplot(quantiles_all) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = vac_status), alpha = 0.25) +
  geom_ribbon(aes(x = year, ymin = lower_50, ymax = upper_50, fill = vac_status), alpha = 0.5) +
  geom_line(aes(x = year, y = median, colour = vac_status),size=0.75) +
  geom_point(aes(x = year, y = median, colour = vac_status), size = 0.75) +
  scale_y_continuous(limits = c(-0.005, ymax), expand = c(0, 0),breaks=seq(0,ymax,by=0.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = strain_isolation_times1, labels = labels2,
                     limits=c(min(strain_isolation_times-0.1),max(strain_isolation_times+0.1))) +
  theme_pubr() +
  scale_fill_manual(values=colour_fills_unvac) +
  scale_color_manual(values=colour_fills_unvac) +
  ylab("Estimated per capita\n incidence (per quarter)") +
  xlab("Time of virus circulation") +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        axis.title=element_text(size=10),
        legend.title = element_blank(),
        legend.text=element_text(size=8,family="sans"),
        legend.position = c(0.7,0.99),
        legend.direction = "horizontal",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key=element_rect(color=NA),
        legend.background = element_blank(),
        legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin=margin(l=10,r=5,t=5))
p


########################
## Plot AR by age group in unvac
########################
### Plot inferred attack rates
# arguments
year_range <- min(strain_isolation_times):max(strain_isolation_times)
# plot_attack_rates_monthly function altered to add additional line for simulated attack rates (to compare against inferred attach rates)
months <- 1:buckets
years <- range(floor(year_range/buckets))
years <- years[1]:years[2]
labels <- c(sapply(years, function(x) paste0(months, "/",x)))
labels1 <- labels[1:length(year_range)]
labels1 <- labels1[seq(1,length(labels1),by=buckets)]
year_break <- year_range[seq(1,length(year_range),by=buckets)]

# determine n_alive for each age group
n_alive <- c(rep(n_indiv,length(year_range)))
DOBs1 <- unique(titre_dat[,c("individual","DOB")])[,2]
ageMask <- create_age_mask(DOBs1, strain_isolation_times)
# Create strain mask
strainMask <- create_strain_mask(titre_dat,strain_isolation_times)
masks <- data.frame(cbind(ageMask, strainMask))
# Number of people that were born before each year and have had a sample taken since that year happened
n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x) nrow(masks[masks$ageMask <=x & masks$strainMask >= x,]))    
age_ids <- unique(wow[,c("individual","Age_1")])
titre_dat_grouped <- merge(titre_dat, age_ids)
titre_dat_grouped$group <- as.numeric(titre_dat_grouped$Age_1)

# inferred attack rates
inf_chain <- all_chains$inf_chain
inf_chain <- pad_inf_chain(inf_chain)

inf_chain1 <- inf_chain[inf_chain$i%in%which(age_ids$Age_1=='<19'),]
data.table::setkey(inf_chain1, "j","samp_no","chain_no")
tmp1 <- inf_chain1[,list(V1=sum(x)),by=key(inf_chain1)]
quantiles1 <- ddply(tmp1, ~j, function(x) quantile(x$V1, c(0.025,0.1,0.5,0.9,0.975)))
quantiles1$age_group <- '<19'

inf_chain2 <- inf_chain[inf_chain$i%in%which(age_ids$Age_1=='19-64'),]
data.table::setkey(inf_chain2, "j","samp_no","chain_no")
tmp2 <- inf_chain2[,list(V1=sum(x)),by=key(inf_chain2)]
quantiles2 <- ddply(tmp2, ~j, function(x) quantile(x$V1, c(0.025,0.1,0.5,0.9,0.975)))
quantiles2$age_group <- '19-64'

inf_chain3 <- inf_chain[inf_chain$i%in%which(age_ids$Age_1=='>64'),]
data.table::setkey(inf_chain3, "j","samp_no","chain_no")
tmp3 <- inf_chain3[,list(V1=sum(x)),by=key(inf_chain3)]
quantiles3 <- ddply(tmp3, ~j, function(x) quantile(x$V1, c(0.025,0.1,0.5,0.9,0.975)))
quantiles3$age_group <- '>64'

quantiles <- rbind(quantiles1,quantiles2,quantiles3)
colnames(quantiles) <- c("j","lower","lower_50","median","upper_50","upper","age_group")
#quantiles[c("lower","median","upper")] <- quantiles[c("lower","median","upper")]/n_alive[1]
quantiles$year <- year_range[quantiles$j]
quantiles$age_group <- factor(quantiles$age_group,levels=c("<19","19-64",">64"))

n_alive_group <- get_n_alive_group(titre_dat_grouped, strain_isolation_times,melt_dat = TRUE)
n_alive_group$age_group <- levels(age_ids$Age_1)[n_alive_group$group]
n_alive_group$age_group <- factor(n_alive_group$age_group)
quantiles <- merge(quantiles, n_alive_group,by=c("j","age_group"))
quantiles[c("lower","lower_50","median","upper_50","upper")] <- quantiles[c("lower","lower_50","median","upper_50","upper")]/quantiles$n_alive

# Labels
i <- 1:(length(strain_isolation_times) + 1)
labels <- as.Date("01/01/2009", format="%d/%m/%Y") + ((i-1)*365/4 + 365/8)
labels1 <- as.Date("01/01/2009", format="%d/%m/%Y") + ((i-2)*365/4)
labels2 <- c("2009-Q1","2009-Q2","2009-Q3","2009-Q4",
             "2010-Q1","2010-Q2","2010-Q3","2010-Q4",
             "2011-Q1","2011-Q2","2011-Q3","2011-Q4")

### plot
colnames(quantiles)[which(colnames(quantiles) == "age_group")] <- "Age group:"

p_age <- ggplot(quantiles) + 
  geom_ribbon(aes(x=year, ymin=lower,ymax=upper,fill=`Age group:`),alpha=0.25) +   
  geom_ribbon(aes(x=year, ymin=lower_50,ymax=upper_50,fill=`Age group:`),alpha=0.5) +   
  geom_line(aes(x=year,y=median,colour=`Age group:`),size=0.75)+ 
  geom_point(aes(x=year,y=median,colour=`Age group:`),size=0.75)+ 
  scale_y_continuous(expand=c(0,0),limits=c(-0.005,ymax),breaks=seq(0,ymax,by=0.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = strain_isolation_times1, labels = labels2,
                     limits=c(min(strain_isolation_times-0.1),max(strain_isolation_times+0.1))) +
  scale_fill_manual(values=colour_fills_age) +
  scale_color_manual(values=colour_fills_age) +
  theme_pubr() +
  theme()+
  ylab("Estimated per capita\n incidence (per quarter)") +
  xlab("Time of infection")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.key=element_rect(color=NA),
        legend.background = element_blank(),
        legend.title = element_text(size=8,family="sans"),
        legend.text=element_text(size=8,family="sans"),
        legend.position = c(0.8, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_blank(),
        legend.direction = "horizontal",
        legend.margin = margin(6, 6, 6, 6),
        axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=margin(l=10,r=5,t=5)
  )
p_age



### Figure 3 ###
# 3a. unvaccinated AR
# 3b. vaccinated AR
# 3c. age-stratifed AR
# 3d. model parameters

p1 <- p + theme(axis.text.x=element_text(size=8),
                axis.text.y=element_text(size=8),
                axis.title.x=element_text(size=10),
                axis.title.y=element_text(size=10))
p2 <- p_age + theme(axis.text.x=element_text(size=8),
                    axis.text.y=element_text(size=8),
                    axis.title.x=element_text(size=10),
                    axis.title.y=element_text(size=10))

p3d <- ggplot(to_use[to_use$individual == 1,])+
  geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),fill="gray90")+
  geom_ribbon(aes(x=virus,ymin=lower_50, ymax=upper_50),fill="gray70")+
  geom_line(aes(x=virus, y=median))+
  geom_point(aes(x=virus, y=titre))+
  coord_cartesian(ylim=c(0,8))+
  ylab("log titre") +
  xlab("Time of virus circulation") +
  theme_bw() +
  facet_wrap(~individual)

theta_chain <- as.data.frame(all_chains$theta_chain)
theta_chain <- melt(theta_chain,id.vars=c("samp_no","chain_no"))
theta_chain$variable <- as.character(theta_chain$variable)

#ranges <- data.frame("variable"=c("mu","tau","sigma1"),"lower"=c(1.5,0,0.06),"upper"=c(3,0.1,0.13))

par_key <- c("mu"="μ","tau"="τ","sigma1"="σ",
             "wane"="ω","error"="ε","mu_short"="μ_s")
par_key <- c("mu"="μ[l]","mu_short"="μ[s]","wane"="ω","error"="ε","total_infections"="ΣZ[i]")
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



densities <- ddply(theta_chain,~variable, function(x) {
  tmp <- density(x$value)
  ret <- data.frame(tmp$x,tmp$y)
})

densities <- densities[densities$variable != "ΣZ[i]",]
sum_inf_table <- table(theta_chain[theta_chain$variable == "ΣZ[i]","value"])
densities_sum_inf <- data.frame("variable"="ΣZ[i]","tmp.x"=as.numeric(names(sum_inf_table)),"tmp.y"=as.numeric(sum_inf_table)/sum(sum_inf_table))
densities <- rbind(densities, densities_sum_inf)

probs <- c(0,0.025,0.975,1)
quantiles <- ddply(theta_chain, ~variable, function(x) quantile(x$value, probs))
factorised <- c()
median_bits_x <- NULL
median_bits_y <- NULL
for(par in par_key[par_key != "ΣZ[i]]"]){
  tmp_den <- densities[densities$variable == par,]
  factorised <- c(factorised, factor(findInterval(tmp_den$tmp.x,as.numeric(quantiles[quantiles$variable == par,2:4]))))
  tmp_den_func <- approxfun(density(theta_chain[theta_chain$variable == par,"value"]))
  tmp_median <- quantile(theta_chain[theta_chain$variable == par,"value"], 0.5)
  y_tmp <- tmp_den_func(tmp_median)
  median_bits_x[par] <- tmp_median
  median_bits_y[par] <- y_tmp
}
densities$quant <- factorised
densities[densities$variable == "ΣZ[i]","quant"] <- 1
median_segments <- data.frame(variable=par_key,x=median_bits_x,y=median_bits_y)

# 
# p3d2 <- ggplot() +
#   #geom_blank(data=ranges,aes(xmin=lower,xmax=upper)) +
#   geom_ribbon(data=densities, aes(ymin=0,ymax=tmp.y,x=tmp.x),fill="grey80",col="black") +
#   geom_ribbon(data=densities[densities$quant %in% c(3),], aes(ymin=0,ymax=tmp.y,x=tmp.x),fill="grey50",col="black") +
#   geom_linerange(data=median_segments,aes(x=x,ymin=0,ymax=y)) +
#   scale_y_continuous(expand=c(0,0,0.05,0)) +
#   ylab("Posterior density") +
#   xlab("Value") +
#   theme_classic() +
#   facet_wrap(~variable, labeller=label_parsed,scales="free",ncol=2) +
#   theme(
#     strip.background = element_blank(),
#     strip.text=element_text(hjust=1),
#     strip.text.x = element_text(size=12,family="sans",hjust=0.5)
#   )
# 
# p3 <- p3d2 + theme(axis.text.x=element_text(size=8),
#                    axis.text.y=element_text(size=8),
#                    axis.title.x=element_text(size=10),
#                    axis.title.y=element_text(size=10))
# p3


hacked_plots <- NULL
for(par in par_key){
  hacked_plots[[par]] <-  ggplot() + 
    #geom_blank(data=ranges,aes(xmin=lower,xmax=upper)) +
    geom_ribbon(data=densities[densities$variable == par,], aes(ymin=0,ymax=tmp.y,x=tmp.x),fill="grey80",col="black") +
    geom_ribbon(data=densities[densities$quant %in% c(3) & densities$variable == par,], 
                aes(ymin=0,ymax=tmp.y,x=tmp.x),fill="grey60",col="black") +
    geom_linerange(data=median_segments[median_segments$variable == par,],aes(x=x,ymin=0,ymax=y)) +
    scale_y_continuous(expand=c(0,0,0.05,0)) +
    ylab("") +
    xlab("") +
    theme_pubr() +
    facet_wrap(~variable, labeller=label_parsed,scales="free",ncol=2) +
    theme(
      strip.background = element_blank(),
      strip.text=element_text(hjust=1),
      strip.text.x = element_text(size=12,family="sans",hjust=0.5),
      axis.text.x=element_text(size=8),
      axis.text.y=element_text(size=8),
      axis.title=element_blank(),
      plot.margin=margin(l=5,t=0,b=0,r=-4)
    ) 
}

par <- "ΣZ[i]"
hacked_plots[["μ[s]"]] <- hacked_plots[["μ[s]"]] + theme(plot.margin=margin(l=-5))
hacked_plots[["ΣZ[i]"]] <-  ggplot() + 
  #geom_blank(data=ranges,aes(xmin=lower,xmax=upper)) +
  geom_bar(data=densities[densities$variable == par,], aes(y=tmp.y,x=tmp.x),fill="grey80",col="black",stat="identity") +
  geom_vline(data=median_segments[median_segments$variable == par,],aes(xintercept=x)) +
  scale_y_continuous(expand=c(0,0,0.05,0)) +
  ylab("") +
  xlab("") +
  theme_pubr() +
  facet_wrap(~variable, labeller=label_parsed,scales="free",ncol=2) +
  theme(
    axis.title=element_blank(),
    strip.background = element_blank(),
    strip.text=element_text(hjust=1),
    strip.text.x = element_text(size=12,family="sans",hjust=0.5),
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=8),
    plot.margin=margin(l=-5,t=0,b=0,r=15)
  ) 

p3 <- plot_grid(hacked_plots[[1]],hacked_plots[[2]],hacked_plots[[3]],hacked_plots[[5]],align="hv",axis="tlbr")
p3
library(grid)
library(gridExtra)
y.grob <- textGrob("Posterior density", 
                   gp=gpar(fontsize=10), rot=90)

x.grob <- textGrob("Value", 
                   gp=gpar(fontsize=10))
p3 <- grid.arrange(arrangeGrob(p3, left = y.grob, bottom = x.grob))
library(cowplot)


left_column <- plot_grid(p1,p2, labels = c('A','B'), align = 'v',axis="l", ncol = 1,label_x = 0.02)
right_column <- plot_grid(titre_pred_p,p3,labels=c("C","D"),align="v",ncol=1,rel_heights = c(1.2,1))
overall_p <- plot_grid(left_column,right_column, ncol=2, rel_widths = c(1.2,1))
overall_p


svg("Fig3_simple.svg",width=7.5,height=5.5)
plot(overall_p)
dev.off()
cairo_pdf("Fig3_simple.pdf",width=7.5,height=5.5)
plot(overall_p)
dev.off()


      
