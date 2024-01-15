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
serosolver <- FALSE

## doRNG handles seeding in foreach loop
set.seed(1234)

setwd("E:/James/Google Drive/Influenza/serosolver/methods_paper/PLOS Comp Biol/Results/case_study_2/")
## We'll be parallelising a few chains
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

filename <- "case_study_2"
filenames <- paste0(filename, "_",1:5)
sample_year <- 2009
buckets <- 1
dat_file <- system.file("extdata", "Fluscape_HI_data.csv", package = "serosolver")
raw_dat <- read.csv(dat_file,stringsAsFactors=FALSE)
raw_dat$individual <- 1:nrow(raw_dat)
melted_dat <- reshape2::melt(raw_dat,id.vars=c("individual","Age"))
colnames(melted_dat) <- c("individual","DOB","virus","titre")
melted_dat$DOB <- 2009 - melted_dat$DOB
melted_dat$virus <- as.character(melted_dat$virus)
melted_dat$virus <- as.numeric(sapply(melted_dat$virus, function(x) strsplit(x,split = "HI.H3N2.")[[1]][2]))
melted_dat[melted_dat$titre == 0,"titre"] <- 5
melted_dat$titre <- log2(melted_dat$titre/5)
melted_dat$samples <- sample_year
titre_dat <- melted_dat
titre_dat <- plyr::ddply(titre_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))

antigenic_map_file <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
antigenicMap <- read.csv(antigenic_map_file,stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- fit_dat[fit_dat$inf_times >= 1968 & fit_dat$inf_times <= sample_year,]
strain_isolation_times <- unique(fit_dat$inf_times)

par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(par_tab_path, stringsAsFactors=FALSE)
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1.3,2.7)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"fixed"] <- 1 # mu, tau, and sigma are fixed
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"values"] <- 0 # set value of mu and tau to 0
par_tab[par_tab$names == "MAX_TITRE","values"] <- 8 # set max titre to 8

mcmc_pars <- c("iterations"=500000,"target_acceptance_rate_theta"=0.44,"target_acceptance_rate_inf_hist"=0.44,"adaptive_frequency"=1000,"thin"=100,"adaptive_iterations"=100000, 
               "save_block"=1000, "thin2"=1000,"proposal_inf_hist_indiv_prop"=0.5,"proposal_ratio"=2, "burnin"=0, "proposal_inf_hist_time_prop"=1, 
               "proposal_inf_hist_distance"=3,"proposal_inf_hist_adaptive"=1,"proposal_inf_hist_indiv_swap_ratio"=0.5,"proposal_inf_hist_group_swap_ratio"=0.5,"proposal_inf_hist_group_swap_prop"=1)


if(serosolver) {
  # starting parameter values  
  res <- foreach(x = filenames) %dopar% {
    ## Not all random starting conditions return finite likelihood, so for each chain generate random
    ## conditions until we get one with a finite likelihood
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
      f <- create_posterior_func_fast(start_tab,titre_dat,fit_dat,version=2) # function in posteriors.R
      start_prob <- sum(f(start_tab$values, start_inf))
    }
    
    res <- serosolver(par_tab = start_tab, titre_dat = titre_dat,antigenic_map = fit_dat,start_inf_hist=start_inf, 
                    mcmc_pars = mcmc_pars,
                    filename = paste0("chains_main/",x), CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2, temp=1,
                    fast_version=TRUE)
  }
  beepr::beep(4)
}

## Read in the MCMC chains automatically
all_chains <- load_mcmc_chains(location="chains_main",thin=1,burnin=100000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)


## Get the MCMC chain list
list_chains <- all_chains$theta_list_chains
## Look at diagnostics for the free parameters
list_chains1 <- lapply(list_chains, function(x) x[,c("mu","sigma1","error",
                                                     "tau","total_infections","lnlike","prior_prob")])
gelman.diag(as.mcmc.list(list_chains1))
gelman.plot(as.mcmc.list(list_chains1))
effectiveSize(as.mcmc.list(list_chains1))

## Plot the MCMC trace
color_scheme_set("viridis")
mcmc_trace(list_chains1)

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
inputDat <- data.table(titre_dat)
setkey(inputDat, "individual")
age_dat <- merge(n_inf_chain, unique(inputDat[,c("individual","DOB")]))

age_dat$age <- 2009 - age_dat$DOB
age_dat$infs_per_life <- age_dat$median_infs/age_dat$age
age_dat$age_group <- cut(age_dat$age, breaks = c(0,19,40,65,90))

age_dist <- ggplot(age_dat) + 
  geom_boxplot(aes(group=age_group,y=infs_per_life*10,x=age_group),fill="grey80") +
  scale_y_continuous(breaks=0:7,lim=c(0,7),expand=c(0,0))+
  theme_pubr() +
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=8),
        plot.margin=margin(l=15,r=5))+
  ylab("Sero-responses per 10\n years alive (posterior medians)") +
  xlab("Age group")
age_dist


#############
## AR plot
#############
inf_chain <- all_chains$inf_chain

## Find samples that were in both theta and inf hist chains
chain <- as.data.frame(all_chains$theta_chain)
intersect_samps <- intersect(unique(inf_chain$samp_no), unique(chain$samp_no))
chain <- chain[chain$samp_no %in% intersect_samps,]
which_mle <- chain[which.max(chain$lnlike),c("samp_no","chain_no")]

## Take subset of chain, as do not need all samples
samps <- unique(inf_chain[,c("samp_no","chain_no")])
n_samps <- sample(1:nrow(samps), 1000)
samps <- samps[n_samps,]
samps <- rbind(samps, which_mle)
samps$samp_no1 <- 1:nrow(samps)
samps$chain_no1 <- 1
## Append the MLE estimate, note that this is max(samp_no)

inf_chain <- merge(inf_chain, samps, by=c("samp_no","chain_no"))
inf_chain <- inf_chain[,c("samp_no1","chain_no1","i","j","x")]
colnames(inf_chain)[1:2] <- c("samp_no","chain_no")
inf_chain <- pad_inf_chain(inf_chain)
## Rename columns to be more informative
colnames(inf_chain) <- c("samp_no","chain_no","individual","year","infected","group")

## Data on which strains belong to which cluster
cluster_file <- system.file("extdata", "fonville_clusters.csv", package = "serosolver")
clusters <- read.csv(cluster_file, stringsAsFactors=FALSE)
clusters <- clusters[clusters$year <= sample_year,]

## j=1 corresponds to the year 1968
inf_chain$year <- inf_chain$year + 1967

## Merge cluster data and infection history data
inf_chain <- merge(inf_chain, clusters[,c("year","cluster1")],by="year")

## Calculate ages and age groups of all individuals
titre_dat$age <- max(strain_isolation_times) - titre_dat$DOB
titre_dat$age_group <- cut(titre_dat$age,breaks=c(0,20,100),include.lowest=TRUE)
ages <- unique(titre_dat[,c("individual","age_group","DOB","age")])

## Merge infection histories with individual data
inf_chain<- merge(inf_chain, data.table(ages), by=c("individual"))

## Alive status for each individual for each time,
## only interested in individuals that were alive 
## when a virus circulated
inf_chain$alive <- inf_chain$DOB <= inf_chain$year
inf_chain <- inf_chain[inf_chain$alive,]

## Find out number of infections per year
inf_chain$potential_infection <- 1
setkey(inf_chain, "samp_no","chain_no","year")
inf_chain_ar <- inf_chain[,list(no_infected=sum(infected),
                                potential_infection=sum(potential_infection)),
                          by=key(inf_chain)]

setkey(inf_chain_ar, "year")
y <- inf_chain_ar[,list(median_ar=median(no_infected/potential_infection),
                        lower_quantile=quantile(no_infected/potential_infection,0.025),
                        upper_quantile=quantile(no_infected/potential_infection,0.975),
                        lower_quantile_50=quantile(no_infected/potential_infection,0.1),
                        upper_quantile_50=quantile(no_infected/potential_infection,0.8)),
                  by=key(inf_chain_ar)]

inf_chain_mle <- inf_chain_ar[inf_chain_ar$samp_no == max(inf_chain_ar$samp_no),]

p_ar <- ggplot(data=y) + 
  geom_ribbon(aes(x=year, ymin=lower_quantile,ymax=upper_quantile), fill="#56B4E9",alpha=0.2) +
  geom_ribbon(aes(x=year, ymin=lower_quantile_50,ymax=upper_quantile_50), fill="#56B4E9",alpha=0.5) +
  geom_line(aes(x=year,y=median_ar),col="#56B4E9") +
  geom_point(aes(x=year,y=median_ar),col="#56B4E9") +
  #geom_line(data=inf_chain_mle, aes(x=year, y=no_infected/potential_infection), linetype="dashed", size=0.8,col="forestgreen")+
  scale_x_continuous(limits=c(1967.5,sample_year+0.5),expand=c(0,0),breaks=c(1968,seq(1970,2010,by=5),2009))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  theme_pubr() +
  ylab("Annual attack rate") +
  xlab("Year of infection") +
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
        plot.margin=margin(l=10,r=5,t=5,b=10))
p_ar


##########
## AR by age by cluster
###########

## Find out number of infections per cluster for each individual
## Get age at the time of infection for each infection
inf_chain$age_at_inf <- inf_chain$year - inf_chain$DOB + 1

## Put into age groups for age at time of infection - young or old
inf_chain$age_group_at_inf <- cut(inf_chain$age_at_inf, breaks = c(0,20,100))

## Get average age at time of infection for each year
setkey(inf_chain, "samp_no","chain_no","year")
average_age <- inf_chain[,list(average_age=mean(age_at_inf)),by=key(inf_chain)]

inf_chain_infected <- inf_chain[inf_chain$infected == 1,]
setkey(inf_chain_infected, "samp_no","chain_no","year")
average_age_infected <- inf_chain_infected[,list(average_age_infected=mean(age_at_inf)),by=key(inf_chain_infected)]

all_average_ages <- merge(average_age, average_age_infected, by=c("samp_no","chain_no","year"))
all_average_ages$diff <- all_average_ages$average_age_infected - all_average_ages$average_age

average_age_quants <- ddply(all_average_ages, .(year), function(x) quantile(x$diff,c(0.025,0.1, 0.5,0.9,0.975)))
colnames(average_age_quants) <- c("year","lower","lower80","median","upper80","upper")
age_diff_plot <- ggplot(average_age_quants) +
  geom_hline(yintercept=0,linetype="dotted") +
  geom_ribbon(aes(x=year,ymin=lower,ymax=upper),alpha=0.2) +
  geom_ribbon(aes(x=year,ymin=lower80,ymax=upper80),alpha=0.5) +
  geom_line(aes(x=year,y=median)) + 
  scale_y_continuous(limits=c(-30,30)) +
  theme_pubr() +
  ylab("Average age of infected - average age") +
  xlab("Year of infection") +
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
        plot.margin=margin(l=10,r=5,t=5,b=10))
age_diff_plot

inf_chain_average <- inf_chain[,list(infected1=sum(infected), 
                                     potential_infection=sum(potential_infection)),
                               by=key(inf_chain)]

## Get total number of infections per cluster for each individual, keep age at inf
setkey(inf_chain, "individual","samp_no","chain_no","cluster1","age_group_at_inf")
inf_chain_cluster <- inf_chain[,list(infected1=sum(infected), 
                                     potential_infection=sum(potential_infection)),
                               by=key(inf_chain)]

## Convert total infections to any infections, and potential infections
## to at least one possible infection
inf_chain_cluster_once <- inf_chain_cluster
inf_chain_cluster_once[inf_chain_cluster_once$infected1 >=1,"infected1"] <- 1
inf_chain_cluster_once[inf_chain_cluster_once$potential_infection >=1,"potential_infection"] <- 1

## Get total number of at least 1 infections per cluster
setkey(inf_chain_cluster_once, "samp_no","chain_no","cluster1","age_group_at_inf")
inf_chain_cluster_once <- inf_chain_cluster_once[,list(total_infected=sum(infected1), 
                                                       total_potential_infection=sum(potential_infection)),by=key(inf_chain_cluster_once)]

## And the same for raw number of infections
setkey(inf_chain_cluster, "samp_no","chain_no","cluster1","age_group_at_inf")
inf_chain_cluster_multiple <- inf_chain_cluster[,list(total_infected=sum(infected1), 
                                                      total_potential_infection=sum(potential_infection)),
                                                by=key(inf_chain_cluster)]

## For which years was each cluster circulating?
year_ranges <- ddply(clusters, ~cluster1, function(x){
  c(min(x$year), max(x$year))
})
year_ranges$width <- year_ranges$V2 - year_ranges$V1 + 1
colnames(year_ranges) <- c("cluster1","start_year","end_year","width")

y <- inf_chain_cluster_once
y <- merge(y, year_ranges,by="cluster1")
y$cluster1 <- as.factor(y$cluster1)
colnames(y)[which(colnames(y)=="width")] <- "Years of circulation"
y$ar <- y$total_infected/y$total_potential_infection
y_summaries <- ddply(y, .(cluster1,age_group_at_inf,`Years of circulation`), 
                     function(x) quantile(x$ar, c(0.025,0.5,0.975)))

ggplot(y_summaries) + 
  geom_bar(aes(x=cluster1, fill=age_group_at_inf, y=`50%`),stat="identity",position="dodge") +
  geom_errorbar(aes(x=cluster1,group=age_group_at_inf, ymin=`2.5%`,ymax=`97.5%`),position="dodge") +
  scale_y_continuous(expand=c(0,0)) +
  theme_pubr()

p_once <- ggplot(y) +
  geom_violin(aes(x=cluster1,y=total_infected/total_potential_infection,
                  fill=`Years of circulation`, group=cluster1),
              adjust=1.2, 
              draw_quantiles=c(0.025,0.5,0.975), scale="width"
  ) + 
  scale_fill_gradient2(low="#F0E442", mid= "#56B4E9", high="#0072B2",midpoint=4,name="Number of years in circulation") +
  guides(fill=guide_colourbar(title.position="right",title.hjust=0.5,direction="vertical",label.position="left",
                              barwidth=0.5,barheight = 10, frame.colour="black",ticks=FALSE)) +
  ylab("Proportion of inviduals\n infected at least once") +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(yintercept=0,size=0.5) +
  xlab("Cluster index") +
  theme_pubr() +
  theme(
    axis.title=element_text(size=10),
    axis.text=element_text(size=8),
    strip.text=element_text(size=8,face="bold"),
    strip.background = element_blank(),
    legend.text=element_text(size=8),
    legend.position="right",
    legend.margin=margin(l=-5),
    legend.title=element_text(angle=-90,size=10),
    plot.margin=margin(l=5),
    panel.grid=element_blank()) +
  facet_wrap(~age_group_at_inf, ncol=1)
print(p_once)


#################
## Titre plot
#################

inf_chain <- all_chains$inf_chain
inf_chain1 <- inf_chain[inf_chain$chain_no == 1,]

# Plot inferred antibody titres
chain <- as.data.frame(all_chains$theta_chain)
chain1 <- chain[chain$chain_no == 1,]
rand_indivs <- 1:5
titre_dat$group <- 1
titre_preds <- get_titre_predictions(chain = chain1, infection_histories = inf_chain1, titre_dat = titre_dat, 
                                     individuals = rand_indivs,nsamp = 500,
                                     antigenic_map = fit_dat, par_tab = par_tab,expand_titredat = FALSE)
to_use <- titre_preds$predicted_observations
model_preds <- titre_preds$predictions
to_use$individual <- rand_indivs[to_use$individual]

x_labels <- seq(1970,2010,by=5)
use_indiv <- 1

titre_pred_p <- ggplot(to_use[to_use$individual == use_indiv,]) +
  geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",size=0.2)+
  geom_ribbon(data=model_preds[model_preds$individual ==use_indiv,], 
              aes(x=virus,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",size=0.2) + 
  geom_line(data=model_preds[model_preds$individual ==use_indiv,], aes(x=virus, y=median),linetype="dotted",color="grey10")+
  geom_rect(ymin=9,ymax=11,xmin=0,xmax=9000,fill="grey70")+
  geom_rect(ymin=-2,ymax=0,xmin=0,xmax=9000,fill="grey70")+
  scale_x_continuous(expand=c(0,0),breaks=x_labels) +
 # geom_vline(xintercept=seq(1970,2010,by=5),col="grey70",size=0.02,alpha=0.5) +
#  geom_hline(yintercept=0:10,col="grey70",size=0.02,alpha=0.5) +
  geom_point(data=titre_dat[titre_dat$individual ==use_indiv,], aes(x=virus, y=titre),fill="black",size=1,shape=23,col="black")+
   ylab("log titre") +
  xlab("Year of virus circulation") +
  theme_pubr()+
  theme(strip.background = element_blank(),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        legend.margin = margin(-1,-1,-3,-1),
        #legend.position="none",
        strip.text=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        plot.margin=margin(r=5,t=5,l=5))+
  coord_cartesian(ylim=c(0,10),xlim=c(1967.5,2009.5)) +
  scale_y_continuous(breaks=seq(0,10,by=2)) 
titre_pred_p

###########
## Posterior densities
###########
theta_chain <- as.data.frame(all_chains$theta_chain)
theta_chain <- melt(theta_chain,id.vars=c("samp_no","chain_no"))
theta_chain$variable <- as.character(theta_chain$variable)

#ranges <- data.frame("variable"=c("mu","tau","sigma1"),"lower"=c(1.5,0,0.06),"upper"=c(3,0.1,0.13))

par_key <- c("mu"="mu[l]","tau"="tau","sigma1"="sigma[l]",
             "error"="epsilon","total_infections"="sum(Z[i])")
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

densities <- densities[densities$variable != "sum(Z[i])",]
tmp_theta_chain <- theta_chain[theta_chain$variable == "sum(Z[i])",]
tmp_theta_chain$bucket <- cut(tmp_theta_chain$value, breaks=seq(0,2000,by=50))
sum_inf_table <- table(tmp_theta_chain[,"bucket"])
densities_sum_inf <- data.frame("variable"="sum(Z[i])","tmp.x"=as.character(names(sum_inf_table)),"tmp.y"=as.numeric(sum_inf_table)/sum(sum_inf_table))
densities_sum_inf$x_min <- seq(0,1950,by=50)
densities_sum_inf$x_max <- seq(50,2000,by=50)

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
median_segments <- data.frame(variable=par_key,x=median_bits_x,y=median_bits_y)

hacked_plots <- NULL
for(par in par_key[par_key != "sum(Z[i])"]){
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
      axis.text.x=element_text(size=6),
      axis.text.y=element_text(size=6),
      axis.title=element_blank(),
      plot.margin=margin(l=0,t=0,b=0,r=0)
    ) 
}

par <- "sum(Z[i])"
hacked_plots[["mu[s]"]] <- hacked_plots[["mu[s]"]] + theme(plot.margin=margin(l=-5))
hacked_plots[["sum(Z[i])"]] <-  ggplot() + 
  #geom_blank(data=ranges,aes(xmin=lower,xmax=upper)) +
  geom_rect(data=densities_sum_inf, aes(ymin=0,ymax=tmp.y,xmin=x_min,xmax=x_max),fill="grey80",col="black",stat="identity") +
  geom_vline(data=median_segments[median_segments$variable == par,],aes(xintercept=x)) +
  scale_y_continuous(expand=c(0,0,0.05,0)) +
  scale_x_continuous(limits=c(750,1750)) +
  ylab("") +
  xlab("") +
  theme_pubr() +
  facet_wrap(~variable, labeller=label_parsed,scales="free",ncol=2) +
  theme(
    axis.title=element_blank(),
    strip.background = element_blank(),
    strip.text=element_text(hjust=1),
    strip.text.x = element_text(size=12,family="sans",hjust=0.5),
    axis.text.x=element_text(size=6),
    axis.text.y=element_text(size=6),
    plot.margin=margin(l=-5,t=0,b=0,r=15)
  ) 

p3 <- plot_grid(hacked_plots[[1]],hacked_plots[[2]],hacked_plots[[3]],align="hv",axis="tlbr",ncol=1)
p3

y.grob <- textGrob("Posterior density", 
                   gp=gpar(fontsize=10), rot=90)

x.grob <- textGrob("Value", 
                   gp=gpar(fontsize=10))
p3 <- grid.arrange(arrangeGrob(p3, left = y.grob, bottom = x.grob))

#top_row <- plot_grid(p_ar, age_dist,rel_widths = c(1.2,1),labels=c("A","B"),align="h")
#bot_row <- plot_grid(titre_pred_p, p3, p_once,nrow=1,rel_widths=c(2,1.2,3),labels=c("C","D","E"))
top_row <- plot_grid(p_ar,labels=c("A"),label_x = 0,label_y=1)
bot_row <- plot_grid(age_dist, titre_pred_p, labels=c("B","C"),label_x = 0)
fig6 <- plot_grid(top_row,bot_row,nrow=2,rel_heights = c(1.2,1))

svg("Fig6.svg",width=6.5,height=5.5)
plot(fig6)
dev.off()
cairo_pdf("Fig6.pdf",width=6.5,height=5.5)
plot(fig6)
dev.off()
