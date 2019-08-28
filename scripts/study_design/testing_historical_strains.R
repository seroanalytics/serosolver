code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

## Specify paramters controlling the MCMC procedure
mcmc_pars <- c("iterations"=500000,"adaptive_period"=500000, "burnin"=50000,
               "hist_sample_prob"=0.5,"thin"=50,"thin_hist"=200,
               "swap_propn"=0.5,"hist_switch_prob"=0.2,"year_swap_propn"=0.5)

## Antigenic map for cross reactivity parameters
antigenic_map <- read.csv(file.path(code.dir,"data/fonville_map_approx.csv"),
                         stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, buckets=1)

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968,]
strain_isolation_times <- unique(fit_dat$inf_years)

## Read in parameter table
par_tab <- read.csv(file.path(code.dir,"inputs/par_tab_study_design.csv"))


  
#Run the MCMC 
for(filename in filename_vec){
  ## Read in ages
  ages <- read.csv(paste("data/study_design/",filename,"_8_ages.csv",sep=""), stringsAsFactors=FALSE)
  for(h in 1:8){
    titre_dat <- read.csv(paste('data/study_design/', filename,"_", h,"_dat.csv", sep=""), stringsAsFactors=FALSE)
    
    titre_dat <- merge(titre_dat, ages)
    
    tryCatch(res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat, 
                    antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                    mvr_pars = NULL, start_inf_hist = NULL, filename=paste('chains/',filename,"_",h,sep=""),
                    CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = NULL,
                    version = 2,  
                    fast_version = TRUE),
             error = function(e) write.table(paste(filename,"_",h,sep=""), paste(filename,"_",h,".txt",sep="")))
  }
  
}



# Check chain dimensions
for(filename in filename_vec){
  ages <- read.csv(paste("data/study_design/",filename,"_8_ages.csv",sep=""), stringsAsFactors=FALSE)
  p <- list()
  for(h in 1:8){
    chain <- read.csv(paste("~/Documents/GitHub/serosolver/chains/",filename,"_",h,"_chain.csv",sep=""))
    print(paste(dim(chain), filename, "_", h))
  } 
}


filename <- c('all_ages_historic')

# Plot attack rates 
#for(filename in filename_vec){
  ages <- read.csv(paste("data/study_design/",filename,"_8_ages.csv",sep=""), stringsAsFactors=FALSE)
  p <- list()
  for(h in 1:8){
    titre_dat <- read.csv(paste("data/study_design/", filename,"_", h,"_dat.csv", sep=""), stringsAsFactors=FALSE)
    titre_dat <- merge(titre_dat, ages)
    inf_chain <- data.table::fread(paste("chains/", filename,"_", h,"_infection_histories.csv",sep=""))
    if( h==8) {
      inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"] + 500000),]
    }else{
      inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
      
    }
    
    p[[h]] <- plot_attack_rates(infection_histories = inf_chain, titre_dat, strain_isolation_times, colour_by_taken = FALSE, by_val = 10) +
          geom_hline(yintercept=0.15, linetype="dashed", color = "gray")
  } 
  
  p_all <- do.call(grid.arrange,c(p, nrow =2))
  
  ggsave(paste(filename,".png",sep=""),p_all, width = 20, height = 6)
#}


# Plot attack rates 
h_vec <- 1:8
  ages <- read.csv(paste("data/study_design/",filename,"_8_ages.csv",sep=""), stringsAsFactors=FALSE)
  p <- list()
  for(k in 1:4){
    h <- h_vec[k]
    titre_dat <- read.csv(paste("data/study_design/", filename,"_", h,"_dat.csv", sep=""), stringsAsFactors=FALSE)
    titre_dat <- merge(titre_dat, ages)
    inf_chain <- data.table::fread(paste("chains/", filename,"_", h,"_infection_histories.csv",sep=""))
    if( h==8) {
      inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"] + 500000),]
    }else{
      inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
      
    }
    
    p[[k]] <- plot_attack_rates(infection_histories = inf_chain, titre_dat, strain_isolation_times, colour_by_taken = FALSE, by_val = 10) +
      geom_hline(yintercept=0.15, linetype="dashed", color = "gray")
  } 
  
  p_all <- do.call(grid.arrange,c(p, nrow =2))
  
  ggsave(paste(filename,".png",sep=""),p_all, width = 20, height = 6)


# Plot likelihood
for(filename in filename_vec){
  png(paste(filename,"_likelihood.png"),width=5000,height=2000,res=300,units='px')
  par(mfrow=c(2,4))
  for(h in 1:8){
    chain <- read.csv(paste("~/Documents/GitHub/serosolver/chains/",filename,"_",h,"_chain.csv",sep=""))
    
    if( h==8) {
      chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"] + 500000),]
    }else{
      chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
    }
    try(plot(chain$likelihood,type='l',main=paste("Scenario ",h),ylab='likelihood'))
  }
  dev.off()
}





png("boxplot.png",width=5000,height=2000,res=300,units='px')
filename <- "all_ages_historic"
par(mfrow=c(2,4))
for(h in 1:8){
  chain1 <-  read.csv(paste("~/Documents/GitHub/serosolver/chains/",filename,"_",h,"_chain.csv",sep=""))
  chain1 <- chain1[chain1$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
  chain1 <- chain1[,c('mu', 'mu_short','wane','sigma1','sigma2','error')]
  
  
  which_pars<-c('mu', 'mu_short','wane','sigma1','sigma2','error')
  true.val<-par_tab$values[par_tab$names%in%which_pars]
  mat<-matrix(nrow=dim(chain1)[1],ncol=6)
  for(j in 1:6){
    chain.adj<-(chain1[,j]-true.val[j])/true.val[j]
    mat[,j]<-chain.adj
  }
  
  #boxplot of parameters
  names<-c(expression(mu[1]),expression(mu[2]),expression(omega),expression(sigma[1]),expression(sigma[2]),expression(epsilon))
  boxplot(mat,ylab="(Estimated - True) / True",col='gray',main=paste("Scenario ",h),ylim=c(min(mat,mat2),max(mat,mat2)),names=names)
  abline(h=0,col="gray")
}
dev.off()
