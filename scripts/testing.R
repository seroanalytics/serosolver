sourceCpp("~/Documents/Fluscape/serology-model/src/infection_model.cpp")
source("~/Documents/Fluscape/serology-model/R/simulate_data.R")
source("~/Documents/Fluscape/serology-model/R/infection_history_funcs.R")

group <- 1
n_indiv <- 10
parTab <- read.csv("~/Documents/Fluscape/serology-model/inputs/text_partab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serology-model/data/antigenic_map.csv")
strainIsolationTimes <- unique(antigenicMap$inf_years)
samplingTimes <- unique(antigenicMap$inf_years)

dat <- simulate_data(parTab,group,n_indiv,strainIsolationTimes, samplingTimes,antigenicMap, 0.5,0.1,5,80)
View(dat[[3]])






create_post_func <- function(parTab, data, samples, strainIsolationTimes, antigenicMap, PRIOR_FUNC){
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  present <- samples$present
  sampleTimes <- unique(samples$sample)
  titres <- data$titre
  
  antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap))
  
  f <- function(pars, infectionHistories){
    names(pars) <- mynames
    return(group_likelihood(pars,infectionHistories, present, sampleTimes,strainIsolationTimes, 
                            1-pars["sigma1"]*antigenicMap,1-pars["sigma2"]*antigenicMap, titres))
  }
  f
}



n_strains <- length(pInf)
n_indiv <- length(ages)
indivs <- 1:n_indiv
infectionHistories <- matrix(0,ncol=n_strains,nrow=n_indiv)
attackRates <- rlnorm(length(pInf), meanlog=log(pInf)- infSD^2/2,sdlog=infSD)create_post_func <- function(parTab, data, infectionHistories, samples, strainIsolationTimes, antigenicMap1, antigenicMap2, PRIOR_FUNC){
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  present <- samples$present
  sampleTimes <- unique(samples$sample)
  titres <- data$titre
  f <- function(pars){
    names(pars) <- mynames
    return(group_likelihood(pars,infectionHistories, present, sampleTimes,strainIsolationTimes, 
                            1-(antigenicMap1*pars["sigma1"]),1-(antigenicMap2*pars["sigma2"]), titres))
  }
  f
}

mcmcPars <- c(iterations=10000,adaptive_period=5000,opt_freq=1000,save_block=100,burnin=0,popt=0.44,thin=1)
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab$fixed[i] == 0){
    startTab$values[i] <- runif(1,startTab$lower_bound[i],startTab$upper_bound[i])
  }
}
wow <- antibodyKinetics::run_MCMC(startTab, dat, mcmcPars, "test",create_post_func,NULL,NULL,0.1, 
                                  infectionHistories=infectionHistories, samples=samples, 
                                  strainIsolationTimes=strainIsolationTimes,antigenicMap1=c(antigenicMap),
                                  antigenicMap2=c(antigenicMap))

chain <- read.csv(wow$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))
