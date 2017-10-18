setwd("~/Documents/Fluscape/serosolver")


group <- 1
n_indiv <- 10
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/text_partab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/antigenic_map.csv")
strainIsolationTimes <- unique(antigenicMap$inf_years)
samplingTimes <- seq(max(strainIsolationTimes)-5,max(strainIsolationTimes),by=1)

dat <- simulate_data(parTab,group,n_indiv,strainIsolationTimes, samplingTimes,antigenicMap, 0.1,0.1,5,80)
samples <- dat[[2]]
infectionHistories <- dat[[3]]
data <- dat[[1]]
ages <- dat[[4]]
ages <- data.frame(individuals=1:n_indiv,ages=ages)
mcmcPars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=100000,
              "save_block"=100,"thin2"=5,"histSampleProb"=0.1,"switch_sample"=2)

startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab$fixed[i] == 0){
    startTab$values[i] <- runif(1,startTab$lower_bound[i],startTab$upper_bound[i])
  }
}

f <- create_post_func(parTab, finalDat,  fit_dat, NULL)
sum(f(parTab$values, infectionHistories))

covMat <- diag(rep(1,nrow(parTab)))
mvrPars <- list(covMat, 2.38/sqrt(nrow(parTab[parTab$fixed == 0,])),w=0.01)
#mvrPars <- list(covMat, 0.03, 0.8)
devtools::load_all()

subsetIndivs <- sample(unique(finalDat$indiv),100)
tmpDat <- finalDat[finalDat$indiv %in% subsetIndivs,]
tmpSamples <- samples[samples$indiv %in% subsetIndivs,]

output <- run_MCMC(startTab, tmpDat, mcmcPars, "test1",create_post_func, NULL, NULL, 0.2, fit_dat, ages)
chain <- read.csv(output$chain_file)
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]

bestPars1 <- zikaProj::get_best_pars(chain)
plot(coda::as.mcmc(chain))


chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat, 2.38/sqrt(nrow(parTab[parTab$fixed == 0,])),w=0.8)
startTab1 <- startTab
startTab1$values <- bestPars

mcmcPars1 <- mcmcPars
mcmcPars1["popt"] <- 0.234
output1 <- run_MCMC(startTab1, data, mcmcPars, "test1",create_post_func, mvrPars, NULL, 0.2, antigenicMap, samples, 
                   NULL)
chain1 <- read.csv(output1$chain_file)
chain1 <- chain1[chain1$sampno >= mcmcPars["adaptive_period"],c(which(parTab$fixed == 0) + 1,ncol(chain))]
plot(coda::as.mcmc(chain1))



infectionHistories1 <- data.table::fread("test1_infectionHistories.csv")
infectionHistories <- infectionHistories1[infectionHistories1$sampno %in% 
                                             unique(infectionHistories1$sampno)[seq(1,length(unique(infectionHistories1$sampno)),by=100)],]
#infectionHistories <- infectionHistories[infectionHistories$sampno == chain[which.max(chain$lnlike),"sampno"]-1,1:47]
i <- 1
tmp <- infectionHistories1[infectionHistories1$individual == i,1:47]
infCounts <- apply(tmp[,1:47],2,function(x) table(factor(x, levels=c(0,1))))
infSD <- apply(tmp[,1:47],2,sd)
y <- (infCounts[2,]/colSums(infCounts))
#print(infSD)
x <- infectionHistories[i,]
wow <- data.frame(real=x,estimated=y)
View(wow)
plot(wow$estimated~seq(1968,2014,by=1))

#startTab$values <- as.numeric(chain[which.max(chain$lnlike),2:(ncol(chain)-1)])
#mcmcPars["popt"] <- 0.234
#output1 <- run_MCMC(parTab, data, mcmcPars, "test",create_post_func, mvrPars, NULL, 0.2, antigenicMap, samples, NULL)

