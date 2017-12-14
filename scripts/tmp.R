devtools::load_all()

parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map.csv")
antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))
pars <- parTab$values
mynames <- parTab$names
names(pars) <- parTab$names

## Work out short and long term boosting cross reactivity
antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
antigenicMapShort[antigenicMapShort < 0] <- 0

samps <- seq(1970*12,2010*12,by=1)
strains <- unique(antigenicMap$inf_years)

n_indiv <- 1000
pars <- parTab$values
names(pars) <- parTab$names
data <- simulate_data(parTab,1,n_indiv, strains, samps, nsamps=2,antigenicMap,0,0,5,80)
infectionHistories <- data$infectionHistories
y <- simulate_group(n_indiv,pars, infectionHistories,strains,samps, 2, antigenicMapLong,antigenicMapShort)
titreDat <- data$data
testedStrains <- seq(1970*2,2010*2,by=5)*6
#testedStrains <- sample(strains, 200)
#titreDat <- y
titreDat <- titreDat[titreDat$virus %in% testedStrains,]
titreDat$run <- 1
#titreDat$titre <- floor(titreDat$titre)
mcmcPars <- c("iterations"=10000,"popt"=0.44,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=5000,"save_block"=100)
f <- create_post_func1(parTab,titreDat,NULL,antigenicMap=antigenicMap,infectionHistories=infectionHistories)
liks <- NULL
for(i in 1:100){
  pars["mu"] <- (i-1)/10
  liks[i] <- f(pars)
}
plot(liks)
which.max(liks)


#parTab[parTab$names == "error","fixed"] <- 1
startTab <- parTab
for(i in 1:nrow(startTab)) startTab[i,"values"] <- runif(1,startTab[i,"lower_bound"],startTab[i,"upper_bound"])
startTab[startTab$names == "MAX_TITRE","values"] <- 8
#startTab[startTab$names == "error","values"] <- 1
#startTab[startTab$names == "error","fixed"] <- 1

res <- lazymcmc::run_MCMC(startTab, data=titreDat, mcmcPars,"test",create_post_func1, mvrPars=NULL,PRIOR_FUNC=NULL,antigenicMap=antigenicMap,infectionHistories=infectionHistories)
chain <- read.csv(res$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))

covMat <- cov(chain[chain$sampno > mcmcPars["adaptive_period"],2:(ncol(chain)-1)])
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
bestPars <- get_best_pars(chain)
startTab$values <- bestPars
res1 <- lazymcmc::run_MCMC(startTab, data=titreDat, mcmcPars,"test1",create_post_func, mvrPars=mvrPars,PRIOR_FUNC=NULL,antigenicMap=antigenicMap,infectionHistories=infectionHistories)
