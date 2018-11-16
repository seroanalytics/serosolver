pars <- c("mu"=2,"mu_short"=2,"tau"=0.05,"wane"=1,"sigma1"=0.1,"sigma2"=0.03,"MAX_TITRE"=8,"error"=1)

buckets <- 4
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
antigenicMapMelted <- c(outputdmatrix.fromcoord(fit_dat[,c("x_coord","y_coord")]))
antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
antigenicMapShort[antigenicMapShort < 0] <- 0

strainIsolationTimes <- unique(fit_dat$inf_years)
strainIndices <- match(strainIsolationTimes, strainIsolationTimes) - 1
nYears <- length(strainIsolationTimes)
infectionHistory <- sample(c(0,1),nYears,prob=c(0.98,0.02),replace=TRUE)

samplingTimes <- sample(strainIsolationTimes, 2)
samplingTimes <- order(samplingTimes)

sampledTitres <- sample(strainIsolationTimes, 10)
sampledTitres <- sampledTitres[order(sampledTitres)]
sampledTitres <- c(sampledTitres, sampledTitres)
dataIndices <- c(10,10)

measuredMapIndices <- match(sampledTitres, fit_dat$inf_year) - 1

microbenchmark::microbenchmark(
y <- titre_data_individual(pars, infectionHistory, strainIsolationTimes, strainIndices, samplingTimes,
                           dataIndices, measuredMapIndices, antigenicMapLong, antigenicMapShort, length(strainIsolationTimes))
)
