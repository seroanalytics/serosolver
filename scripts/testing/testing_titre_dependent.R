setwd("~/Documents/Fluscape/serosolver")
#devtools::load_all()

theta <- c("mu"=1.8,"mu_short"=2,"tau"=0.04,"wane"=1,"sigma1"=0.1,"sigma2"=0.03,"boost_limit"=3,"gradient"=0.2,
           "error"=1)

parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",
                   stringsAsFactors=FALSE)

antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, 1)
fit_dat <- fit_dat[fit_dat$inf_years >= 2000 & fit_dat$inf_years <= 2015,]
strainIsolationTimes <- unique(fit_dat$inf_years)
measuredStrains <- c(2000,2003,2007,2008,2010,2011,2013,2014)
measuredStrains <- strainIsolationTimes

antigenicMapMelted <- c(outputdmatrix.fromcoord(fit_dat[,c("x_coord","y_coord")]))
antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])

#infectionHistory <- sample(c(0,0,0,1),length(strainIsolationTimes),replace=TRUE)
infectionHistory <- c(0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0)
samplingTime <- 2016
measurementMapIndices <- match(measuredStrains, strainIsolationTimes) - 1
infectionMapIndices <- seq_along(strainIsolationTimes) - 1
numberStrains <- length(strainIsolationTimes)

y <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                               infectionMapIndices, samplingTime, measurementMapIndices, antigenicMapLong,
                               antigenicMapShort,numberStrains)
plot(y,ylim=c(0,8))
