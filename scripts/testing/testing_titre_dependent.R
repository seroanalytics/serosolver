setwd("~/Documents/Fluscape/serosolver")
Rcpp::compileAttributes()
devtools::load_all()
devtools::test()

theta <- c("mu"=2,"mu_short"=2.7,"tau"=0,"wane"=0,"sigma1"=0.1,"sigma2"=0.03,"boost_limit"=6,"gradient"=0.01,
           "error"=1,"wane_type"=1,"titre_dependent"=1,"kappa"=0.9,"t_change"=6)
theta1 <- theta[names(theta) != "gradient"]
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_titre.csv",
                   stringsAsFactors=FALSE)

antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, 1)
fit_dat <- fit_dat[fit_dat$inf_years >= 1968 & fit_dat$inf_years <= 2015,]
strainIsolationTimes <- unique(fit_dat$inf_years)
measuredStrains <- c(2000,2003,2007,2008,2010,2011,2013,2014)
measuredStrains <- strainIsolationTimes

fit_dat$x_coord <- seq(1,nrow(fit_dat))
fit_dat$y_coord <- seq(1,nrow(fit_dat))

antigenicMapMelted <- c(outputdmatrix.fromcoord(fit_dat[,c("x_coord","y_coord")]))


antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])

#infectionHistory <- sample(c(0,0,0,1),length(strainIsolationTimes),replace=TRUE)
infectionHistory <- sample(c(0,1),length(measuredStrains),replace=TRUE,prob = c(0.9,0.1))

samplingTime <- 2016
measurementMapIndices <- match(measuredStrains, strainIsolationTimes) - 1
infectionMapIndices <- seq_along(strainIsolationTimes) - 1
numberStrains <- length(strainIsolationTimes)
t1 <- 2000
t2 <- 2005
t3 <- 2015
filled <- which(infectionHistory == 1)
theta["wane"] <- 0.1
theta["wane_type"] <- 1

boosting_vec_indices <- rep(seq(1,4),each=12)-1
additional_arguments <- list("boosting_vec_indices"=boosting_vec_indices,"mus"=c(1,3,5,5))
y <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                      infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                      antigenicMapShort,numberStrains,DOB=0,additional_arguments)


y <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                               infectionMapIndices, 2002, measurementMapIndices, antigenicMapLong,
                               antigenicMapShort,numberStrains)
y1 <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                           infectionMapIndices, 2003, measurementMapIndices, antigenicMapLong,
                           antigenicMapShort,numberStrains)
y2 <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                            infectionMapIndices, 2010, measurementMapIndices, antigenicMapLong,
                            antigenicMapShort,numberStrains)
plot(y,type='l',ylim=c(0,8))
lines(y1,col="red")
lines(y2,col="blue")


y <- infection_model_indiv_many(theta, infectionHistory, strainIsolationTimes,
                           infectionMapIndices, t1, measurementMapIndices, antigenicMapLong,
                           antigenicMapShort,numberStrains,1000)

#y1 <- infection_model_indiv(theta1, infectionHistory, strainIsolationTimes,
#                            infectionMapIndices, samplingTime, measurementMapIndices, antigenicMapLong,
#                            antigenicMapShort,numberStrains)
par(mfrow=c(3,1))
plot(y~strainIsolationTimes,ylim=c(0,8),type='l', main=t1)
abline(v=strainIsolationTimes[which(infectionHistory==1)])
y1 <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                           infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                           antigenicMapShort,numberStrains)
y_old <- infection_model_indiv_OLD(theta1, infectionHistory, strainIsolationTimes,
                                   infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                                   antigenicMapShort,numberStrains)
identical(y1, y_old)

plot(y1~strainIsolationTimes,ylim=c(0,8),type='l',main=t2)
#points(y1~strainIsolationTimes,col="blue")
abline(v=strainIsolationTimes[which(infectionHistory==1)])

y2 <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                            infectionMapIndices, t3, measurementMapIndices, antigenicMapLong,
                            antigenicMapShort,numberStrains)
plot(y2~strainIsolationTimes,ylim=c(0,8),type='l',main=t3)
#points(y1~strainIsolationTimes,col="blue")
abline(v=strainIsolationTimes[which(infectionHistory==1)])


library(microbenchmark)
infectionHistory <- rep(0, length(infectionHistory))
#infectionHistory[20] <- 1

infectionHistory <- sample(c(0,1),length(measuredStrains),replace=TRUE,prob = c(0.5,0.5))
theta1 <- theta
theta1["titre_dependent"] <- 0
theta["wane_type"] <- 0
theta2 <- c("mu"=1.8,"mu_short"=2,"tau"=0,"wane"=1,"sigma1"=0.08,"sigma2"=0.004,"error"=1,"wane_type"=0)

filled = which(infectionHistory==1)
theta["gradient"] <- 1
infectionHistory <- sample(c(0,1),length(measuredStrains),replace=TRUE,prob = c(0.9,0.1))
y1 <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                            infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                            antigenicMapShort,numberStrains)
y_old <- infection_model_indiv_OLD(theta, infectionHistory, strainIsolationTimes,
                                   infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                                   antigenicMapShort,numberStrains)
identical(y1, y_old)


filled = which(infectionHistory==1)
infectionHistory <- sample(c(0,1),length(measuredStrains),replace=TRUE,prob = c(0.9,0.1))

new_f <- function(){
  y <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                        infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                        antigenicMapShort,numberStrains)
}

old_f <- function(){
  y <- infection_model_indiv_OLD(theta, infectionHistory, strainIsolationTimes,
                            infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                            antigenicMapShort,numberStrains)
}
microbenchmark(new_f,old_f,times=100000)
theta["titre_dependent"] <- 1
theta["gradient"] <- 0.1
microbenchmark(infection_model_indiv_many(theta, infectionHistory, strainIsolationTimes,
                                     infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                                     antigenicMapShort,numberStrains,n=2000),
               infection_model_indiv_OLD_many(theta, infectionHistory, strainIsolationTimes,
                                         infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                                         antigenicMapShort,numberStrains,2000),
               times=500)

y <- infection_model_indiv(theta, infectionHistory, strainIsolationTimes,
                           infectionMapIndices, t2, measurementMapIndices, antigenicMapLong,
                           antigenicMapShort,numberStrains)
