library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

n_indiv <-100
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)

## 1
buckets <- 1
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets,]
antigenicMapMelted1 <- c(outputdmatrix.fromcoord(fit_dat[,c("x_coord","y_coord")]))
strainIsolationTimes <- unique(fit_dat$inf_years)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=FALSE)


titreDat <- dat[[1]]
titreDat1 <- titreDat[titreDat$virus %in% (viruses*buckets),]
infectionHistories <- infHist <- dat[[2]]
f1 <- create_post_func(parTab,titreDat1,fit_dat, NULL, 1)
parTab1 <- parTab
infHist1 <- infHist
f1(parTab1$values, infHist1)


##2 
buckets <- 12
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets,]

antigenicMapMelted2 <- c(outputdmatrix.fromcoord(fit_dat[,c("x_coord","y_coord")]))
strainIsolationTimes <- unique(fit_dat$inf_years)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=TRUE)


titreDat <- dat[[1]]
titreDat2 <- titreDat[titreDat$virus %in% (viruses*buckets),]
infectionHistories <- infHist <- dat[[2]]
f2 <- create_post_func(parTab,titreDat2,fit_dat)
parTab2 <- parTab
infHist2 <- infHist

microbenchmark::microbenchmark(f1(parTab1$values, infHist1),f2(parTab2$values, infHist2))

Rprof(tmp1 <- tempfile())
for(i in 1:1000) f1(parTab1$values, infHist1)
Rprof()
annual <- summaryRprof(tmp1)

Rprof(tmp2 <- tempfile())
for(i in 1:1000) f2(parTab2$values, infHist2)
Rprof()
monthly <- summaryRprof(tmp2)


wow <- microbenchmark::microbenchmark(f1(parTab1$values, infHist1), 
                                      f4(parTab4$values, infHist4),
                                      f5(parTab5$values, infHist5),
                                      f6(parTab6$values, infHist6),
                               f2(parTab2$values,infHist2),
                               f7(parTab7$values, infHist7),
                               f3(parTab3$values,infHist3),
                               f9(parTab9$values,infHist9),
                               f8(parTab8$values,infHist8),
                               f10(parTab10$values,infHist10),
                               f11(parTab11$values,infHist11))
wow1 <- summary(wow)
lq <- wow1$lq
uq <- wow1$uq
mq <- wow1$median

minT <- min(mq)
lq <- lq/minT
uq <- uq/minT
mq <- mq/minT
x <- c(1,2,3,4,6,8,12,18,24,48,96)
plot(lq~x,type='l',col="red",ylim=c(0,100))
lines(uq~x,col="red")
lines(mq~x,col="black")
y <- x^2/34 + 1
plot(mq~x)
lines(y~x)


Rcpp::sourceCpp("scripts/testing/quicker_vector.cpp")

pars <- parTab1$values
names(pars) <- parTab1$names
sigma <- pars["sigma1"]
system.time(
for(i in 1:1000){
  #antigenicMapLong2 <- max(sigma*antigenicMapMelted2,0)
  antigenicMapLong2a <- f_vector(antigenicMapMelted2, pars["sigma1"])
  antigenicMapShort2a <- f_vector(antigenicMapMelted2, pars["sigma2"])
  #antigenicMapLong2[antigenicMapLong2 < 0] <- 0
}
)
system.time(
  for(i in 1:1000){
    antigenicMapLong2b <- 1-pars["sigma1"]*antigenicMapMelted2
    antigenicMapLong2b[antigenicMapLong2b < 0] <- 0
    antigenicMapShortb <- 1-pars["sigma2"]*antigenicMapMelted
    antigenicMapShortb[antigenicMapShortb < 0] <- 0
  }
)

