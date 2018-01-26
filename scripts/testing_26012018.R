library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
buckets <- 1
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
pars <- parTab$values
names(pars) <- parTab$names

infectionHistory <- c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 
                      0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 
                      0, 0, 0, 0, 1, 0, 0)


antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, 1)

strainIsolationTimes <- unique(fit_dat$inf_years)
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1) 
antigenicMapMelted <- c(outputdmatrix.fromcoord(fit_dat[,c("x_coord","y_coord")]))

antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
antigenicMapShort[antigenicMapShort < 0] <- 0

circulationTimes <- strainIsolationTimes
circulationMapIndices <- match(circulationTimes, circulationTimes)-1

indices <- match(strainIsolationTimes,strainIsolationTimes)-1

sampleTimes <- seq(1970,2010,by=10)
measurementTimes <- circulationTimes
measuredTimes <- rep(circulationTimes, length(sampleTimes))
dataIndices <- rep(length(circulationTimes),length(sampleTimes))
measurementIndices <- match(measuredTimes, circulationTimes)-1


#measurementTimes <- sample(circulationTimes, 9)
#measurementTimes <- measurementTimes[order(measurementTimes)]
#dataIndices <- rep(length(measurementTimes), length(sampleTimes))
#measuredTimes <- rep(measurementTimes, length(sampleTimes))
#measurementIndices <- match(measuredTimes, circulationTimes)-1

y1 <- titre_data_individual(pars, infectionHistory, circulationTimes, circulationMapIndices,
                           sampleTimes,dataIndices,measurementIndices,
                           antigenicMapLong,antigenicMapShort,length(circulationTimes)
                           )
wow <- rep(sampleTimes, each=length(measurementTimes))
omg1 <- data.frame("sample"=wow,"virus"=measuredTimes,"titre"=y1)
ggplot(omg) + 
  geom_point(aes(x=virus,y=titre)) + 
  geom_line(data=omg1,aes(x=virus,y=titre),col="red")+
  facet_wrap(~sample) + theme_bw()
