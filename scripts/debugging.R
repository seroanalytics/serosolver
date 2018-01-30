library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/serosolver_own/serosolver")
devtools::load_all()


fit_dat <- read.csv("~/Documents/serosolver_own/serosolver/data/antigenicMap_AK.csv")

parTab <- read.csv("~/Documents/serosolver_own/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(2,12)

samplingTimes <- 2009:2015
strainIsolationTimes <- unique(fit_dat$inf_years)
dat <- simulate_data(parTab, 1, 2, 1,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0.1, 10,75,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=TRUE)
data <- dat[[1]]
#data <- data[data$virus %in% seq(1968,2014,by=2),]
data <- data[complete.cases(data),]
infectionHistories <- dat[[2]]
ages <- dat[[3]]


pars <- parTab$values
mynames <- parTab$names
names(pars) <- parTab$names

antigenicMap <- fit_dat

## Isolate data table as vectors for speed
titres <- data$titre
## The entry of each virus in the data corresponding to the antigenic map
virusIndices <- match(data$virus, antigenicMap$inf_years) - 1
## Unique strains that an individual could be exposed to
strains <- unique(antigenicMap$inf_years)
## The entry of each strain in the antigenic map table
strainIndices <- match(strains, strains) - 1

## Note that the strain isolation times in the data set should have
## corresponding indices in the antigenic map table
strainIsolationTimes <- data$virus
antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))


## Get unique measurement sets for each individual at
## each sampling time for each repeat
samples <- unique(data[,c("individual","samples","run")])
samples <- samples[order(samples$individual, samples$run, samples$samples),]

## Extract vector of sampling times and individual indices for speed
sampleTimes <- samples$samples
individuals <- samples$individual  
n_indiv <- length(unique(samples$individual))

indicesData <- NULL
for(i in 1:nrow(samples)){
    indicesData <- c(indicesData, nrow(samples[data$individual == samples[i,"individual"] &
                                               data$samples == samples[i,"samples"] &
                                               data$run == samples[i,"run"],]))
}

## Get indexing for individuals. This works out which rows in the titre data correspond
## to which individuals
indicesSamples <- c(0)
for(individual in unique(individuals)){
    indicesSamples <- c(indicesSamples, length(individuals[individuals==individual]))
}
indicesSamples <- cumsum(indicesSamples)

indicesDataOverall <- NULL
for(individual in unique(individuals)){
    indicesDataOverall <- c(indicesDataOverall, nrow(data[data$individual == individual,]))
}
indicesDataOverall <- cumsum(c(0,indicesDataOverall))  

indicesOverallDiff <- diff(indicesDataOverall)


names(pars) <- mynames
## Work out short and long term boosting cross reactivity
antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
antigenicMapShort[antigenicMapShort < 0] <- 0

y <- titre_data_group(pars, infectionHistories, strains, strainIndices, sampleTimes,
                      indicesData,indicesDataOverall,indicesSamples, virusIndices, 
                      antigenicMapLong, antigenicMapShort)
data$y <- y
ggplot(data) +
    geom_line(aes(x=virus,y=y),col="red")+
    geom_point(aes(x=virus,y=titre)) +
    facet_grid(individual~samples) +
    theme_bw()
    
