strains <- unique(antigenicMap$inf_years)
numberStrains <- length(strains)
strainIsolationTimes <- strains
strainIsolationTimes <- sample(strains,15)
strainIsolationTimes <- strainIsolationTimes[order(strainIsolationTimes)]
virusIndices <- match(strainIsolationTimes,strains)-1
samplingTime <- 2010
theta <- pars

infectionHistory <- sample(c(0,0,0,1),numberStrains,replace=T)
infectionHistory <- rep(0,numberStrains)
infectionHistory[which(strains==1968)] <- 1
infectionHistory[which(strains==1990)] <- 1
infectionHistory[which(strains==2009)] <- 1


y <- infection_model_indiv(pars, infectionHistory, samplingTime, strainIsolationTimes, virusIndices,
                      antigenicMapLong,antigenicMapShort,numberStrains)

samplingTimes <- c(2007)
dataIndices <- c(48)
a <- sample(strains,21)
b <- sample(strains,18)
c <- sample(strains,17)
a <- a[order(a)]
b <- b[order(b)]
c <- c[order(c)]

strainIsolationTimes <- c(a,b,c)
strainIsolationTimes <- strains
virusIndices <- match(strainIsolationTimes,strains)-1
