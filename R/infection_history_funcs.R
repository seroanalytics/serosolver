#' @export
setup_infection_histories<- function(dat, strainIsolationTimes, ageMask){
    SAMPLE_PROB <- 0.2
    n_indiv <- length(unique(dat$individual))
    n_strain <- length(strainIsolationTimes)
    samplingTimes <- unique(dat$sample)
    infectionHistories <- matrix(0,nrow=n_indiv,ncol=n_strain)
    
    index <- 1
    ## For each individual
    tmpInfHistIndiv <- matrix(0,nrow=length(samplingTimes),ncol=n_strain)
    for(indiv in unique(dat$individual)){
        ## For each sampling time
        index1 <- 1
        for(sampleTime in unique(dat$sample)){
            ## For each strain
            tmpInfHist <- numeric(n_strain)
            index2 <- 1
            for(strain in unique(dat$strain)){
                tmpTitre <- dat[dat$strain == strain &
                                dat$sample == sampleTime &
                                dat$individual == indiv,"titre"]
                tmpInf <- 0
                if(!is.na(tmpTitre) && (tmpTitre >= 4 & runif(1) > SAMPLE_PROB)){
                    tmpInf <- 1
                }
                tmpInfHist[index2] <- tmpInf
                index2 <- index2 + 1
            }
            tmpInfHistIndiv[index1,] <- tmpInfHist
            index1 <- index1 + 1
        }

        infectionHistories[index,] <- as.numeric(colSums(tmpInfHistIndiv) > 0)
        ## Add infection at some point in the last 10 years
        forcedInfection <- which(strainIsolationTimes==sample(strainIsolationTimes[strainIsolationTimes <= sampleTime & strainIsolationTimes >= sampleTime-10],1))
        if(sum(infectionHistories[index, which(strainIsolationTimes < sampleTime)])==0) infectionHistories[index, forcedInfection] <- 1
        if(ageMask[index] > 1) infectionHistories[index, 1:ageMask[index]] <- 0
        index <- index + 1
    }
    infectionHistories

}
