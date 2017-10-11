# - - - - - - - - - - - - - - - -
# Set arbitrary initial condition for infection history - aim is to ensure likelihood is valid
setuphistIC<-function(ii,jj,inf.n,test.list,testyear_index, test_years, inf_years){ # ii=participant | jj=test year
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  spyear=unique(as.numeric(test.jj[3,])) # year of samples taken
  
  hist0=rep(0,inf.n)   
  #hist0[sample(c(1:inf.n),round(0.1*inf.n))]=1
  
  # Check test data available - may be issue if age column added too
  if(length(test.jj[,1])>1){
    
    # Set up test strains
    titredat=as.numeric(test.jj[2,]) # Define titre data
    
    # Use simple cutoff for titres -- set high titres = 1 in history if >=4
    for(i in 1:length(spyear)){
      if(max(titredat[(as.numeric(test.jj[3,])==spyear[i])])>=4 & runif(1)>0.2 ){
        hist0[(inf_years==spyear[i])]=1
      }
    }

  }
  
  min.range = max(1,testyear_index[1]-10) # Add an infection within past 10 years (useful for initial sampling)
  inf_index = inf_years-min(inf_years)+1
  inf_pick = sample(c(1:inf.n)[inf_index>=min.range & inf_index<= testyear_index[1]],1)  # pick strain within plausible region to add
  if(sum(hist0[inf_years < min(test_years)])==0){hist0[inf_pick]=1} # Make sure at least one infection occurs before earliest test year
  hist0
  
}

setupInfectionHistory <- function(dat, strainIsolationTimes, ageMask){
    SAMPLE_PROB <- 0.2
    n_indiv <- length(unique(dat$individual))
    n_strain <- length(strainIsolationTimes)
    samplingTimes <- unique(dat$sample)
    
    infectionHistories <- matrix(0,nrow=n_indiv,ncol=n_strain)
    index <- 1
    ## For each individual
    for(indiv in unique(dat$individual)){
        ## For each sampling time
        tmpInfHistIndiv <- matrix(0,nrow=length(samplingTimes),ncol=n_strain)
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
        infectionHistories[index, 1:ageMask[index]] <- 0
        index <- index + 1
    }
    infectionHistories
}
