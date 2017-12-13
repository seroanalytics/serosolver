create_post_func_dopar <- function(parTab, data, PRIOR_FUNC,
                             ...){
  
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  
  strains <- unique(antigenicMap$inf_years)
  strainIndices <- match(strains, strains) - 1
  antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))
  
  ## Get unique measurement sets for each individual at
  ## each sampling time for each repeat
  samples <- unique(data[,c("individual","samples","run")])
  samples <- samples[order(samples$individual, samples$run, samples$samples),]
  
  individuals <- unique(samples$individual)
  chunks <- length(individuals)/ncore
  indices <- c(0,cumsum(rep(chunks,ncore)))
  
  indiv_list <- lapply(1:ncore, function(i) individuals[(indices[i]+1):indices[i+1]])
  data_list <- lapply(1:ncore, function(i) data[data$individual %in% indiv_list[[i]],])
  samples_list <- lapply(1:ncore, function(i) samples[samples$individual %in% indiv_list[[i]],])
  
  titres <- NULL
  sampleTimes <- NULL
  virusIndices <- NULL
  strainIsolationTimes <- NULL
  indicesData <- NULL
  indicesDataOverall <- NULL
  indicesSamples <- NULL
  
  for(j in 1:ncore){
    tmpIndividuals <- indiv_list[[j]]
    tmpData <- data_list[[j]]
    tmpSamples <- samples_list[[j]]
    ## Extract vector of sampling times and individual indices for speed
    sampleTimes[[j]] <- tmpSamples$samples
    individuals <- tmpSamples$individual
    
    ## 
    titres[[j]] <- tmpData$titre
    virusIndices[[j]] <- match(tmpData$virus, antigenicMap$inf_years) - 1

    ## Note that the strain isolution times in the data set should have
    ## corresponding indices in the antigenic map table
    strainIsolationTimes[[j]] <- tmpData$virus
    
    indicesDataTmp <- NULL
    for(i in 1:nrow(tmpSamples)){
      indicesDataTmp <- c(indicesDataTmp, nrow(tmpSamples[tmpData$individual == tmpSamples[i,"individual"] &
                                                         tmpData$samples == tmpSamples[i,"samples"] &
                                                         tmpData$run == tmpSamples[i,"run"],]))
    }
   indicesData[[j]] <- indicesDataTmp
    ## Get indexing for individuals. This works out which rows in the titre data correspond
    ## to which individuals
    indicesSamplesTmp <- c(0)
    for(individual in unique(individuals)){
      indicesSamplesTmp <- c(indicesSamplesTmp, length(individuals[individuals==individual]))
    }
    indicesSamplesTmp <- cumsum(indicesSamplesTmp)
    indicesSamples[[j]] <- indicesSamplesTmp
    
    indicesDataOverallTmp <- NULL
    for(individual in unique(individuals)){
      indicesDataOverallTmp <- c(indicesDataOverallTmp, nrow(data[data$individual == individual,]))
    }
    indicesDataOverallTmp <- cumsum(c(0,indicesDataOverallTmp))
    indicesDataOverall[[j]] <- indicesDataOverallTmp
  }
  
  r_likelihood <- function(expected, data, theta){
    largeI <- data > theta["MAX_TITRE"]
    smallI <- data <= 0
    restI <- data > 0 & data <= theta["MAX_TITRE"]
    
    large <- pnorm(theta["MAX_TITRE"], expected[largeI],theta["error"],lower.tail=FALSE,log.p=TRUE)
    small <- pnorm(1, expected[smallI],theta["error"],lower.tail=TRUE,log.p=TRUE)
    rest <- log(pnorm(data[restI]+1,expected[restI],theta["error"],lower.tail=TRUE,log.p=FALSE) - 
                  pnorm(data[restI],expected[restI], theta["error"],lower.tail=TRUE,log.p=FALSE))
    return(sum(large, small, rest))
  }
    
   ## The function pointer
  f <- function(pars){
    names(pars) <- mynames
    
    ## Work out short and long term boosting cross reactivity
    antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
    antigenicMapLong[antigenicMapLong < 0] <- 0
    antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
    antigenicMapShort[antigenicMapShort < 0] <- 0
    
    ## Now pass to the C++ function
    res <- foreach(i=1:ncore,.combine='sum',.export=c("titre_data_group")) %do% {
      y <- titre_data_group(pars, infectionHistories[(indices[i]+1):indices[i+1],], strains, strainIndices, sampleTimes[[i]],
                            indicesData[[i]],indicesDataOverall[[i]],indicesSamples[[i]], virusIndices[[i]], 
                            antigenicMapLong, antigenicMapShort)
      r_likelihood(y, titres[[i]], pars)
    }
  }
  f
}

