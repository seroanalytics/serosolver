#' Posterior function pointer
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to passs the data set back and forth.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param PRIOR_FUNC user function of prior for model parameters. Should take parameter values only
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infectionHistories as unnamed arguments
#' @export
create_post_func <- function(parTab, data, antigenicMap,
                             PRIOR_FUNC,
                             ...){
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  
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
  
  r_likelihood <- function(expected, data, theta){
    liks <- numeric(length(expected))
    largeI <- data > theta["MAX_TITRE"]
    smallI <- data <= 0
    restI <- data > 0 & data <= theta["MAX_TITRE"]
    
    liks[largeI] <- pnorm(theta["MAX_TITRE"], expected[largeI],theta["error"],lower.tail=FALSE,log.p=TRUE)
    liks[smallI] <- pnorm(1, expected[smallI],theta["error"],lower.tail=TRUE,log.p=TRUE)
    liks[restI] <- log(pnorm(data[restI]+1,expected[restI],theta["error"],lower.tail=TRUE,log.p=FALSE) - 
                  pnorm(data[restI],expected[restI], theta["error"],lower.tail=TRUE,log.p=FALSE))
    return(liks)
  }
  
  indicesOverallDiff <- diff(indicesDataOverall)
  ## The function pointer
  f <- function(pars, infectionHistories){
      names(pars) <- mynames
      ## Work out short and long term boosting cross reactivity
      antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
      antigenicMapLong[antigenicMapLong < 0] <- 0
      antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
      antigenicMapShort[antigenicMapShort < 0] <- 0

      ## Now pass to the C++ function
      y <- titre_data_group(pars, infectionHistories, strains, strainIndices, sampleTimes,
                                 indicesData,indicesDataOverall,indicesSamples, virusIndices, 
                                 antigenicMapLong, antigenicMapShort)
      liks <- r_likelihood(y, titres, pars)
                                        #return(list(liks, y, titres))
                                        #liks <- numeric(length(y))
                                        #return(liks)
                                        #liks <- dnorm(titres,y,sd=pars["error"],1)
                                        #prior <- dunif(rowSums(infectionHistories),0,ncol(infectionHistories),log=TRUE)
     # prior <- dbinom(rowSums(infectionHistories), ncol(infectionHistories), 0.5,  1)
      #return(rep(0,n_indiv))
      return(sum_buckets(liks, indicesOverallDiff))
  }
  f
}


#' Posterior function pointer
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to passs the data set back and forth.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param PRIOR_FUNC user function of prior for model parameters. Should take parameter values only
#' @param infectionHistories other arguments to pass to the posterior solving function
#' @param ... blah
#' @return a single function pointer that takes only pars and infectionHistories as unnamed arguments
#' @export
create_post_func1 <- function(parTab, data,antigenicMap,
                             PRIOR_FUNC,infectionHistories, ...){
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  
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
  
  r_likelihood <- function(expected, data, theta){
    liks <- numeric(length(expected))
    largeI <- data > theta["MAX_TITRE"]
    smallI <- data <= 0
    restI <- data > 0 & data <= theta["MAX_TITRE"]
    
    liks[largeI] <- pnorm(theta["MAX_TITRE"], expected[largeI],theta["error"],lower.tail=FALSE,log.p=TRUE)
    liks[smallI] <- pnorm(1, expected[smallI],theta["error"],lower.tail=TRUE,log.p=TRUE)
    liks[restI] <- log(pnorm(data[restI]+1,expected[restI],theta["error"],lower.tail=TRUE,log.p=FALSE) - 
                         pnorm(data[restI],expected[restI], theta["error"],lower.tail=TRUE,log.p=FALSE))
    return(liks)
  }
  
  indicesOverallDiff <- diff(indicesDataOverall)
  ## The function pointer
  f <- function(pars){
    names(pars) <- mynames
    ## Work out short and long term boosting cross reactivity
    antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
    antigenicMapLong[antigenicMapLong < 0] <- 0
    antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
    antigenicMapShort[antigenicMapShort < 0] <- 0

    ## Now pass to the C++ function
    y <- titre_data_group(pars, infectionHistories, strains, strainIndices, sampleTimes,
                          indicesData,indicesDataOverall,indicesSamples, virusIndices, 
                          antigenicMapLong, antigenicMapShort)
    liks <- r_likelihood(y, titres, pars)
  
    #liks <- dnorm(titres,y,sd=pars["error"],1)
    lik <- -sum(sum_buckets(liks, indicesOverallDiff))
    
    lik
  }
  f
}
