#' Posterior function pointer
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to passs the data set back and forth.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param PRIOR_FUNC user function of prior for model parameters. Should take parameter values only
#' @param version integer specifying which version of the posterior function to use. This is mainly used for turning off or on the likelihood for debugging the implicit prior. The versions are as follows: 1) Likelihood turned on and no explicit prior. Any prior is implicit from the proposal distribution; 2) No likelihood and no explicit prior; 3) Explicit prior term (as specified by PRIOR_FUNC) and no likelihood. Any implicit prior from the proposal will also be included; 4) Version of the posterior distribution where we infer an explicit force of infection for each epoch; 5+) Explicit likelihood and prior - any implicit prior from the proposal will also be included, so proposals should be symmetric
#' @param ageMask see \code{\link{create_age_mask}} - a vector with one entry for each individual specifying the first epoch of circulation in which an individual could have been exposed
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infectionHistories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @export
create_post_func <- function(parTab, data, antigenicMap,
                             PRIOR_FUNC,version=1,ageMask=NULL,
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
  
  indicesOverallDiff <- diff(indicesDataOverall)
  if(version==1){
      print("likelihood only - prior implicit in proposal")
      ## The function pointer
 
      f <- function(pars, infectionHistories){

          names(pars) <- mynames
          ## Work out short and long term boosting cross reactivity
          #antigenicMapLong <- 1-pars1["sigma1"]*antigenicMapMelted
          #antigenicMapLong[antigenicMapLong < 0] <- 0
          #antigenicMapShort <- 1-pars1["sigma2"]*antigenicMapMelted
          #antigenicMapShort[antigenicMapShort < 0] <- 0

          antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, pars1["sigma1"])
          antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, pars1["sigma2"])
          
          ## Now pass to the C++ function
          y <- titre_data_group(pars, infectionHistories, strains, strainIndices, sampleTimes,
                                indicesData,indicesDataOverall,indicesSamples, virusIndices, 
                                antigenicMapLong, antigenicMapShort)
          liks <- r_likelihood(y, titres, pars)                                        
          ##return(rep(0,n_indiv))
          return(sum_buckets(liks, indicesOverallDiff))
      }
  } else if(version==2){
      print("nothing - posterior implicit in proposal")
      theta_indices <- which(parTab$identity == 1)
      lambda_indices <- which(parTab$identity == 2)
      mynames <- parTab[theta_indices,"names"]
      f <- function(pars1, infectionHistories){
          pars <- pars1[theta_indices]
          lambdas <- pars1[lambda_indices]
          names(pars) <- mynames
          liks <- rep(0,n_indiv)
          liks <- liks + calc_lambda_probs_indiv(lambdas, infectionHistories, ageMask)
          return(liks)
      }
  } else if(version == 3){
      print("explicit prior - no likelihood. proposals should be symmetric")
      f <- function(pars, infectionHistories){
          names(pars) <- mynames
          return(PRIOR_FUNC(pars, infectionHistories))
      }
  } else if(version==4) {
      print("Explicit FOI inference, lambda. Using likelihood and no explicit prior")
      theta_indices <- which(parTab$identity == 1)
      lambda_indices <- which(parTab$identity == 2)
      mynames <- parTab[theta_indices,"names"]
      ## The function pointer
      f <- function(pars1, infectionHistories){
          pars <- pars1[theta_indices]
          lambdas <- pars1[lambda_indices]
        
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
          liks <- sum_buckets(liks, indicesOverallDiff)
          liks <- liks + calc_lambda_probs_indiv(lambdas, infectionHistories, ageMask)
          ##return(rep(0,n_indiv))
                                        #return(sum_buckets(liks, indicesOverallDiff))
          return(liks)
      }
  } else {
      print("likelihood and prior - proposals should be symmetric")
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
          liks <- sum_buckets(liks, indicesOverallDiff)
          liks <- liks + PRIOR_FUNC(pars, infectionHistories)         
          return(liks)
      }
  }
  f
}


#' Posterior function pointer - for use with optim
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to passs the data set back and forth.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param PRIOR_FUNC user function of prior for model parameters. Should take parameter values only
#' @param infectionHistories infection histories to solve model with
#' @param ... other arguments to pass to the posterior solving function
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
      lik <- -sum(sum_buckets(liks, indicesOverallDiff))
      lik
  }
  f
}

#' @export
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

#' Model solving function
#'
#' As in \code{\link{create_post_func}}, but returns the matrix of titres rather than a vector of likelihoods
#' @export
create_model_func <- function(parTab, data, antigenicMap,
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
    return(y)
  }
}

#' Calculate FOI log probability
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param lambdas a vector of FOIs
#' @param infHist the matrix of infection histories
#' @param ageMask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a single log probability
#' @export
calc_lambda_probs <- function(lambdas, infHist, ageMask){
    lik <- 0
    for(i in 1:ncol(infHist)){
        lik <- lik + sum(log((lambdas[i]^infHist[which(ageMask <= i),i] * (1-lambdas[i])^(1-infHist[which(ageMask <= i),i]))))
    }
    lik
}

#' Calculate FOI log probability vector
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param lambdas a vector of FOIs
#' @param infHist the matrix of infection histories
#' @param ageMask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a vector of log probabilities for each individual
#' @export
calc_lambda_probs_indiv <- function(lambdas, infHist, ageMask){
  lik <- numeric(nrow(infHist))
  for(i in 1:ncol(infHist)){
    lik <- lik + log(((lambdas[i]^infHist[,i]) * (1-lambdas[i])^(1-infHist[,i])))* as.numeric(ageMask <= i)
  }
  lik
}



#' FOR DEBUGGING
#'
#' Brute force implementation of calculating the explicit FOI
#' @export
calc_lambda_probs_indiv_brute<- function(lambdas, infHist, ageMask){
  lik <- numeric(nrow(infHist))
  for(j in 1:nrow(infHist)){
      lik[j] <- 0
      age <- ageMask[j]
      for(i in 1:ncol(infHist)){
          if(i >= age){
              lik[j] <- lik[j] + log(lambdas[i]^infHist[j,i] * (1-lambdas[i])^(1-infHist[j,i]))
          }
      }
  }
  lik
}
