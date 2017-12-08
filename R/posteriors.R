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
create_post_func <- function(parTab, data,
                             antigenicMap, PRIOR_FUNC,
                             ...){
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  sampleTimes <- unique(data$samples)
  titres <- data$titre
  strainIsolationTimes <- unique(antigenicMap$inf_years)
  antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))

  ## Get indexing for individuals. This works out which rows in the titre data correspond
  ## to which individuals
  indicesA <- c(0)
  individuals <- data$individual
  for(individual in unique(individuals)){
      indicesA <- c(indicesA, length(individuals[individuals==individual]))
  }
  indicesA <- cumsum(indicesA)

  ## Get indexing for samples
  ## This finds which samples correspond to which individual
  samples <- unique(data[,c("individual","samples")])
  indicesB <- c(0)
  individuals <- samples$individual
  for(individual in unique(individuals)){
      indicesB <- c(indicesB, length(individuals[individuals == individual]))
  }
  indicesB <- cumsum(indicesB)

  ## The function pointer
  f <- function(pars, infectionHistories){
      names(pars) <- mynames

      ## Work out short and long term boosting cross reactivity
      antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
      antigenicMapLong[antigenicMapLong < 0] <- 0
      antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
      antigenicMapShort[antigenicMapShort < 0] <- 0

      ## Now pass to the C++ function
      return(group_likelihood_vector(pars,infectionHistories,
                                     indicesA, indicesB,
                                     samples$samples,strainIsolationTimes, 
                                     antigenicMapLong,antigenicMapShort,
                                     titres))
  }
  f
}
