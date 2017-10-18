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


  ## Get indexing for individuals
  indicesA <- c(0)
  individuals <- data$individual
  for(individual in unique(individuals)){
      indicesA <- c(indicesA, length(individuals[individuals==individual]))
  }
  indicesA <- cumsum(indicesA)

  ## Get indexing for samples
  
  samples <- unique(data[,c("individual","samples")])
  indicesB <- c(0)
  individuals <- samples$individual
  for(individual in unique(individuals)){
      indicesB <- c(indicesB, length(individuals[individuals == individual]))
  }
  indicesB <- cumsum(indicesB)
  
  f <- function(pars, infectionHistories){
      names(pars) <- mynames
      antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
      antigenicMapLong[antigenicMapLong < 0] <- 0
      antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
      antigenicMapShort[antigenicMapShort < 0] <- 0
                                        #return(-100000)
      return(group_likelihood_vector(pars,infectionHistories,
                                     indicesA, indicesB,
                                     samples$samples,strainIsolationTimes, 
                                     antigenicMapLong,antigenicMapShort,
                                     titres))
  }
  f
}
