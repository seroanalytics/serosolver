#' @export
create_post_func <- function(parTab, data, samples,
                             antigenicMap, PRIOR_FUNC,
                             ...){
  pars1 <- parTab$values
  mynames <- parTab$names
  names(pars1) <- parTab$names
  present <- samples$present
  sampleTimes <- unique(samples$samples)
  titres <- data$titre
  strainIsolationTimes <- unique(antigenicMap$inf_years)
  
  antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))
  
  f <- function(pars, infectionHistories){
      names(pars) <- mynames
      antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
      antigenicMapLong[antigenicMapLong < 0] <- 0
      antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
      antigenicMapShort[antigenicMapShort < 0] <- 0
      #return(-100000)
      return(group_likelihood_vector(pars,infectionHistories, present, sampleTimes,strainIsolationTimes, 
                              antigenicMapLong,antigenicMapShort, titres))
  }
  f
}
