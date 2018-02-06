#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' Multivariate proposal function
#'
#' Function used to give multivariate normal proposals for free model parameters.
#' Takes into account parameter covariance and ensures Containment condition with beta, if covMat0 (the iedntity matrix) is specified.
#' @param covMat the 2D covariance matrix for all of the parameters
#' @param covMat0 optional, usually the identity matrix for theta
#' @param useLog flag. If TRUE, propose on log scale
#' @param beta Beta as in Rosenthal and Roberts 2009
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
#' @useDynLib serosolver
mvr_proposal <- function(values, fixed, covMat, covMat0 = NULL, useLog=FALSE, beta=0.05){
    proposed <- values
    # Sample either from a single covariance matrix or weighted average of the identity matrix and
    # given cov matrix, if specified. On either a log or linear scale.
    if(is.null(covMat0)){
        if(!useLog) proposed[fixed] <- MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat)
        else proposed[fixed] <- exp(MASS::mvrnorm(n=1,mu=log(proposed[fixed]),Sigma=(5.6644/length(fixed))*covMat))
    } else {
        if(!useLog){
            proposed[fixed] <-
                (1-beta)*MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat) +
                beta*MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(0.01/length(fixed))*covMat0)
        } else {
            proposed[fixed] <-
                (1-beta)*exp(MASS::mvrnorm(n=1,mu=log(proposed[fixed]),Sigma=(5.6644/length(fixed))*covMat)) +
                beta*exp(MASS::mvrnorm(n=1,mu=log(proposed[fixed]),Sigma=(0.01/length(fixed))*covMat0))
        }
    }
    return(proposed)
}


#' Brute force infection history proposal
#'
#' Performs a flipping/swapping infection history update step for a matrix of infection histories. 50/50 chance of performing a flip or a swap
#' @param newInfHist a matrix of infection histories - rows for individuals, columns for infection epochs. Contents should be 1s and 0s
#' @param sampledIndivs a vector of indices describing rows in the infection history matrix that should be updated
#' @param ageMask a vector (one value for each individual) giving the first infection epoch that an individual could have been exposed in. That is, if an individual was born in the 7th epoch, their entry in ageMask would be 7.
#' @param moveSizes when performing a move step, how far should two epochs be swapped?
#' @param alpha deprecated from beta binomial proposal - kept for compatibility for now
#' @param beta deprecated from beta binomial proposal - kept for compatibility for now
#' @return a matrix of infection histories matching the input newInfHist
#' @export
infection_history_betabinom_symmetric<- function(newInfHist, sampledIndivs, ageMask, moveSizes, alpha=1, beta=1){
    newInf <- newInfHist
    ## For each individual
    for(indiv in sampledIndivs){
        ## Isolate infection history
        x <- newInfHist[indiv, ageMask[indiv]:ncol(newInfHist)]
        maxI <- length(x)
        rand1 <- runif(1)

        ## Flip or swap with prob 50%
        if(rand1 < 1/2){
            ## Choose a location and turn 0 -> 1 or 1 -> 0
            loc <- sample(length(x), 1)            
            x[loc] <- !x[loc]                                       
        } else {
            ## Choose a location
            id1 <- sample(length(x),1)
            moveMax <- moveSizes[indiv]

            ## Choose another location up to moveMax epochs away
            move <- sample(-moveMax:moveMax,1)
            id2 <- id1 + move

            ## Control boundaries
            if(id2 < 1) id2 <- maxI + id2
            if(id2 > maxI) id2 <- id2 - maxI

            ## Swap the contents of these locations
            tmp <- x[id1]
            x[id1] <- x[id2]
            x[id2] <- tmp       
        }
        newInf[indiv,ageMask[indiv]:ncol(newInfHist)]=x
    }
    return(newInf)    
}

#' Beta binomial infection history update
#'
#' Samples a new infection history from a beta binomial distribution for the specified number of individuals. Note that only one epoch is updated each time with either a flip or a swap step.
#' @param newInfHist a matrix of infection histories - rows for individuals, columns for infection epochs. Contents should be 1s and 0s
#' @param sampledIndivs a vector of indices describing rows in the infection history matrix that should be updated
#' @param ageMask a vector (one value for each individual) giving the first infection epoch that an individual could have been exposed in. That is, if an individual was born in the 7th epoch, their entry in ageMask would be 7.
#' @param moveSizes when performing a move step, how far should two epochs be swapped?
#' @param alpha alpha parameter of beta binomial
#' @param beta beta parameter of beta binomial
#' @return a matrix of infection histories matching the input newInfHist
#' @export
infection_history_betabinom <- function(newInfHist, sampledIndivs, ageMask, moveSizes, alpha, beta){
    newInf <- newInfHist
    ## For each individual
    for(indiv in sampledIndivs){
        x <- newInfHist[indiv, ageMask[indiv]:ncol(newInfHist)]
        maxI <- length(x)
        rand1 <- runif(1)

        ## With prob 0.5 swap or move/add
        if(rand1 < 0.5){

            ## Choose a location
            loc <- sample(length(x), 1)            
            x_new <- x_old <- x

            ## This location can either be a 1 or a 0
            x_new[loc] <- 1
            x_old[loc] <- 0
            ##probA <- dbb(sum(x_new), length(x), alpha, beta)/choose(length(x), sum(x_new))
            ##probB <- dbb(sum(x_old), length(x), alpha, beta)/choose(length(x), sum(x_old))
            ##ratio <- probA/(probA + probB)

            ## Add a 1 some proportion of the time depending on how many 1s are already present
            ## (the sum(x[-loc]) gives the number of 1s in the vector less the location to be
            ## updated
            prob1 <- (alpha+sum(x[-loc]))/(alpha+beta+(length(x)-1))
            if(runif(1)<prob1){
                x[loc] <- 1
                ## Otherwise, make a 0
            } else {
                x[loc] <- 0
            }
        } else {
            ## Choose a location
            id1 <- sample(length(x),1)
            moveMax <- moveSizes[indiv]

            ## Propose a location up to moveMax epochs away
            move <- sample(-moveMax:moveMax,1)
            id2 <- id1 + move

            ## Control boundary conditions
            if(id2 < 1) id2 <- maxI + id2
            if(id2 > maxI) id2 <- id2 - maxI

            ## Swap contents
            tmp <- x[id1]
            x[id1] <- x[id2]
            x[id2] <- tmp       
        }
        newInf[indiv,ageMask[indiv]:ncol(newInfHist)]=x
    }
    return(newInf)    
}


#' DEPRECATED - implemented in Cpp for speed
#' @export
infection_history_betabinom_group <- function(newInfHist, sampledIndivs, ageMask, moveSizes, nInfs, alpha, beta){
    newInf <- newInfHist
    for(indiv in sampledIndivs){
        x <- newInf[indiv,ageMask[indiv]:ncol(newInfHist)]
        maxI <- length(x)
        if(runif(1) < 1/2){
            
            k <- nInfs[indiv]
            locs <- sample(length(x), k)
            number_1s <- sum(x[-locs])
            n <- length(x[-locs])

            for(i in 1:k){
                ratio <- (alpha+number_1s)/(alpha+beta+n)
                if(runif(1) < ratio){
                    x[locs[i]] <- 1
                    number_1s <- number_1s + 1
                } else {
                    x[locs[i]] <- 0
                }
                n <- n + 1
            }
        } else {
            for(i in 1:nInfs[indiv]){
                id1 <- sample(length(x),1)
                moveMax <- moveSizes[indiv]
                move <- sample(-moveMax:moveMax,1)
                id2 <- id1 + move
                if(id2 < 1) id2 <- maxI + id2
                if(id2 > maxI) id2 <- id2 - maxI
                tmp <- x[id1]
                x[id1] <- x[id2]
                x[id2] <- tmp
            }
        }
        newInf[indiv,ageMask[indiv]:ncol(newInfHist)]=x        
    }
    return(newInf)    
}


#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size. Random walk may be on a linear or log scale
#' @param values a vector of the parameters to be explored
#' @param lower_bounds a vector of the low allowable bounds for the proposal
#' @param upper_bounds a vector of the upper allowable bounds for the proposal
#' @param steps a vector of step sizes for the proposal
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
#' @useDynLib serosolver
univ_proposal <- function(values, lower_bounds, upper_bounds,steps, index){
    mn <- lower_bounds[index]
    mx <- upper_bounds[index]

    rtn <- values
    
    x <- toUnitScale(values[index],mn,mx)

    ## 5th index is step size
    stp <- steps[index]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

    ## Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

    ## Cyclical boundary conditions
    ##if (x < 0) x <- 1 + x	
    ##if (x > 1) x <- x - 1
    
    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx)
    rtn
}

#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
#' @export
#' @useDynLib serosolver
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    step <- max(0.00001, step)
    return(step)
}

#' @export
ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  # Flat priors on theta => symmetric update probability
  calc.lik = exp(marg_likelihood_star-marg_likelihood)
  calc.lik[calc.lik>1]=1 
  calc.lik
}

#' Propose initial infection histories - OLD VERSION
#'
#' Given a matrix of titre data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param dat the matrix of titres data with columns for individual, sample, and titre
#' @param strainIsolationTimes vector of real times for all strains
#' @param ageMask vector of indices for each individual corresponding to the first index of the strainIsolationTimes vector that each individual could be infected with
#' @param sample_prob given an infection seems likely based on titre, suggest infection with some probability
#' @param titre_cutoff specifies how high the titre must be to imply an infection
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of potential infecting strains
setup_infection_histories_OLD <- function(dat, strainIsolationTimes, ageMask,sample_prob, titre_cutoff=3){
    # 
    SAMPLE_PROB <- sample_prob
    n_indiv <- length(unique(dat$individual))
    n_strain <- length(strainIsolationTimes)
    samplingTimes <- unique(dat$samples)
    sampleTime <- max(samplingTimes)
    infectionHistories <- matrix(0,nrow=n_indiv,ncol=n_strain)
    
    index <- 1
    ## For each individual
    tmpInfHistIndiv <- matrix(0,nrow=length(samplingTimes),ncol=n_strain)
    for(indiv in unique(dat$individual)){
        ## For each sampling time
        tmpInfHist <- numeric(n_strain)
        index2 <- 1
        for(strain in strainIsolationTimes){
            tmpTitre <- max(dat[dat$virus == strain &
                                dat$individual == indiv,"titre"], 0)
            tmpInf <- 0
                                        # If high titre, set infection presence to 1 with some probability (0.2)
            if(!is.na(tmpTitre) && (tmpTitre >= titre_cutoff & runif(1) > SAMPLE_PROB)){
                tmpInf <- 1
            }
            tmpInfHist[index2] <- tmpInf
            index2 <- index2 + 1
        }
        infectionHistories[index,] <- tmpInfHist
        ## Add infection at some point in the last 10 years
        forcedInfection <- which(strainIsolationTimes==sample(strainIsolationTimes[strainIsolationTimes <= sampleTime & strainIsolationTimes >= (sampleTime-10)],1))
        ## Pick strain within plausible region to add
        if(sum(infectionHistories[index, which(strainIsolationTimes < sampleTime)])==0) infectionHistories[index, forcedInfection] <- 1
        ## Make sure that no infections happened before birth
        if(ageMask[index] > 1) infectionHistories[index, 1:ageMask[index]] <- 0
        index <- index + 1
    }
    infectionHistories

}

#' Propose initial infection histories
#'
#' Given a matrix of titre data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' The idea is to move along time in the context of antigenic drift and look at an individual's titre against each strain. Where titres are raised, we suggest an infection. However, to avoid suggesting multiple infections for regions of high antigenic similarity, we place a necessary gap (defined by `space`) between proposed infection times.
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param data the matrix of titres data with columns for individual, sample, and titre
#' @param ages a matrix of ages for each individual. Columns should be 1) the individual's ID; 2) the individual's date of birth (called DOB) in the same time resolution of the model (ie. 1991 for years, 1991*12 for months etc)
#' @param strainIsolationTimes vector of real times for all strains
#' @param space how many epochs must separate proposed infections
#' @param titre_cutoff specifies how high the titre must be to imply an infection
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of potential infecting strains
#' @export
setup_infection_histories_new <- function(data, ages, strainIsolationTimes, space=5, titre_cutoff=2){
  startInf <- NULL
  individuals <- unique(data$individual)

  ## For each individual
  for(individual in individuals){
      ## Isolate that individual's data and date of birth
      dat <- data[data$individual == individual,]
      dob <- as.numeric(ages[ages$individual == individual,"DOB"])
      strains <- unique(dat$virus)

      ## Only look at strains that circulated when an individual was alive
      strains <- strains[strains >= dob]
      infYears <- NULL
      i <- 0
      while(i < length(strains)){ ## Start going through each strain
          i <- i + 1
          strain <- strains[i] ## Get current strain of interest
          measurement <- -1
          dist <- 0
          
          titre <- max(dat[dat$virus == strain, "titre"]) ## Get max titre against this strain
          if(titre >= titre_cutoff){ ## If elevated against this strain, assume an infection
              newInf <- strain
              ## Begin counting up distance
              while(dist < space & i < length(strains)){
                  i <- i + 1
                  strain <- strains[i]
                  dist <- strain - newInf ## How many years since last infection?
                  new_titre <- max(dat[dat$virus == strain, "titre"]) ## Get max titre against next strain along
                  ## If this one is better, replace and reset distance
                  if(new_titre > titre){
                      newInf <- strain
                      titre <- new_titre
                      dist <- 0
                  }
              }
              infYears <- c(infYears, newInf)
              dist <- 0
              measurement <- -1
          }
      }
      infections <- rep(0, length(strainIsolationTimes))
      infections[match(infYears, strainIsolationTimes)] <- 1
      startInf <- rbind(startInf, infections)
  }
  colnames(startInf) <- strainIsolationTimes
  rownames(startInf) <- NULL
  return(startInf)
}

#' @export
create_age_mask <- function(ages, strainIsolationTimes, n_indiv){
    ageMask <- sapply(ages$DOB, function(x){
        if(is.na(x)){
            1
        } else {
            which(as.numeric(x <= strainIsolationTimes) > 0)[1]
        }
    })
    return(ageMask)
}