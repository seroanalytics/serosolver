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

#' Infection history proposal group
#'
#' Proposes new infection histories for a vector of infection histories, where rows represent individuals and columns represent years. Proposals are either removal, addition or switching of infections.
#' Also requires the indices of sampled individuals, the vector of strain isolation times, and a vector of age masks (ie. which index of the strainIsolationTimes vector is the first year in which
#' an individual *could* be infected).
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param newInfectionHistories an n*m matrix of 1s & 0s indicating infection histories, where n is individuals and m i strains
#' @param sampledIndivs the indices of sampled individuals to receive proposals
#' @param ageMask the vector of indices for each individual specifiying which index of strainIsolationTimes is the first strain each individual could have seen
#' @param nInfs the number of infections to move/add/remove
#' @return a new matrix matching newInfectionHistories in dimensions with proposed moves
#' @export
infection_history_proposal_group <-function(newInfectionHistories,sampledIndivs,ageMask,nInfs=1){
    newInf <- newInfectionHistories
    for(indiv in sampledIndivs){ # Resample subset of individuals
        rand1 = runif(1)
        x=newInfectionHistories[indiv,ageMask[indiv]:ncol(newInfectionHistories)] # Only resample years individual was alive
        
        if(rand1<1/3){
            infectID= which(x>0)
            n <- min(nInfs, length(infectID))
            if(n>0){
                x[sample(infectID,n)]=0 # Why double? DEBUG
            }
        }
        ## Add infection
        if(rand1>1/3 & rand1<2/3){
            ninfecID=which(x==0)
            n <- min(nInfs, length(ninfecID))
            if(n>0){
                x[sample(ninfecID,n)]=1
            }
        }
        ## Move infection position
        if(rand1>2/3){
            infectID=which(x > 0)
            ninfecID=which(x == 0)
            n <- min(nInfs, length(infectID), length(ninfecID))
            if(n){
                x[sample(infectID,n)]=0
                x[sample(ninfecID,n)]=1
            }
        }
        newInf[indiv,ageMask[indiv]:ncol(newInfectionHistories)]=x # Only =1 if individual was alive
    } # end loop over individuals
    return(newInf)
}
#' Infection history proposal
#'
#' Proposes new infection histories for a vector of infection histories, where rows represent individuals and columns represent years. Proposals are either removal, addition or switching of infections.
#' Also requires the indices of sampled individuals, the vector of strain isolation times, and a vector of age masks (ie. which index of the strainIsolationTimes vector is the first year in which
#' an individual *could* be infected).
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param newInfectionHistories an n*m matrix of 1s & 0s indicating infection histories, where n is individuals and m i strains
#' @param sampledIndivs the indices of sampled individuals to receive proposals
#' @param strainIsolationTimes the vector of strain isolation times in real time
#' @param ageMask the vector of indices for each individual specifiying which index of strainIsolationTimes is the first strain each individual coul dhave seen
#' @return a new matrix matching newInfectionHistories in dimensions with proposed moves
#' @export
infection_history_proposal <-function(newInfectionHistories,sampledIndivs,strainIsolationTimes,ageMask){
    newInf <- newInfectionHistories
    for(indiv in sampledIndivs){ # Resample subset of individuals
        rand1=runif(1)
        x=newInfectionHistories[indiv,ageMask[indiv]:length(strainIsolationTimes)] # Only resample years individual was alive
        ## Remove infection
        if(rand1<1/3){
            infectID= which(x>0)
            if(length(infectID)>0){
                x[sample(infectID,1)]=0 # Why double? DEBUG
            }
        }
        ## Add infection
        if(rand1>1/3 & rand1<2/3){
            ninfecID=which(x==0)
            if(length(ninfecID)>0){
                x[sample(ninfecID,1)]=1
            }
        }
        ## Move infection position
        if(rand1>2/3){
            infectID=which(x > 0)
            ninfecID=which(x == 0)
            if(length(infectID)>0 & length(ninfecID)>0){
                x[sample(infectID,1)]=0
                x[sample(ninfecID,1)]=1
            }
        }
        
        newInf[indiv,ageMask[indiv]:length(strainIsolationTimes)]=x # Only =1 if individual was alive
    } # end loop over individuals
    return(newInf)
}

#' Individual infection history sample - for testing
#'
#' @export
sample_indiv <- function(x){
    rand1=runif(1)
    
    ## Remove infection
    if(rand1<1/3){
        infectID= which(x>0)
        if(length(infectID)>0){
            x[sample(infectID,1)]=0 # Why double? DEBUG
        }
    }
    ## Add infection
    if(rand1>1/3 & rand1<2/3){
        ninfecID=which(x==0)
        if(length(ninfecID)>0){
            x[sample(ninfecID,1)]=1
        }
    }
    ## Move infection position
    if(rand1>2/3){
        infectID=which(x > 0)
        ninfecID=which(x == 0)
        if(length(infectID)>0 & length(ninfecID)>0){
            x[sample(infectID,1)]=0
            x[sample(ninfecID,1)]=1
        }
    }
    x
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

#' Propose initial infection histories
#'
#' Given a matrix of titre data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param dat the matrix of titres data with columns for individual, sample, and titre
#' @param strainIsolationTimes vector of real times for all strains
#' @param ageMask vector of indices for each individual corresponding to the first index of the strainIsolationTimes vector that each individual could be infected with
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of potential infecting strains
#' @export
setup_infection_histories<- function(dat, strainIsolationTimes, ageMask){
    # 
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
                # If high titre, set infection presence to 1 with some probability (0.2)
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
        forcedInfection <- which(strainIsolationTimes==sample(strainIsolationTimes[strainIsolationTimes <= sampleTime & strainIsolationTimes >= (sampleTime-10)],1))
        ## Pick strain within plausible region to add
        if(sum(infectionHistories[index, which(strainIsolationTimes < sampleTime)])==0) infectionHistories[index, forcedInfection] <- 1
        ## Make sure that no infections happened before birth
        if(ageMask[index] > 1) infectionHistories[index, 1:ageMask[index]] <- 0
        index <- index + 1
    }
    infectionHistories

}
