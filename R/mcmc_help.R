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

#' @export
infection_history_proposal_smart <- function(infectionHistories, sampledIndivs, ageMask, noJumps){
    newInf <- infectionHistories
    maxI <- ncol(newInf)
    for(i in seq_along(sampledIndivs)){
        indiv <- sampledIndivs[i]
        x = newInf[indiv, ageMask[indiv]:maxI]
        for(j in 1:noJumps[i]){
            rand1 = runif(1)
            if(rand1 < 1/2){
                id1 <- sample(length(x),1)
                x[id1] <- !x[id1]
            } else {
                id1 <- sample(length(x),1)
                moveMax <- 5
                ##       
                move <- sample(-moveMax:moveMax,1)
                id2 <- id1 + move
                                        #
                if(id2 < 1) id2 <- maxI + id2
                if(id2 > maxI) id2 <- id2 - maxI
                tmp <- x[id1]
                x[id1] <- x[id2]
                x[id2] <- tmp       
            }
           
        }
        newInf[indiv,ageMask[indiv]:maxI]=x # Only =1 if individual was alive
    } # end loop over individuals
    return(newInf)
}


#' @export
infection_history_betabinom_symmetric<- function(newInfHist, sampledIndivs, ageMask, moveSizes, alpha, beta){
    newInf <- newInfHist
    for(indiv in sampledIndivs){
        x <- newInfHist[indiv, ageMask[indiv]:ncol(newInfHist)]
        maxI <- length(x)
        rand1 <- runif(1)
        if(rand1 < 1/2){
            loc <- sample(length(x), 1)            
            x[loc] <- !x[loc]                                       
        } else {
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
        newInf[indiv,ageMask[indiv]:ncol(newInfHist)]=x
    }
    return(newInf)    
}

#' @export
infection_history_betabinom<- function(newInfHist, sampledIndivs, ageMask, moveSizes, alpha, beta){
    newInf <- newInfHist
    for(indiv in sampledIndivs){
#        message(cat("Indiv: ", indiv,sep="\t"))
        x <- newInfHist[indiv, ageMask[indiv]:ncol(newInfHist)]
        maxI <- length(x)
        rand1 <- runif(1)
        if(rand1 < 1){
            loc <- sample(length(x), 1)            
            x_new <- x_old <- x
            x_new[loc] <- 1
            x_old[loc] <- 0
            ##probA <- dbb(sum(x_new), length(x), alpha, beta)/choose(length(x), sum(x_new))
            ##probB <- dbb(sum(x_old), length(x), alpha, beta)/choose(length(x), sum(x_old))
            ##ratio <- probA/(probA + probB)
            prob1 <- (alpha+sum(x[-loc]))/(alpha+beta+(length(x)-1))
            if(runif(1)<prob1){
                x[loc] <- 1
            } else {
                x[loc] <- 0
            }
        } else {
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
        if(length(x) != length(newInf[indiv, ageMask[indiv]:ncol(newInfHist)])){
            message(cat("Indiv: ", indiv, sep="\t"))
            message(cat("Length x: ", length(x), sep="\t"))
            message(cat("Length y: ", length(newInf[indiv, ageMask[indiv]:ncol(newInfHist)]), sep="\t"))
            
        }

        newInf[indiv,ageMask[indiv]:ncol(newInfHist)]=x
    }
    return(newInf)    
}



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


#' Infection history proposal
#'
#' Proposes new infection histories for a vector of infection histories, where rows represent individuals and columns represent years. Proposals are either removal, addition or switching of infections.
#' Also requires the indices of sampled individuals, the vector of strain isolation times, and a vector of age masks (ie. which index of the strainIsolationTimes vector is the first year in which
#' an individual *could* be infected).
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param newInfectionHistories an n*m matrix of 1s & 0s indicating infection histories, where n is individuals and m i strains
#' @param sampledIndivs the indices of sampled individuals to receive proposals
#' @param strainIsolationTimes the vector of strain isolation times in real time
#' @param ageMask the vector of indices for each individual specifiying which index of strainIsolationTimes is the first strain each individual could have seen
#' @return a new matrix matching newInfectionHistories in dimensions with proposed moves
#' @export
infection_history_proposal_OLD<-function(newInfectionHistories,sampledIndivs,strainIsolationTimes,ageMask){
    newInf <- newInfectionHistories
    acceptance_distribution <- rep(1, nrow(newInf))
    for(indiv in sampledIndivs){ # Resample subset of individuals
        rand1=runif(1)
        x=newInfectionHistories[indiv,ageMask[indiv]:length(strainIsolationTimes)] # Only resample years individual was alive
        accept <- 0
      
        ## Remove infection
        if(rand1<1/3){
            infectID= which(x>0)
                                        # Number of 1s in first place
            n_1 <- length(infectID)
            if(n_1 > 0){
                x[sample(infectID,1)]=0 # Why double? DEBUG
                n_0 <- length(which(x==0))
                p_forward <- 1/3 * 1/n_1
                p_back <- 1/3 * 1/n_0
                accept <- (p_back/p_forward)
            }
        }
        ## Add infection
        if(rand1>1/3 & rand1<2/3){
            ninfecID=which(x==0)
            n_0 <- length(ninfecID)
            if(n_0>0){
                x[sample(ninfecID,1)]=1
                n_1 <- length(which(x==1))
                p_forward <- 1/3 * 1/n_0
                p_back <- 1/3 * 1/n_1
                accept <- (p_back/p_forward)
            }
        }
        ## Move infection position
        if(rand1>2/3){
            infectID=which(x > 0)
            ninfecID=which(x == 0)
            n_1 <- length(infectID)
            n_0 <- length(ninfecID)
            if(n_1 > 0 & n_0 > 0){
                x[sample(infectID,1)]=0
                x[sample(ninfecID,1)]=1
                p_forward <- 1/3 * n_1 * n_0
                
                n_1b <- length(which(x > 0))
                n_0b <- length(which(x==0))
                p_back <- 1/3 * n_1b * n_0b
                accept <- (p_back/p_forward)
            }
        }
        acceptance_distribution[indiv] <- accept
        newInf[indiv,ageMask[indiv]:length(strainIsolationTimes)]=x # Only =1 if individual was alive
    } # end loop over individuals
return(list(newInf, acceptance_distribution))
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
setup_infection_histories<- function(dat, strainIsolationTimes, ageMask,sample_prob, titre_cutoff=3){
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

#' @export
setup_infection_histories_new <- function(data, sages, trainIsolationTimes, space=5, titre_cutoff=2){
  startInf <- NULL
  individuals <- unique(data$individual)
  for(individual in individuals){
    dat <- data[data$individual == individual,]
    dob <- as.numeric(ages[ages$individual == individual,"DOB"])
    strains <- unique(dat$virus)
    strains <- strains[strains >= dob]
    infYears <- NULL
    i <- 0
    while(i < length(strains)){ ## Start going through each strain
      i <- i + 1
      strain <- strains[i] ## Get current strain of interest
      #message(cat("Strain under consideration: ", strain, sep="\t"))
      measurement <- -1
      dist <- 0
      
      titre <- max(dat[dat$virus == strain, "titre"]) ## Get max titre against this strain
      if(titre >= titre_cutoff){ ## If elevated against this strain, assume an infection
       # message(cat("Elevated titre against: ", strain,sep="\t"))
        newInf <- strain
        ## Begin counting up distance
        while(dist < space & i < length(strains)){
          i <- i + 1
          strain <- strains[i]
          #message(cat("Strain under consideration: ", strain, sep="\t"))
          dist <- strain - newInf ## How many years since last infection?
          new_titre <- max(dat[dat$virus == strain, "titre"]) ## Get max titre against next strain along
          ## If this one is better, replace and reset distance
          if(new_titre > titre){
            #message(cat("Better titre against: ", strain, "; at: ",new_titre,sep="\t"))
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
