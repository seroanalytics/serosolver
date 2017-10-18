#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' @param covMat the 2D covariance matrix for all of the parameters
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
#' @useDynLib sero-solver
mvr_proposal <- function(values, fixed, covMat, useLog=FALSE){
    proposed <- values
    if(!useLog) proposed[fixed] <- MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat)
    else proposed[fixed] <- exp(MASS::mvrnorm(n=1,mu=log(proposed[fixed]),Sigma=(5.6644/length(fixed))*covMat))
    return(proposed)
}


#' Age mask is a value indicating from which strain index an individual may be infected
#' @export
infection_history_proposal <-function(newInfectionHistories,sampledIndivs,strainIsolationTimes,ageMask){
    for(i in sampledIndivs){ # Resample subset of individuals
        rand1=runif(1)
        x=newInfectionHistories[i,ageMask[i]:length(strainIsolationTimes)] # Only resample years individual was alive
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
        
        newInfectionHistories[i,ageMask[i]:length(strainIsolationTimes)]=x # Only =1 if individual was alive
    } # end loop over individuals
    newInfectionHistories
}

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
#' @useDynLib sero-solver
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
#' @useDynLib sero-solver
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}


#' Protect function
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
#' @param f the function to be protected
#' @return the protected function
#' @export
#' @useDynLib sero-solver
protect <- function(f){
    function(...){
        tryCatch(f(...),error=function(e){
            message("caught error: ", e$message)
            -100000
        })
    }
}

#' @export
ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  # Flat priors on theta => symmetric update probability
  calc.lik = exp(marg_likelihood_star-marg_likelihood)
  calc.lik[calc.lik>1]=1 
  calc.lik
}


toUnitScale <- function(x, min, max){
    return((x-min)/(max-min))
}
fromUnitScale <- function(x,min,max){
    return(min + (max-min)*x)
}
