# - - - - - - - - - - - - - - - - 
# Acceptance probability in MH algorithm
ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  # Flat priors on theta => symmetric update probability
  calc.lik = exp(marg_likelihood_star-marg_likelihood)
  calc.lik[calc.lik>1]=1 
  calc.lik
}

#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' @param covMat the 2D covariance matrix for all of the parameters
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
#' @useDynLib serology-model
mvr_proposal <- function(values, fixed, covMat){
    proposed <- values
    proposed[fixed] <- MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat)
    return(proposed)
}


## Age mask is a value indicating from which strain index an individual may be infected
infection_history_proposal <-function(newInfectionHistories,sampledIndivs,strainIsolationTimes,ageMask){
    for(i in sampledIndivs){ # Resample subset of individuals
        rand1=runif(1)
        x=newInfectionHistories[i,ageMask[i]:length(strainIsolationTimes)] # Only resample years individual was alive
        infvector=c(1:length(x))
        infvector2=rev(infvector)
        
        ## Remove infection
        if(rand1<1/3){
            infectID=infvector[(as.numeric(x)>0)]
            if(length(infectID)>0){
                x[sample(c(infectID),1)]=0 # Why double? DEBUG
            }
        }
        ## Add infection
        if(rand1>1/3 & rand1<2/3){
            ninfecID=infvector[(as.numeric(x)==0)]
            if(length(ninfecID)>0){
                x[sample(c(ninfecID),1)]=1
            }
        }
        ## Move infection position
        if(rand1>2/3){
            infectID=infvector[(as.numeric(x)>0)]
            ninfecID=infvector[(as.numeric(x)==0)]
            
            if(length(infectID)>0 & length(ninfecID)>0){
                x[sample(c(infectID),1)]=0
                x[sample(c(ninfecID),1)]=1
            }
        }
        newInfectionHistories[i,ageMask[i]:length(strainIsolationTimes)]=x # Only =1 if individual was alive
    } # end loop over individuals
    newInfectionHistories
}
