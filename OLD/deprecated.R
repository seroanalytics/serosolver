
#' DEPRECATED
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
