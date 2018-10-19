#' @export
infection_history_proposal_gibbs_R <- function(pars, infHist, indivSampPropn,
                                               n_years_samp, ageMask,
                                               n_alive,
                                               alpha, beta,
                                               strains, strainIndices, sampleTimes,
                                               indicesData, indicesDataOverall,
                                               indicesSamples, virusIndices, antigenicMapLong, antigenicMapShort,
                                               titres){
    n_indivs <- nrow(infHist)
    n_strains <- ncol(infHist)

    newInfHist <- infHist

    for(i in 1:n_indivs){
        if(runif(1) < indivSampPropn){
            
        #message(cat("Indiv: ", i))
            sample_years <- seq(ageMask[i],n_strains,by=1)
            #message(cat("Sample years: " ,sample_years,sep="\t"))
            n_samp_max <- length(sample_years)
            n_samp_max <- min(n_years_samp, n_samp_max)
            samps <- seq(1,n_samp_max,by=1)
            
            locs <- sample(samps, n_samp_max,replace=FALSE)
            #locs <- sample(ageMask[i]:ncol(infHist),n_samp_max)
            #message(cat("Locs: ", locs,sep="\t"))
                                        ##locs <- sample(1:ncol(infHist),n_years_samp)
            for(j in 1:length(locs)){                
             #   indivHist <- newInfHist[i,]
                
                year <- sample_years[locs[j]]
                #year <- locs[j]
               # message(cat("Loc: ", locs[j]))
               # message(cat("Year: ", year))
                m <- sum(newInfHist[,year]) - newInfHist[i,year]
                #m <- sum(newInfHist[-i,year])
               # message(cat("Infections: ",m))
                n <- n_alive[year] - 1
                #n <- n_indivs - 1
               # message(cat("Alive: ",n))
                ratio <- (m+alpha)/(n+alpha+beta)

                rand1 <- runif(1)
                if(rand1 < ratio){
                    new_entry <- 1
                } else {
                    new_entry <- 0
                }

                if(new_entry!= newInfHist[i,year]){
                    old_prob <- new_prob <- 1
                    log_prob <- min(0, new_prob - old_prob)
                    rand1 <- runif(1)
                    if(log(rand1) < log_prob){
                        newInfHist[i,year] <- new_entry
                    }
                }
            }
        }
    }
    return(newInfHist)
}

inf_hist_swap <- function(infHist, ageMask, strainMask, swapPropn, moveSize){
    y1 <- sample(1:ncol(infHist), 1)
    move <- 0
    while(move == 0) move <- sample((-moveSize):moveSize,1)
    y2 <- y1 + move
    if(y2 < 1) y2 = y2 + ncol(infHist)
    if(y2 > ncol(infHist)) y2 = y2 - floor(y2/ncol(infHist))*ncol(infHist)
    smallYear <- min(y1, y2)
    bigYear <- max(y1, y2)
    indivs <- 1:nrow(infHist)
    alive_indivs <- indivs[intersect(which(ageMask <= smallYear), which(strainMask >= bigYear))]
    samp_indivs <- sample(alive_indivs, floor(length(alive_indivs)*swapPropn))
    tmp <- infHist[samp_indivs, y1]
    infHist[samp_indivs, y1] <- infHist[samp_indivs, y2]
    infHist[samp_indivs, y2] <- tmp
    return(infHist)
    
}

#' @export
inf_hist_swap_lambda <- function(infHist, lambdas, ageMask, strainMask, swapPropn, moveSize,n_alive){
    y1 <- sample(1:ncol(infHist), 1)
    move <- 0
    while(move == 0) move <- sample((-moveSize):moveSize,1)
    y2 <- y1 + move
    if(y2 < 1) y2 = y2 + ncol(infHist)
    if(y2 > ncol(infHist)) y2 = y2 - floor(y2/ncol(infHist))*ncol(infHist)
    smallYear <- min(y1, y2)
    bigYear <- max(y1, y2)
    indivs <- 1:nrow(infHist)
    alive_indivs <- indivs[intersect(which(ageMask <= smallYear), which(strainMask >= bigYear))]
    samp_indivs <- sample(alive_indivs, floor(length(alive_indivs)*swapPropn))
    
    tmp <- infHist[samp_indivs, y1]

    no_infs_y1 <- sum(tmp)
    no_infs_y2 <- sum(infHist[samp_indivs, y2])

    total_infs_y1 <- sum(infHist[,y1])
    total_infs_y2 <- sum(infHist[,y2])

    new_total_infs_y1 <- total_infs_y1 - no_infs_y1 + no_infs_y2
    new_total_infs_y2 <- total_infs_y2 - no_infs_y2 + no_infs_y1

    lost_foi_y1 <- gained_foi_y1 <- 0
    if(total_infs_y1 > 0){
        lost_foi_y1 <- no_infs_y1/total_infs_y1
    }
    if(new_total_infs_y1 > 0){
        gained_foi_y1 <- no_infs_y2/new_total_infs_y1
    }
    lambdas[y1] <- lambdas[y1] - lambdas[y1]*lost_foi_y1 + lambdas[y1]*gained_foi_y1
    
    lost_foi_y2 <- gained_foi_y2 <- 0
    if(total_infs_y2 > 0){
        lost_foi_y2 <- no_infs_y2/total_infs_y2
    }
    if(new_total_infs_y2 > 0){
        gained_foi_y2 <- no_infs_y1/new_total_infs_y2
    }
    lambdas[y2] <- lambdas[y2] - lambdas[y2]*lost_foi_y2 + lambdas[y2]*gained_foi_y2
#    calc_change_lambda <- function(delta_y, n, lambda, infs){
#        new_lambda <- lambda
#        if(delta_y < 0){
#            new_lambda <- lambda + (lambda/infs)*delta_y
#        }
#        if(delta_y > 0){
#            new_lambda <- lambda + ((1-lambda)/(n-infs))*delta_y
#        }
#        return(new_lambda)           
#    }
    
    
    infHist[samp_indivs, y1] <- infHist[samp_indivs, y2]
    infHist[samp_indivs, y2] <- tmp
    
    return(list(infHist, lambdas))
}
