#' Swap infection history years
#'
#' Swaps the entire contents of two columns of the infection history matrix, adhering to age and strain limitations.
#' @param infection_history matrix of 1s and 0s to swap around
#' @param age_mask the first index in infection_history that each individual (row) could be infected in
#' @param strain_mask the last index in infection_history that each individual (row) could be infected in ie. the time of the latest blood sample
#' @param swap_propn what proportion of infections should be swapped?
#' @param move_size How many time points away should be chosen as candidate swaps?
#' @return the same infection_history matrix, but with two columns swapped
#' @export
inf_hist_swap <- function(infection_history, age_mask, strain_mask, swap_propn, move_size){
    ## Choose a column
    y1 <- sample(1:ncol(infection_history), 1)

    ## Propose another column some random distance, but not 0, away
    move <- 0
    while(move == 0) move <- sample((-move_size):move_size,1)

    ## Need to adjust if we've proposed too far away
    y2 <- y1 + move
    while(y2 < 1) y2 = y2 + ncol(infection_history)
    if(y2 > ncol(infection_history)) y2 = y2 - floor(y2/ncol(infection_history))*ncol(infection_history)

    ## Get the first and last year chronologically
    small_year <- min(y1, y2)
    big_year <- max(y1, y2)

    ## Find individuals that are alive/sampled in both years and choose the lesser of swap_propn*n_indivs and
    ## the number that are actually able to be infected in both years
    indivs <- 1:nrow(infection_history)
    alive_indivs <- indivs[intersect(which(age_mask <= small_year), which(strain_mask >= big_year))]
    samp_indivs <- sample(alive_indivs, floor(length(alive_indivs)*swap_propn))

    ## Swap contents
    tmp <- infection_history[samp_indivs, y1]
    infection_history[samp_indivs, y1] <- infection_history[samp_indivs, y2]
    infection_history[samp_indivs, y2] <- tmp
    return(infection_history)
    
}
#' Swap infection history years with lambda term
#'
#' Swaps the entire contents of two columns of the infection history matrix, adhering to age and strain limitations. Also swaps the values of lambda that correspond to these years
#' @param infection_history matrix of 1s and 0s to swap around
#' @param lambdas vector of force of infection parameters for each column
#' @param age_mask the first index in infection_history that each individual (row) could be infected in
#' @param strain_mask the last index in infection_history that each individual (row) could be infected in ie. the time of the latest blood sample
#' @param swap_propn what proportion of infections should be swapped?
#' @param move_size How many time points away should be chosen as candidate swaps?
#' @param n_alive number of individuals alive in each entry of lambdas
#' @return a list: the same infection_history matrix, but with two columns swapped; also the swapped lambdas
#' @seealso \code{\link{inf_hist_swap}}
#' @export
inf_hist_swap_lambda <- function(infection_history, lambdas, age_mask, strain_mask, swap_propn, move_size,n_alive){
    ## This first bit of code is the same as inf_hist_swap
    y1 <- sample(1:ncol(infection_history), 1)
    move <- 0
    while(move == 0) move <- sample((-move_size):move_size,1)
    y2 <- y1 + move
    
    while(y2 < 1) y2 = y2 + ncol(infection_history)
    if(y2 > ncol(infection_history)) y2 = y2 - floor(y2/ncol(infection_history))*ncol(infection_history)
    small_year <- min(y1, y2)
    big_year <- max(y1, y2)
    
    indivs <- 1:nrow(infection_history)
    
    alive_indivs <- indivs[intersect(which(age_mask <= small_year), which(strain_mask >= big_year))]
    samp_indivs <- sample(alive_indivs, floor(length(alive_indivs)*swap_propn))
    
    tmp <- infection_history[samp_indivs, y1]

    ## Get number of infections amongst selected individuals in the two years
    no_infs_y1 <- sum(tmp)
    no_infs_y2 <- sum(infection_history[samp_indivs, y2])

    ## Total number of infections in these years
    total_infs_y1 <- sum(infection_history[,y1])
    total_infs_y2 <- sum(infection_history[,y2])

    ## How many infections will there be in these two years after the swap?
    new_total_infs_y1 <- total_infs_y1 - no_infs_y1 + no_infs_y2
    new_total_infs_y2 <- total_infs_y2 - no_infs_y2 + no_infs_y1

    ## Adjust the FOI parameters proportional to the number of infections gained/lost
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

    ## And finish the swap
    infection_history[samp_indivs, y1] <- infection_history[samp_indivs, y2]
    infection_history[samp_indivs, y2] <- tmp
    
    return(list(infection_history, lambdas))
}


#' R implementation of the infection history gibbs proposal, now in C++ code
#' @export
infection_history_proposal_gibbs_R <- function(pars, infection_history, indivSampPropn,
                                               n_years_samp, age_mask,
                                               n_alive,
                                               alpha, beta,
                                               strains, strainIndices, sampleTimes,
                                               indicesData, indicesDataOverall,
                                               indicesSamples, virusIndices, antigenic_map_long, antigenic_map_short,
                                               titres){
    n_indivs <- nrow(infection_history)
    n_strains <- ncol(infection_history)

    newInfHist <- infection_history

    for(i in 1:n_indivs){
        if(runif(1) < indivSampPropn){
            
        #message(cat("Indiv: ", i))
            sample_years <- seq(age_mask[i],n_strains,by=1)
            #message(cat("Sample years: " ,sample_years,sep="\t"))
            n_samp_max <- length(sample_years)
            n_samp_max <- min(n_years_samp, n_samp_max)
            samps <- seq(1,n_samp_max,by=1)
            
            locs <- sample(samps, n_samp_max,replace=FALSE)
            #locs <- sample(age_mask[i]:ncol(infection_history),n_samp_max)
            #message(cat("Locs: ", locs,sep="\t"))
                                        ##locs <- sample(1:ncol(infection_history),n_years_samp)
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
