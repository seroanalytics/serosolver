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
    move <- sample((-moveSize):moveSize,1)
    y2 <- y1 + move
    if(y2 < 0) y2 = y2 + ncol(infHist)
    if(y2 >= ncol(infHist)) y2 = y2 - floor(y2/ncol(infHist))*ncol(infHist)
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
