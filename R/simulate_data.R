#' Simulate group data
#'
#' Simulates a full set of titre data for n_indiv individuals with known theta and infectionHistories. Each individual gets nsamps random samples from sampleTimes, and infections can occur at any of strainIsolationTimes
#' @param n_indiv the number of individuals to simulate
#' @param theta the named parameter vector
#' @param infectionHistories the matrix of 1s and 0s giving presence/absence of infections for each individual
#' @param strainIsolationTimes the vector of strain circulation times (ie. possible infection times)
#' @param sampleTimes the vector of times at which samples could be taken
#' @param nsamps the number of samples per individual
#' @param antigenicMapLong the collapsed antigenic map for long term cross reactivity, after multiplying by sigma1
#' @param antigenicMapShort the collapsed antigenic map for short term cross reactivity, after multiplying by sigma2
#' @return a data frame with columns individual, samples, virus and titre of simulated data
#' @export
#' @seealso \code{\link{simulate_individual}}
simulate_group <- function(n_indiv, theta, infectionHistories,
                           strainIsolationTimes, sampleTimes,
                           nsamps, antigenicMapLong, antigenicMapShort){
    dat <- NULL
    ## For each individual
    for(i in 1:n_indiv){
        ## Choose random sampling times
        samps <- sample(sampleTimes, nsamps)
        samps <- samps[order(samps)]
        virusSamples <- rep(strainIsolationTimes, length(samps))
        dataIndices <- rep(length(strainIsolationTimes), length(samps))
        virusIndices <- match(virusSamples, strainIsolationTimes)-1
        y <- as.data.frame(simulate_individual(theta, infectionHistories[i,],
                                               samps, dataIndices, virusSamples,
                                               virusIndices,
                                               antigenicMapLong,
                                               antigenicMapShort,
                                               strainIsolationTimes))
        ## Record individual ID
        y$indiv <- i
        colnames(y) <- c("samples","virus","titre","individual")
        ## Combine data
        dat <- rbind(dat, y[,c("individual","samples","virus","titre")])
    }
    return(dat)
    
}

#' Simulate individual data
#'
#' Simulates a full set of titre data for an individual with known theta and infectionHistory. 
#' @param theta the named parameter vector
#' @param infectionHistory the vector of 1s and 0s giving presence/absence of infections
#' @param sampleTimes the vector of times at which samples are taken
#' @param strainIsolationTimes the vector of strain circulation times (ie. possible infection times)
#' @param antigenicMapLong the collapsed antigenic map for long term cross reactivity, after multiplying by sigma1
#' @param antigenicMapShort the collapsed antigenic map for short term cross reactivity, after multiplying by sigma2
#' @return a data frame with columns samples, virus and titre of simulated data
#' @export
#' @seealso \code{\link{infection_model_indiv}}
simulate_individual <- function(theta,
                                infectionHistory,
                                samplingTimes,
                                dataIndices,
                                strainIsolationTimes,
                                virusIndices,
                                antigenicMapLong,
                                antigenicMapShort,
                                strains){
    numberStrains <- length(strains)
    dat <- matrix(ncol=3, nrow=length(strainIsolationTimes))

    titres <- titre_data_individual(theta, infectionHistory, strains, seq_along(strains)-1, samplingTimes,
                                    dataIndices, match(strainIsolationTimes, strains)-1,
                                    antigenicMapLong, antigenicMapShort, numberStrains)
   
    dat[,1] <- rep(samplingTimes, dataIndices)
    dat[,2] <- strainIsolationTimes
    #dat[,3] <- add_noise(titres,theta)
    dat[,3] <- titres
    #dat[,3] <- rnorm(length(titres),titres,sd=theta["error"])
    return(dat)
}

#' Add noise
#'
#' Adds truncated noise to titre data
#' @param y the titre
#' @param theta a vector with MAX_TITRE and error parameters
#' @return a noisy titre
#' @export
add_noise <- function(y, theta){
    ## Draw from normal
    noise_y <- floor(rnorm(length(y), mean=y, sd=theta["error"]))

    ## If outside of bounds, truncate
    noise_y[noise_y < 0] <- 0
    noise_y[noise_y > theta["MAX_TITRE"]] <- theta["MAX_TITRE"]
    return(noise_y)    
}

#' Simulate attack rates
#'
#' Given a number of possible infection years, simulates attack rates from a log normal distribution with specified mean and standard deviation.
#' @param infectionYears the number of infection years
#' @param meanPar the mean of the log normal
#' @param sdPar the sd of the log normal
#' @param largeFirstYear simulate an extra large attach rate in the first year?
#' @param bigYearMean if large first year, what mean to use?
#' @return a vector of attack rates
#' @export
simulate_attack_rates <- function(infectionYears,meanPar=0.15,sdPar=0.5, 
                                  largeFirstYear=FALSE,bigYearMean=0.5){
    attack_year <- rlnorm(infectionYears, meanlog=log(meanPar)-sdPar^2/2,sdlog=sdPar)
    if(largeFirstYear) attack_year[1] <- rlnorm(1,meanlog=log(bigYearMean)-(sdPar/2)^2/2,sdlog=sdPar/2)
    return(attack_year)
}

#' Simulate infection histories
#'
#' Given a vector of infection probabilities and potential infection times, simulates infections for each element of ages (ie. each element is an individual age. Only adds infections for alive individuals)
#' @param pInf a vector of attack rates (infection probabilities) for each year
#' @param infSD the standard deviation of attack rate draws
#' @param strainIsolationTimes the vector of possible infection times
#' @param samplingTimes vector of potential sampling times
#' @param ages a vector of ages for each individual
#' @return a matrix of infection histories for each individual in ages
#' @export
simulate_infection_histories <- function(pInf, infSD, strainIsolationTimes, samplingTimes, ages){
    n_strains <- length(pInf) # How many strains
    n_indiv <- length(ages) # How many individuals
    indivs <- 1:n_indiv
    ## Empty matrix
    infectionHistories <- matrix(0,ncol=n_strains,nrow=n_indiv)
    ## Simulate attack rates
    attackRates <- rlnorm(length(pInf), meanlog=log(pInf)- infSD^2/2,sdlog=infSD)

    ## Should this be necessary?
    attackRates[attackRates > 1] <- 1

    ## For each strain (ie. each infection year)
    for(i in 1:n_strains){
        ## Find who was alive (all we need samplingTimes for is its max value)
        alive <- (max(samplingTimes) - ages) <= strainIsolationTimes[i]

        ## Sample a number of infections for the alive individuals, and set these entries to 1
        infectionHistories[sample(indivs[alive], round(length(indivs[alive])*attackRates[i])),i] <- 1
    }
    return(infectionHistories)    
}

#' Simulate full data set
#'
#' Simulates a full data set for a given set of parameters etc.
#' @param parTab the full parameter table controlling parameter ranges and values
#' @param group which group index to give this simulated data
#' @param n_indiv number of individuals to simulate
#' @param strainIsolationTimes vector of strain ciruclation times
#' @param samplingTimes possible sampling times for the individuals
#' @param nsamps the number of samples each individual has
#' @param antigenicMap the raw antigenic map with colnames x_coord, y_coord and inf_years
#' @param sampleSensoring DEPRECATED - what proportion of samples are randomly missing?
#' @param titreSensoring what proportion of titres are randomly missing?
#' @param ageMin minimum age to simulate
#' @param ageMax maximum age to simulate
#' @param simInfPars vector of parameters to pass to \code{\link{simulate_attack_rates}}
#' @return a list with: 1) the data frame of titre data as returned by \code{\link{simulate_group}}; 2) a matrix of infection histories as returned by \code{\link{simulate_infection_histories}}; 3) a vector of ages
#' @export
simulate_data <- function(parTab, group=1,n_indiv,
                          strainIsolationTimes, samplingTimes, nsamps=2,
                          antigenicMap,
                          sampleSensoring=0, titreSensoring=0,
                          ageMin=5,ageMax=80,
                          simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1)){
    ## Extract parameters
    pars <- parTab$values
    names(pars) <- parTab$names

    ## Create antigenic map for short and long term boosting
    antigenicMap1 <- outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")])
    
    antigenicMapLong <- 1 - pars["sigma1"]*c(antigenicMap1)
    antigenicMapShort <- 1 - pars["sigma2"]*c(antigenicMap1)

    antigenicMapLong[antigenicMapLong < 0] <- 0
    antigenicMapShort[antigenicMapShort < 0] <- 0

    ## Simulate ages
    ages <- floor(runif(n_indiv,ageMin,ageMax))

    ## Simulate attack rates
    pInf <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])

    ## Simulate infection histories
    infHist <- simulate_infection_histories(pInf, simInfPars["logSD"], strainIsolationTimes, samplingTimes, ages)

    ## Simulate titre data
    y <- simulate_group(n_indiv, pars, infHist, strainIsolationTimes, samplingTimes,
                        nsamps, antigenicMapLong,antigenicMapShort)

    ## Randomly censor titre values
    y$titre <- y$titre*sample(c(0,1),nrow(y),prob=c(titreSensoring,1-titreSensoring),replace=TRUE)
    
    y <- cbind(group=group,y)
    DOB <- max(samplingTimes) - ages
    y <- y[,c("individual","virus","samples","titre","group")]
    return(list(data=y, infectionHistories=infHist, ages=DOB))
}

#' Create useable antigenic map
#'
#' Creates an antigenic map from an input data frame that can be used to calculate cross reactivity. This will end up being an NxN matrix, where there are N strains circulating.
#' @param anti.map.in can either be a 1D antigenic line to calculate distance from, or a two dimensional matrix with x and y coordinates on an antigenic map
#' @return the euclidean antigenic distance between each pair of viruses in anti.map.in
#' @export
outputdmatrix.fromcoord <- function(anti.map.in){ #anti.map.in can be vector or matrix - rows give inf_years, columns give location
                                        # Calculate antigenic distances
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
                                        # If 1D antigenic 'line' defined, calculate distances directory from input
        (dmatrix=sapply(anti.map.in,function(x){y=abs(anti.map.in-x); y   })) 
    } else { # If 2D antigenic map defined, calculate distances directory from input
        (dmatrix=apply(anti.map.in,1,
                       function(x){
                           y=sqrt(colSums(apply(anti.map.in,1,function(y){(y-x)^2})));
                           y 
                       }))
    }
}
