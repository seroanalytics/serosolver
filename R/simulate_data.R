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
#' @param repeats number of repeated samples for each year
#' @param mus default NULL, optional vector of boosting parameters for each strain
#' @param mu_indices default NULL, optional vector giving the index of `mus` that each strain uses the boosting parameter from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @param measurement_bias default NULL, optional vector of measurement bias shift parameters for each strain
#' @param measurement_indices default NULL, optional vector giving the index of `measurement_bias` that each strain uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @return a data frame with columns individual, samples, virus and titre of simulated data
#' @export
#' @seealso \code{\link{simulate_individual}}
simulate_group <- function(n_indiv, theta, infectionHistories,
                           strainIsolationTimes, sampleTimes,
                           nsamps, antigenicMapLong, antigenicMapShort,
                           repeats=1, mus=NULL, mu_indices=NULL,
                           measurement_bias=NULL, measurement_indices=NULL,
                           addNoise=TRUE){
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
                                               strainIsolationTimes,
                                               repeats, mus, mu_indices,
                                               measurement_bias, measurement_indices,
                                               addNoise))
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
#' @param dataIndices see \code{\link{create_post_func}}
#' @param strainIsolationTimes see \code{\link{create_post_func}}
#' @param virusIndices see \code{\link{create_post_func}}
#' @param strainIsolationTimes the vector of strain circulation times (ie. possible infection times)
#' @param antigenicMapLong the collapsed antigenic map for long term cross reactivity, after multiplying by sigma1
#' @param antigenicMapShort the collapsed antigenic map for short term cross reactivity, after multiplying by sigma2
#' @param strains vector of all possible circulating strains
#' @param repeats number of repeated samples for each year
#' @param mus default NULL, optional vector of boosting parameters for each strain
#' @param mu_indices default NULL, optional vector giving the index of `mus` that each strain uses the boosting parameter from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @param measurement_bias default NULL, optional vector of measurement bias shift parameters for each strain
#' @param measurement_indices default NULL, optional vector giving the index of `measurement_bias` that each strain uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clus
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
                                strains,
                                repeats=1,
                                mus=NULL,
                                mu_indices=NULL,
                                measurement_bias=NULL,
                                measurement_indices=NULL,
                                addNoise=TRUE){
    numberStrains <- length(strains)
    dat <- matrix(ncol=3, nrow=length(strainIsolationTimes)*repeats)
    
    if(is.null(mu_indices)){    
        titres <- titre_data_individual(theta, infectionHistory, strains, seq_along(strains)-1, samplingTimes,
                                        dataIndices, match(strainIsolationTimes, strains)-1,
                                        antigenicMapLong, antigenicMapShort, numberStrains)
    } else {
        titres <- titre_data_individual_mus(theta, mus, infectionHistory, strains, seq_along(strains)-1,
                                        mu_indices, samplingTimes,
                                        dataIndices, match(strainIsolationTimes, strains)-1,
                                        antigenicMapLong, antigenicMapShort, numberStrains)
    }
    titres <- rep(titres, each=repeats)
    samplingTimes <- rep(samplingTimes, dataIndices)
    samplingTimes <- rep(samplingTimes, each=repeats)
    dat[,1] <- samplingTimes
    dat[,2] <- rep(strainIsolationTimes,each=repeats)
    if(addNoise){
        if(!is.null(measurement_indices)){
            dat[,3] <- add_noise(titres,theta, measurement_bias, measurement_indices[match(dat[,2], strains)])
        } else {
            dat[,3] <- add_noise(titres,theta, NULL, NULL)
        }
    } else {
        dat[,3] <- titres
    }
    #dat[,3] <- titres
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
add_noise <- function(y, theta, measurement_bias=NULL, indices=NULL){
    ## Draw from normal
    if(!is.null(measurement_bias)) {
        noise_y <- floor(rnorm(length(y), mean=y+measurement_bias[indices], sd=theta["error"]))
    } else {
        noise_y <- floor(rnorm(length(y), mean=y, sd=theta["error"]))
    }

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
#' @return a list with a matrix of infection histories for each individual in ages and the true attack rate for each epoch
#' @export
simulate_infection_histories <- function(pInf, infSD, strainIsolationTimes, samplingTimes, ages){
    n_strains <- length(pInf) # How many strains
    n_indiv <- length(ages) # How many individuals
    indivs <- 1:n_indiv
    ## Empty matrix
    infectionHistories <- matrix(0,ncol=n_strains,nrow=n_indiv)
    
    ## Simulate attack rates
    attackRates <- pInf
    #attackRates <- rlnorm(length(pInf), meanlog=log(pInf)- infSD^2/2,sdlog=infSD)

    ## Should this be necessary?
    #attackRates[attackRates > 1] <- 1
    ARs <- numeric(n_strains)
    ## For each strain (ie. each infection year)
    for(i in 1:n_strains){
        ## Find who was alive (all we need samplingTimes for is its max value)
        alive <- (max(samplingTimes) - ages) <= strainIsolationTimes[i]
        
        ## Sample a number of infections for the alive individuals, and set these entries to 1
        y <- round(length(indivs[alive])*attackRates[i])
        #y <- rbinom(1, length(indivs[alive]),attackRates[i])
        ARs[i] <- y/length(indivs[alive])
        x <- sample(indivs[alive], y)
        infectionHistories[x,i] <- 1
    }
    return(list(infectionHistories,ARs))    
}

#' Generates attack rates from an SIR model with fixed beta/gamma, specified final attack rate and the number of time "buckets" to solve over ie. buckets=12 returns attack rates for 12 time periods
#' @export
generate_ar_annual <- function(AR, buckets){
  SIR_odes <- function(t, x, params) {
    S <- x[1]
    I <- x[2]
    R <- x[3]
    inc <- x[4]
    
    beta <- params[1]
    gamma <- params[2]
    dS <- -beta*S*I
    dI <- beta*S*I - gamma*I
    dR <- gamma*I
    dinc <- beta*S*I
    list(c(dS,dI,dR, dinc))
  }
  R0 <- 1.2
  gamma <- 1/5
  beta <- R0*gamma
  t <- seq(0,360,by=0.1)
  results <- as.data.frame(deSolve::ode(y=c(S=1,I=0.0001,R=0, inc=0),
                                        times=t, func=SIR_odes,
                                        parms=c(beta,gamma)))
  incidence <- diff(results$inc)
  incidence <- incidence*AR/sum(incidence)
  group <- 360*10/buckets
  monthly_risk <- colSums(matrix(incidence, nrow=group))
  return(monthly_risk)
}

#' @export
simulate_ars_buckets <- function(infectionYears, buckets, meanPar=0.15,sdPar=0.5, 
                                 largeFirstYear=FALSE,bigYearMean=0.5){
    n <- ceiling(length(infectionYears)/buckets)
    attack_year <- rlnorm(n, meanlog=log(meanPar)-sdPar^2/2,sdlog=sdPar)
    if(largeFirstYear) attack_year[1] <- rlnorm(1,meanlog=log(bigYearMean)-(sdPar/2)^2/2,sdlog=sdPar/2)
    ars <- NULL
    
    for(i in seq_along(attack_year)){
        ars <- c(ars, generate_ar_annual(attack_year[i],buckets))
    }
    
    ars <- ars[1:length(infectionYears)]
    return(ars)
}
#' @export
simulate_ars_spline <- function(infectionYears, buckets, meanPar=0.15,sdPar=0.5, largeFirstYear=FALSE,bigYearMean=0.5, knots,theta){
    infectionYears <- infectionYears[seq(1,length(infectionYears),by=buckets)]/buckets
    n <- length(infectionYears)
    attack_year <- rlnorm(n, meanlog=log(meanPar)-sdPar^2/2,sdlog=sdPar)
    if(largeFirstYear) attack_year[1] <- rlnorm(1,meanlog=log(bigYearMean)-(sdPar/2)^2/2,sdlog=sdPar/2)
    ars <- generate_lambdas(attack_year, knots, theta, n, buckets) 
    return(ars)
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
#' @param useSIR boolean specifying whether to sample from an SIR model or just log normal distribution for each epoch. If TRUE, uses an SIR model
#' @param useSpline if TRUE, calculates a spline to generate sub annual FOI from
#' @param pInf if given, provides infection probabilities to sample from rather than simulating
#' @param repeats number of repeats for each year
#' @param mu_indices default NULL, optional vector giving the index of `mus` that each strain uses the boosting parameter from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @return a list with: 1) the data frame of titre data as returned by \code{\link{simulate_group}}; 2) a matrix of infection histories as returned by \code{\link{simulate_infection_histories}}; 3) a vector of ages
#' @export
simulate_data <- function(parTab, group=1,n_indiv,buckets=12,
                          strainIsolationTimes, samplingTimes, nsamps=2,
                          antigenicMap,
                          sampleSensoring=0, titreSensoring=0,
                          ageMin=5,ageMax=80,
                          simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),
                          theta=rep(0.1,5),
                          knots=c(0.33,0.66),
                          useSIR=FALSE,
                          useSpline=TRUE,
                          pInf=NULL,
                          repeats=1,
                          mu_indices=NULL,
                          measurement_indices=NULL,
                          addNoise=TRUE){
    ## Extract parameters
    pars <- parTab$values
    names(pars) <- parTab$names
    if(!is.null(mu_indices)){
        mus <- parTab[parTab$identity == 3,"values"]
        pars <- parTab[parTab$identity == 1,"values"]
        names(pars) <- parTab[parTab$identity==1,"names"]
    }

    if(!is.null(measurement_indices)){
        measurement_bias <- parTab[parTab$identity == 4,"values"]
        pars <- parTab[parTab$identity == 1,"values"]
        names(pars) <- parTab[parTab$identity==1,"names"]
    }
    lambda_pars <- parTab[parTab$identity == 2,"values"]

    ## Create antigenic map for short and long term boosting
    antigenicMap1 <- outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")])
    
    antigenicMapLong <- 1 - pars["sigma1"]*c(antigenicMap1)
    antigenicMapShort <- 1 - pars["sigma2"]*c(antigenicMap1)

    antigenicMapLong[antigenicMapLong < 0] <- 0
    antigenicMapShort[antigenicMapShort < 0] <- 0

    ## Simulate ages
    ages <- floor(runif(n_indiv,ageMin,ageMax))
    
    ## Simulate attack rates
    if(is.null(pInf)){
        if(useSIR){
            pInf <- simulate_ars_buckets(strainIsolationTimes, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
        } else if(useSpline){
            pInf <- simulate_ars_spline(strainIsolationTimes, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"], knots,theta)
        } else {
            pInf <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
        }
    }
    ## Simulate infection histories    
    tmp <- simulate_infection_histories(pInf, simInfPars["logSD"], strainIsolationTimes,
                                        samplingTimes, ages)
    
    infHist <- tmp[[1]]
    ARs <- tmp[[2]]

    ## Simulate titre data
    y <- simulate_group(n_indiv, pars, infHist, strainIsolationTimes, samplingTimes,
                        nsamps, antigenicMapLong,antigenicMapShort,repeats,
                        mus, mu_indices, measurement_bias, measurement_indices,addNoise)

    ## Randomly censor titre values
    y$titre <- y$titre*sample(c(NA,1),nrow(y),prob=c(titreSensoring,1-titreSensoring),replace=TRUE)
    y$run <- 1
    y$group <- 1
    
    DOB <- max(samplingTimes) - ages
    ages <- data.frame("individual"=1:n_indiv, "DOB"=DOB)
    attackRates <- data.frame("year"=strainIsolationTimes,"AR"=ARs)
    return(list(data=y, infectionHistories=infHist, ages=ages, attackRates=attackRates, lambdas=pInf))
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

#' @export
generate_antigenic_map <- function(antigenicDistances, buckets=1){
    ## Following assumptions:
  ## 1. X31 == 1969
  ## 2. PE2009 is like the strain circulating in 2010
  virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
                 "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets
  antigenicDistances$Strain <- virus_key[antigenicDistances$Strain]
  fit <- smooth.spline(antigenicDistances$X,antigenicDistances$Y,spar=0.3)
  x_line <- lm(data = antigenicDistances, X~Strain)
  Strain <- seq(1968*buckets,2016*buckets-1,by=1)
  x_predict <- predict(x_line,data.frame(Strain))
  y_predict <- predict(fit, x=x_predict)
  fit_dat <- data.frame(x=y_predict$x,y=y_predict$y)
  fit_dat$strain <- Strain
  colnames(fit_dat) <- c("x_coord","y_coord","inf_years")
  return(fit_dat)
}

#' @export
euc_distance <- function(i1, i2, fit_dat){
  return(sqrt((fit_dat[i1,"x_coord"] - fit_dat[i2,"x_coord"])^2 + (fit_dat[i1,"y_coord"] - fit_dat[i2,"y_coord"])^2))
}
