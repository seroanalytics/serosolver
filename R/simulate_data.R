#' @export
simulate_group <- function(n_indiv, theta, infectionHistories,
                           samples, strainIsolationTimes,
                           antigenicMapLong, antigenicMapShort){
    dat <- NULL
    for(i in 1:n_indiv){
        y <- as.data.frame(simulate_individual(theta, infectionHistories[i,],
                                               samples[samples$individual == i,"present"],
                                               unique(samples$sample), strainIsolationTimes, antigenicMapLong,
                                               antigenicMapShort))
        y$indiv <- i
        colnames(y) <- c("sample","strain","titre","individual")
        dat <- rbind(dat, y[,c("individual","sample","strain","titre")])
    }
    return(dat)
    
}

#' @export
simulate_individual <- function(theta, infectionHistory,
                                samples, samplingTimes,
                                strainIsolationTimes,
                                antigenicMapLong, antigenicMapShort){
    dat <- matrix(ncol=3, nrow=length(samplingTimes)*length(infectionHistory))
    titres <- NULL
    dates <- NULL
    for(i in 1:length(samplingTimes)){
        if(samples[i] > 0){
            y <- infection_model_indiv(theta, infectionHistory, samplingTimes[i],
                                       strainIsolationTimes, antigenicMapLong, antigenicMapShort)
            y <- add_noise(y, theta)
            titres <- c(titres, y)
        } else {
            titres <- c(titres, rep(NA, length(strainIsolationTimes)))
        }
        dates <- c(dates, rep(samplingTimes[i],length(strainIsolationTimes)))
    }
    dat[,1] <- dates
    dat[,2] <- strainIsolationTimes
    dat[,3] <- titres
    return(dat)
}
#' @export
add_noise <- function(y, theta){
    noise_y <- floor(rnorm(length(y), mean=y, sd=theta["error"]))
    noise_y[noise_y < 0] <- 0
    noise_y[noise_y > theta["MAX_TITRE"]] <- theta["MAX_TITRE"]
    return(noise_y)    
}

#' @export
simulate_attack_rates <- function(infectionYears,meanPar=0.15,sdPar=0.5, 
                                  largeFirstYear=FALSE,bigYearMean=0.5){
    attack_year <- rlnorm(infectionYears, meanlog=log(meanPar)-sdPar^2/2,sdlog=sdPar)
    if(largeFirstYear) attack_year[1] <- rlnorm(1,meanlog=log(bigYearMean)-(sdPar/2)^2/2,sdlog=sdPar/2)
    return(attack_year)
}
#' @export
simulate_infection_histories <- function(pInf, infSD, strainIsolationTimes, samplingTimes, ages){
    n_strains <- length(pInf)
    n_indiv <- length(ages)
    indivs <- 1:n_indiv
    infectionHistories <- matrix(0,ncol=n_strains,nrow=n_indiv)
    attackRates <- rlnorm(length(pInf), meanlog=log(pInf)- infSD^2/2,sdlog=infSD)

    ## Should this be necessary?
    attackRates[attackRates > 1] <- 1
    for(i in 1:n_strains){
        alive <- (max(samplingTimes) - ages) <= strainIsolationTimes[i]
        infectionHistories[sample(indivs[alive], round(length(indivs[alive])*attackRates[i])),i] <- 1
    }
    return(infectionHistories)    
}
#' @export
simulate_data <- function(parTab, group=1,n_indiv,
                          strainIsolationTimes, samplingTimes,
                          antigenicMap,
                          sampleSensoring=0, titreSensoring=0,
                          ageMin=5,ageMax=80,
                          simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1)){
    pars <- parTab$values
    names(pars) <- parTab$names

    antigenicMap1 <- outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")])
    
    ## Create matrix to describe which individuals had which samples
    samples <- expand.grid(group=group, individual=1:n_indiv, samples=samplingTimes, present=1)
    samples$present <- sample(c(0,1),nrow(samples),prob=c(sampleSensoring,1-sampleSensoring),replace=TRUE)

    antigenicMapLong <- 1 - pars["sigma1"]*c(antigenicMap1)
    antigenicMapShort <- 1 - pars["sigma2"]*c(antigenicMap1)

    antigenicMapLong[antigenicMapLong < 0] <- 0
    antigenicMapShort[antigenicMapShort < 0] <- 0
    
    ages <- floor(runif(n_indiv,ageMin,ageMax))

    pInf <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
    infHist <- simulate_infection_histories(pInf, simInfPars["logSD"], strainIsolationTimes, samplingTimes, ages)

    y <- simulate_group(n_indiv, pars, infHist, samples, strainIsolationTimes, antigenicMapLong,antigenicMapShort)
    
    y$titre <- y$titre*sample(c(0,1),nrow(y),prob=c(titreSensoring,1-titreSensoring),replace=TRUE)
    y <- cbind(group=group,y)
    
    return(list(data=y,samples=samples, infectionHistories=infHist, ages=ages))
}

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
