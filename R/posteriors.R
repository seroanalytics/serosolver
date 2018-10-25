#' Posterior function pointer
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to pass the data set back and forth.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param titreDat the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param PRIOR_FUNC user function of prior for model parameters. Should take parameter values only
#' @param version integer specifying which version of the posterior function to use. This is mainly used for turning off or on the likelihood for debugging the implicit prior. The versions are as follows: 1) Likelihood turned on and no explicit prior. Any prior is implicit from the proposal distribution; 2) No likelihood and no explicit prior; 3) Explicit prior term (as specified by PRIOR_FUNC) and no likelihood. Any implicit prior from the proposal will also be included; 4) Version of the posterior distribution where we infer an explicit force of infection for each epoch; 5) Explicit likelihood and prior - any implicit prior from the proposal will also be included, so proposals should be symmetric; 6) just the negative log likelihood, but function only expects vector of pars (not infectionHistories), which must be specified with ...; 7+) Returns the predicted titres (model solving function); 99) Returns a gibbs-sampled infection history matrix
#' @param ageMask see \code{\link{create_age_mask}} - a vector with one entry for each individual specifying the first epoch of circulation in which an individual could have been exposed
#' @param temp temperature for parallel tempering
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infectionHistories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @export
create_posterior_func <- function(parTab,
                                  titreDat,
                                  antigenicMap,
                                  version=1,
                                  ageMask=NULL,
                                  measurement_indices_by_time=NULL,
                                  mu_indices=NULL,
                                  n_alive=NULL,
                                  ...){
    ## Sort data in same way
    titreDat <- titreDat[order(titreDat$individual, titreDat$run, titreDat$samples, titreDat$virus),]
    
    ## Isolate data table as vectors for speed
    titres <- titreDat$titre
    
    ## The entry of each virus in the data corresponding to the antigenic map
    ## Note 0 indexed for C++ code
    virusIndices <- match(titreDat$virus, antigenicMap$inf_years) - 1
    ## Unique strains that an individual could be exposed to
    strains <- unique(antigenicMap$inf_years)
    n_strains <- length(strains)
    
    ## The entry of each strain in the antigenic map table
    strainIndices <- match(strains, strains) - 1

    ## Generate distances between each pair of viruses from the antigenic map
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))

    ## Get unique measurement sets for each individual at
    ## each sampling time for each repeat
    samples <- unique(titreDat[,c("individual","samples","run")])
    #samples <- samples[order(samples$individual, samples$run, samples$samples),]

    ## Extract vector of sampling times and individual indices for speed
    sampleTimes <- samples$samples
    individuals <- samples$individual  
    n_indiv <- length(unique(samples$individual))

    ## For each individual, sample and run, get the corresponding number of rows from the overall
    ## data matrix. 
    indicesData <- NULL
    for(i in 1:nrow(samples)){
        indicesData <- c(indicesData, nrow(samples[titreDat$individual == samples[i,"individual"] &
                                                   titreDat$samples == samples[i,"samples"] &
                                                   titreDat$run == samples[i,"run"],]))
    }
    
    ## Get indexing for individuals. This works out which rows in the titre data correspond
    ## to which individuals.
    ## Firstly, how many runs/samples correspond to each individual?
    indicesSamples <- c(0)
    for(individual in unique(individuals)){
        indicesSamples <- c(indicesSamples, length(individuals[individuals==individual]))
    }
    indicesSamples <- cumsum(indicesSamples)
    
    ## How many rows in the titre data correspond to each individual?
    indicesDataOverall <- NULL
    for(individual in unique(individuals)){
        indicesDataOverall <- c(indicesDataOverall, nrow(titreDat[titreDat$individual == individual,]))
    }
    indicesDataOverall <- cumsum(c(0,indicesDataOverall))    
    indicesOverallDiff <- diff(indicesDataOverall)

    if(!is.null(titreDat$DOB)){
        DOBs <- unique(titreDat[,c("individual","DOB")])[,2]
    } else {
        DOBs <- rep(min(strains), n_indiv)
    }
    if(is.null(ageMask)){
        if(!is.null(titreDat$DOB)){
            ageMask <- create_age_mask(DOBs, strains)
        } else {
            ageMask <- rep(1, n_indiv)
        }
    }
    ageMask <- create_age_mask(DOBs, strains)
    strainMask <- create_strain_mask(titreDat,strains)
    masks <- data.frame(cbind(ageMask, strainMask))
    if (is.null(n_alive)) {
        n_alive <- sapply(seq(1,length(strains)), function(x)
            nrow(masks[masks$ageMask <=x & masks$strainMask >= x,]))
    }

#########################################################
    ## Extract parameter type indices from parTab, to split up
    ## similar parameters in model solving functions
    option_indices <- which(parTab$type == 0)
    theta_indices <- which(parTab$type %in% c(0,1))
    lambda_indices <- which(parTab$type == 2)
    measurement_indices_parTab <- which(parTab$type == 3)
    weights_indices <- which(parTab$type == 4) ## For functional form version
    knot_indices <- which(parTab$type == 5)
    mu_indices_parTab <- which(parTab$type == 6)
#########################################################
    
    parNames_theta <- parTab[theta_indices,"names"]

    ## Find which options are being used in advance for speed
    explicit_lambda <- (length(lambda_indices) > 0)
    spline_lambda <- (length(knot_indices) > 0)
    use_measurement_bias <- (length(measurement_indices_parTab) > 0) & !is.null(measurement_indices_by_time)
    titre_shifts <- NULL
    expected_indices <- NULL
    measurement_bias <- NULL
    use_strain_dependent <- (length(mu_indices) > 0) & !is.null(mu_indices)
    additional_arguments <- NULL
    
    if (use_measurement_bias) {
        expected_indices <- measurement_indices_by_time[match(titreDat$virus,strains)]
    }

    if (use_strain_dependent) {
        additional_arguments <- list("boosting_vec_indices"=mu_indices-1,
                                     "mus"=rep(2,length(strains)))
    }
    
    if (version==1) {
        f <- function(pars, infectionHistories){
            lambdas <- pars[lambda_indices]
            theta <- pars[theta_indices]
            weights <- pars[weights_indices]
            knots <- pars[knot_indices]
            mus <- pars[mu_indices_parTab]
            
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_parTab]
                to_add <- measurement_bias[expected_indices]
            }

            if (use_strain_dependent) {
                additional_arguments[["mus"]] <- mus
            }
            names(theta) <- parNames_theta
            ## Work out short and long term boosting cross reactivity - C++ function
            antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
            antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])
            
            ## Now pass to the C++ function
            y <- titre_data_group(theta, infectionHistories, strains, strainIndices, sampleTimes,
                                  indicesData,indicesDataOverall,indicesSamples, virusIndices, 
                                  antigenicMapLong, antigenicMapShort, DOBs, additional_arguments)
            liks <- r_likelihood(y, titres, theta, expected_indices, measurement_bias)
            liks <- sum_buckets(liks, indicesOverallDiff)
            
            if (explicit_lambda) {
                liks <- liks + calc_lambda_probs_indiv(lambdas, infectionHistories, ageMask, strainMask)
            }
            if (spline_lambda) {
                liks <- liks + calc_lambda_probs_monthly(lambdas,  knots, weights, infectionHistories, ageMask)
            }          
            return(liks)
        }
    } else if (version == 2) {
        f <- function(pars, infectionHistories){
            lambdas <- pars[lambda_indices]          
            theta <- pars[theta_indices]
            weights <- pars[weights_indices]
            knots <- pars[knot_indices]
            mus <- pars[mu_indices_parTab]
            
            liks <- rep(-100000,n_indiv)
            
            if (explicit_lambda) {
                liks <- liks + calc_lambda_probs_indiv(lambdas, infectionHistories, ageMask, strainMask)
            }
            if (spline_lambda) {
                liks <- liks + calc_lambda_probs_monthly(lambda,  knots, weights, infectionHistories, ageMask)
            }
            
            return(liks)
        }
    } else if (version == 99) {
        ## Gibbs proposal on infection histories
        f <- function(pars, infectionHistories, alpha, beta, indivPropn, nYears,
                      swapPropn=0.5,swapDistance=1, temp=1){
            lambdas <- pars[lambda_indices]
            theta <- pars[theta_indices]
            weights <- pars[weights_indices]
            knots <- pars[knot_indices]
            mus <- pars[mu_indices_parTab]
            names(theta) <- parNames_theta
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_parTab]
                #print(measurement_indices_parTab)
                #print(parTab$names[measurement_indices_parTab])
                #print(measurement_bias)
                titre_shifts <- measurement_bias[expected_indices]
            }

            if (use_strain_dependent) {
                additional_arguments[["mus"]] <- mus
            }
            names(theta) <- parNames_theta
            ## Work out short and long term boosting cross reactivity - C++ function
            antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
            antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])
            ## Now pass to the C++ function
            #print(tracemem(titre_shifts))
            #for(i in 2:length(indicesDataOverall)){
            #    print(tracemem(titre_shifts[indicesDataOverall[i-1]:indicesDataOverall[i]]))
            #    
            #}
            new_infectionHistories <- infection_history_proposal_gibbs(theta,
                                                                       infectionHistories,
                                                                       indivPropn,
                                                                       nYears,
                                                                       ageMask,
                                                                       strainMask,
                                                                       n_alive,
                                                                       swapPropn,
                                                                       swapDistance,
                                                                       alpha,
                                                                       beta,
                                                                       strains,
                                                                       strainIndices,
                                                                       sampleTimes,
                                                                       indicesData,
                                                                       indicesDataOverall,
                                                                       indicesSamples,
                                                                       virusIndices, 
                                                                       antigenicMapLong,
                                                                       antigenicMapShort,
                                                                       titres,
                                                                       titre_shifts,
                                                                       additional_arguments,
                                                                       DOBs,
                                                                       temp
                                                                       )
            return(new_infectionHistories)
        }
    } else {
        print("Model solving function")
        f <- function(pars, infectionHistories){
            lambdas <- pars[lambda_indices]
            theta <- pars[theta_indices]
            weights <- pars[weights_indices]
            knots <- pars[knot_indices]
            mus <- pars[mu_indices_parTab]
            
            if (use_strain_dependent) {
                additional_arguments[["mus"]] <- mus
            }
            
            names(theta) <- parNames_theta

            ## Work out short and long term boosting cross reactivity - C++ function
            antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
            antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])
            
            ## Now pass to the C++ function
            y <- titre_data_group(theta, infectionHistories, strains, strainIndices, sampleTimes,
                                  indicesData,indicesDataOverall,indicesSamples, virusIndices, 
                                  antigenicMapLong, antigenicMapShort, DOBs, additional_arguments)
            
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_parTab]
                titre_shifts <- measurement_bias[expected_indices]
                y <- y + titre_shifts
            }
            
            return(y)
        }
    }
    f
}

#' @export
r_likelihood <- function(expected, data, theta, expected_indices=NULL, measurement_shifts=NULL){
    if(!is.null(expected_indices) & !is.null(measurement_shifts)){
        expected <- expected + measurement_shifts[expected_indices]
    }
    
    liks <- numeric(length(expected))
    largeI <- data > theta["MAX_TITRE"]
    smallI <- data <= 0
    restI <- data > 0 & data <= theta["MAX_TITRE"]
    
    liks[largeI] <- pnorm(theta["MAX_TITRE"], expected[largeI],theta["error"],lower.tail=FALSE,log.p=TRUE)
    liks[smallI] <- pnorm(1, expected[smallI],theta["error"],lower.tail=TRUE,log.p=TRUE)
    liks[restI] <- log(pnorm(data[restI]+1,expected[restI],theta["error"],lower.tail=TRUE,log.p=FALSE) - 
                       pnorm(data[restI],expected[restI], theta["error"],lower.tail=TRUE,log.p=FALSE))
    return(liks)
}


#' Calculate FOI log probability
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param lambdas a vector of FOIs
#' @param infHist the matrix of infection histories
#' @param ageMask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a single log probability
#' @export
calc_lambda_probs <- function(lambdas, infHist, ageMask, strainMask){
    lik <- 0
    for(i in 1:ncol(infHist)){
        use_indivs <- intersect(which(ageMask <= i), which(strainMask >= i))
        lik <- lik + sum(log((lambdas[i]^infHist[use_indivs,i] *
                              (1-lambdas[i])^(1-infHist[use_indivs,i]))))
    }
    lik
}

#' Calculate FOI log probability vector
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param lambdas a vector of FOIs
#' @param infHist the matrix of infection histories
#' @param ageMask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a vector of log probabilities for each individual
#' @export
calc_lambda_probs_indiv <- function(lambdas, infHist, ageMask, strainMask){
    lik <- numeric(nrow(infHist))
    for(i in 1:ncol(infHist)){
        lik <- lik + log(((lambdas[i]^infHist[,i]) * (1-lambdas[i])^(1-infHist[,i])))* as.numeric(ageMask <= i)*as.numeric(strainMask >= i)
    }
    lik
}

#' @export
calc_lambda_probs_monthly <- function(foi, knots, theta, infHist, ageMask){
    lambdas <- generate_lambdas(foi, knots, theta, length(foi), 12)
    lik <- numeric(nrow(infHist))
    for(i in 1:ncol(infHist)){
        #lik <- lik + log(((lambdas[i]^infHist[,i]) * (1-lambdas[i])^(1-infHist[,i]))) + as.numeric(ageMask <= i)
        lik <- lik + infHist[,i]*log(lambdas[i]) + (1-infHist[,i])*log(lambdas[i]) + log(as.numeric(ageMask <= i))
    }
    lik
}


#' FOR DEBUGGING
#'
#' Brute force implementation of calculating the explicit FOI
#' @export
calc_lambda_probs_indiv_brute<- function(lambdas, infHist, ageMask){
  lik <- numeric(nrow(infHist))
  for(j in 1:nrow(infHist)){
      lik[j] <- 0
      age <- ageMask[j]
      for(i in 1:ncol(infHist)){
          if(i >= age){
              lik[j] <- lik[j] + log(lambdas[i]^infHist[j,i] * (1-lambdas[i])^(1-infHist[j,i]))
          }
      }
  }
  lik
}


#' Generate FOI lambdas
#'
#' @export
generate_lambdas <- function(foi, knots, theta, nYears, buckets, degree=2){
    x <- seq(0,buckets-1,by=1)/buckets
    nKnots <- length(knots) + degree + 1
    allDat <- NULL
    index <- 1
    tmp <- genSpline_y(x, knots, degree, theta)
    tmp <- tmp/sum(tmp)
    allDat <- numeric(length(tmp)*length(foi))
    for(i in 1:nYears){        
        allDat[index:(index + length(tmp)-1)] <-  tmp*foi[i]
        index <- index + length(tmp)    
    }
    allDat
}

#' @export
genSpline_y <- function(x, knots, degree, theta, intercept=TRUE) {
  
    basis <- bs(x = x, knots = knots, degree = degree,
              Boundary.knots = c(0,1), intercept = intercept)
    
    y.spline <- basis %*% theta
    return(as.vector(y.spline))
}


#' @export
prob_shifts <- function(rhos, pars){
    rho_mean <- pars["rho_mean"]
    rho_sd <- pars["rho_sd"]
    ## l_mean <- log(mu_mean) - (mu_sd^2)/2
    return(sum(dnorm(rhos, rho_mean, rho_sd,log=TRUE)))
    #return(sum(dlnorm(mus, l_mean, mu_sd, log=TRUE)))
}

#' @export
create_prob_shifts <- function(parTab){
    parNames <- parTab$names
    rho_indices <- which(parTab$type == 3)
    
    f <- function(pars){
        names(pars) <- parNames
        rhos <- pars[rho_indices]
        rho_mean <- pars["rho_mean"]
        rho_sd <- pars["rho_sd"]
        return(sum(dnorm(rhos, rho_mean, rho_sd,log=TRUE)))
    }
    f
}

create_prior_mu <- function(parTab){
    ## Extract parameter type indices from parTab, to split up
    ## similar parameters in model solving functions    
    option_indices <- which(parTab$type == 0)
    theta_indices <- which(parTab$type %in% c(0,1))
    lambda_indices <- which(parTab$type == 2)
    measurement_indices_parTab <- which(parTab$type == 3)
    weights_indices <- which(parTab$type == 4) ## For functional form version
    knot_indices <- which(parTab$type == 5)
    mu_indices_parTab <- which(parTab$type == 6)

    parNames_theta <- parTab[theta_indices,"names"]
    
    f <- function(pars){
        mus <- pars[mu_indices_parTab]
        pars <- pars[theta_indices]
        names(pars) <- parNames_theta
        return(prob_mus(mus, pars))
    }
}

#' @export
prob_mus <- function(mus, pars){
    #return(0)
    mu_mean <- pars["mu_mean"]
    mu_sd <- pars["mu_sd"]
    return(sum(dnorm(mus, mu_mean,mu_sd,log=TRUE)))
    #location <- log(mu_mean^2 / sqrt(mu_sd^2 + mu_mean^2))
    #shape <- sqrt(log(1 + (mu_sd^2/mu_mean^2)))
    #l_mean <- log(mu_mean) - (mu_sd^2)/2
    #p <- sum(dnorm(log(mus),mu_mean,mu_sd,log=TRUE))
    #p_mean <- 0.6
    #p_sd <- 0.5
    #p_mu <- log(p_mean/sqrt(1 + (p_sd/p_mean)^2))
    #p_sigma <- sqrt(log(1 + (p_sd/p_mean)^2))

    #p_lik <- log(p_sigma*2.506628) - 0.5*((mu_mean - p_mu)/p_sigma)^2

    #return(p+p_lik)
    #return(sum(log(dtruncnorm(mus, a=0,mean=mu_mean, sd=mu_sd))))
    #return(sum(dlnorm(mus, location, shape, log=TRUE)))
    #mean_log_y <- mean(log(mus))
    #sd_log_y <- sd(log(mus))
    #sigmaOfLogY <- dunif(mu_sd, 0.001*sd_log_y,1000*sd_log_y)
    #muOfLogY <- dnorm(mu_mean, mean_log_y, 1/(10*sd_log_y)^2)
    #return(sum(dlnorm(mus, mu_mean, 1/mu_sd^2, log=TRUE)) + sigmaOfLogY + muOfLogY)
     return(sum(dlnorm(mus, mu_mean, mu_sd, log=TRUE)))
} 
