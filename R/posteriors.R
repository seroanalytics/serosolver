#' Posterior function pointer
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to pass the data set back and forth.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param titreDat the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param version which version of the posterior function to solve (corresponds mainly to the infection history prior). Mostly just left to 1, but there is one special case where this should be set to 4 for the gibbs sampler. This is only really used by \code{\link{run_MCMC}} to place the infection history prior on the total number of infections across all years and individuals when version = 4
#' @param solve_likelihood usually set to TRUE. If FALSE< does not solve the likelihood and instead just samples/solves based on the model prior
#' @param ageMask see \code{\link{create_age_mask}} - a vector with one entry for each individual specifying the first epoch of circulation in which an individual could have been exposed
#' @param measurement_indices_by_time if not NULL, then use these indices to specify which measurement bias parameter index corresponds to which time
#' @param mu_indices if not NULL, then use these indices to specify which boosting parameter index corresponds to which time
#' @param n_alive if not NULL, uses this as the number alive in a given year rather than calculating from the ages. This is needed if the number of alive individuals is known, but individual birth dates are not
#' @param function_type integer specifying which version of this function to use. Specify 1 to give a posterior solving function; 2 to give the gibbs sampler for infection history proposals; otherwise just solves the titre model and returns predicted titres
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infectionHistories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @export
create_posterior_func <- function(parTab,
                                  titreDat,
                                  antigenicMap,
                                  version=1,
                                  solve_likelihood=TRUE,
                                  ageMask=NULL,
                                  measurement_indices_by_time=NULL,
                                  mu_indices=NULL,
                                  n_alive=NULL,
                                  function_type=1,
                                  ...){
    ## Sort data in same way
    titreDat <- titreDat[order(titreDat$individual, titreDat$run, titreDat$samples, titreDat$virus),]
    
    ## Isolate data table as vectors for speed
    titres <- titreDat$titre

    ## Setup data vectors and extract
    setup_dat <- setup_titredat_for_posterior_func(titreDat, antigenicMap, ageMask, n_alive)

    individuals <- setup_dat$individuals
    antigenicMapMelted <- setup_dat$antigenicMapMelted
    strain_isolation_times <- setup_dat$strain_isolation_times
    infection_strain_indices <- setup_dat$infection_strain_indices
    sample_times <- setup_dat$sample_times
    rows_per_indiv_in_samples <- setup_dat$rows_per_indiv_in_samples
    nrows_per_individual_in_data <- setup_dat$nrows_per_individual_in_data
    cum_nrows_per_individual_in_data <- setup_dat$cum_nrows_per_individual_in_data
    nrows_per_blood_sample <- setup_dat$nrows_per_blood_sample
    measured_strain_indices <- setup_dat$measured_strain_indices
    n_alive <- setup_dat$n_alive
    ageMask <- setup_dat$ageMask
    strainMask <- setup_dat$strainMask
    n_indiv <- setup_dat$n_indiv
    DOBs <- setup_dat$DOBs

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
    
    if (function_type==1) {
        if(solve_likelihood){
                f <- function(pars, infection_history_mat){
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
                    antigenic_map_long <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
                    antigenic_map_short <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])
                    
                    ## Now pass to the C++ function
                    y <- titre_data_group(theta, infection_history_mat, strain_isolation_times,
                                          infection_strain_indices, sample_times,
                                          rows_per_indiv_in_samples,
                                          cum_nrows_per_individual_in_data,
                                          nrows_per_blood_sample,
                                          measured_strain_indices, 
                                          antigenic_map_long, antigenic_map_short,
                                          DOBs, additional_arguments)
                    liks <- r_likelihood(y, titres, theta, expected_indices, measurement_bias)
                    liks <- sum_buckets(liks, nrows_per_individual_in_data)

                    if (explicit_lambda) {
                        liks <- liks + calc_lambda_probs_indiv(lambdas, infection_history_mat, ageMask, strainMask)
                    }
                    if (spline_lambda) {
                        liks <- liks + calc_lambda_probs_monthly(lambdas,  knots, weights,
                                                                 infection_history_mat, ageMask)
                    }          
                    return(liks)
                }
        } else {           
            f <- function(pars, infection_history_mat){
                lambdas <- pars[lambda_indices]          
                theta <- pars[theta_indices]
                weights <- pars[weights_indices]
                knots <- pars[knot_indices]
                mus <- pars[mu_indices_parTab]
                
                liks <- rep(-100000,n_indiv)
                
                if (explicit_lambda) {
                    liks <- liks + calc_lambda_probs_indiv(lambdas, infection_history_mat, ageMask, strainMask)
                }
                if (spline_lambda) {
                    liks <- liks + calc_lambda_probs_monthly(lambda,  knots, weights, infection_history_mat, ageMask)
                }
                
                return(liks)
            }
        }
    } else if (function_type == 2) {
        if (version == 4) {
            n_alive_total <- sum(n_alive)
        } else {
            n_alive_total <- -1
        }
        
        ## Gibbs proposal on infection histories
        f <- function(pars, infection_history_mat,
                      alpha, beta,
                      indivPropn, nYears,
                      swapPropn=0.5,swapDistance=1, temp=1){
            lambdas <- pars[lambda_indices]
            theta <- pars[theta_indices]
            weights <- pars[weights_indices]
            knots <- pars[knot_indices]
            mus <- pars[mu_indices_parTab]
            names(theta) <- parNames_theta
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_parTab]
                titre_shifts <- measurement_bias[expected_indices]
            }

            if (use_strain_dependent) {
                additional_arguments[["mus"]] <- mus
            }
            names(theta) <- parNames_theta
            ## Work out short and long term boosting cross reactivity - C++ function
            antigenic_map_long <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])
            ## Now pass to the C++ function
            new_infection_history_mat <- infection_history_proposal_gibbs(theta,
                                                                       infection_history_mat,
                                                                       indivPropn,
                                                                       nYears,
                                                                       ageMask,
                                                                       strainMask,
                                                                       n_alive,
                                                                       swapPropn,
                                                                       swapDistance,
                                                                       alpha,
                                                                       beta,
                                                                       strain_isolation_times,
                                                                       infection_strain_indices,
                                                                       sample_times,
                                                                       rows_per_indiv_in_samples,
                                                                       cum_nrows_per_individual_in_data,
                                                                       nrows_per_blood_sample,
                                                                       measured_strain_indices, 
                                                                       antigenic_map_long,
                                                                       antigenic_map_short,
                                                                       titres,
                                                                       titre_shifts,
                                                                       additional_arguments,
                                                                       DOBs,
                                                                       solve_likelihood,
                                                                       n_alive_total,
                                                                       temp
                                                                       )
            return(new_infection_history_mat)
        }
    } else {
        f <- function(pars, infection_history_mat){
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
            antigenic_map_long <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenicMapMelted, theta["sigma2"])
            
            ## Now pass to the C++ function
            y <- titre_data_group(theta, infection_history_mat, strain_isolation_times,
                                  infection_strain_indices, sample_times,
                                  rows_per_indiv_in_samples,
                                  cum_nrows_per_individual_in_data,
                                  nrows_per_blood_sample,
                                  measured_strain_indices, 
                                  antigenic_map_long, antigenic_map_short,
                                  DOBs, additional_arguments)
            
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
