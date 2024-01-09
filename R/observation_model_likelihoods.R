#' Likelihood function given discrete data (normal)
#'
#' Calculates the likelihood of observing a set of discrete measurements given a corresponding set of predicted antibody levels
#' @param expected vector of expected antibody levels
#' @param data vector of observed discrete measurements
#' @param expected_indices the indices of the measurement_shifts vector that each predicted antibody level needs adding to it
#' @param measurement_shifts the vector of measurement shifts for each cluster to add to the predicted antibody levels
#' @return a vector with the likelihood of making each observation given the predictions
#' @export
r_likelihood <- function(expected, data, theta, expected_indices = NULL, measurement_shifts = NULL) {
  if (!is.null(expected_indices) & !is.null(measurement_shifts)) {
    expected <- expected + measurement_shifts[expected_indices]
  }
  
  ## Vectorise, calculate boundaries separately
  liks <- numeric(length(expected))
  large_i <- data >= theta["max_titre"]
  small_i <- data < 1
  rest_i <- data >= 1 & data < theta["max_titre"]
  
  liks[large_i] <- pnorm(theta["max_titre"], expected[large_i], theta["obs_sd"], lower.tail = FALSE, log.p = TRUE)
  liks[small_i] <- pnorm(1, expected[small_i], theta["obs_sd"], lower.tail = TRUE, log.p = TRUE)
  liks[rest_i] <- log(pnorm(data[rest_i] + 1, expected[rest_i], theta["obs_sd"], lower.tail = TRUE, log.p = FALSE) -
                        pnorm(data[rest_i], expected[rest_i], theta["obs_sd"], lower.tail = TRUE, log.p = FALSE))
  return(liks)
}


#' Likelihood function given continuous data (normal)
#'
#' Calculates the likelihood of observing a set of continuous and bounded measurements given a corresponding set of predicted antibody levels
#' @param expected vector of expected antibody levels
#' @param data vector of observed continuous measurements
#' @param expected_indices the indices of the measurement_shifts vector that each predicted antibody levels needs adding to it
#' @param measurement_shifts the vector of measurement shifts for each cluster to add to the predicted antibody levels
#' @return a vector with the likelihood of making each observation given the predictions
#' @export
r_likelihood_continuous <- function(expected, data, theta, expected_indices = NULL, measurement_shifts = NULL) {
  if (!is.null(expected_indices) & !is.null(measurement_shifts)) {
    expected <- expected + measurement_shifts[expected_indices]
  }
  
  ## Vectorise, calculate boundaries separately
  liks <- numeric(length(expected))
  large_i <- data >= theta["max_titre"]
  small_i <- data <= theta["min_titre"]
  rest_i <- data > theta["min_titre"] & data < theta["max_titre"]
  
  liks[large_i] <- pnorm(theta["max_titre"], expected[large_i], theta["obs_sd"], lower.tail = FALSE, log.p = TRUE)
  liks[small_i] <- pnorm(theta["min_titre"], expected[small_i], theta["obs_sd"], lower.tail = TRUE, log.p = TRUE)
  liks[rest_i] <- dnorm(data[rest_i], expected[rest_i], theta["obs_sd"], log = TRUE)
  return(liks)
}


#' Prior on measurement shifts
#'
#' Assumes measurement shifts are drawn from a normal distribution (random effects) with given standard deviation and mean. Code is commented out to assume log normally distributied
#' @param rhos vector of measurement shifts
#' @param pars vector of parameters, including rho_mean and rho_sd for the normal distribution
#' @return a single prior probability
#' @family priors
#' @export
prob_shifts <- function(rhos, pars) {
  rho_mean <- pars["rho_mean"]
  rho_sd <- pars["rho_sd"]
  return(sum(dnorm(rhos, rho_mean, rho_sd, log = TRUE)))
}


#' Measurement shift creation
#'
#' Creates a function to calculate the prior probability of a set of measurement shifts. Assumes normally distributed random effects
#' @param par_tab the parameter table as in \code{\link{create_posterior_func}}
#' @return a function pointer to solve the measurement shifts prior
#' @family priors
#' @export
create_prob_shifts <- function(par_tab) {
  par_names <- par_tab$names
  rho_indices <- which(par_tab$type == 3)
  
  f <- function(pars) {
    names(pars) <- par_names
    rhos <- pars[rho_indices]
    rho_mean <- pars["rho_mean"]
    rho_sd <- pars["rho_sd"]
    return(sum(dnorm(rhos, rho_mean, rho_sd, log = TRUE)))
  }
  f
}