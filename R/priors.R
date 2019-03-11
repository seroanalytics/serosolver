#' Beta binomial infection history prior
#'
#' Calculates the beta binomial infection history prior for each individual, taking into account the number of years they could be infected.
#' @param pars the named parameter vector with alpha and beta
#' @param infection_history the infection hsitory matrix
#' @param age_mask the age mask, giving the first index of the infection_history matrix that each individual can be exposed to. One entry per individual
#' @return a prior probability for each individual
#' @family priors
#' @export
infection_history_prior <- function(pars, infection_history, age_mask) {
  N <- ncol(infection_history) - age_mask + 1
  a <- pars["alpha"]
  b <- pars["beta"]
  priors <- numeric(nrow(infection_history))
  for (i in 1:length(priors)) {
    priors[i] <- log(dbb_prior(sum(infection_history[i, ]), N[i], a, b))
  }
  return(priors)
}

#' Beta binomial density
#' @export
density_beta_binom <- function(x, N, u, v) {
  (beta(x + u, N - x + v) / beta(u, v)) * choose(N, x)
}

#' Beta binomial prior used here (no choose constant)
#' @export
dbb_prior <- function(x, N, u, v) {
  (beta(x + u, N - x + v) / beta(u, v))
}

#' Beta prior on an infection
db <- function(x, a, b) {
  x^(a - 1) * (1 - x^(b - 1)) / beta(a, b)
}

#' Infection history prior R
#'
#' R implementation of the infection history prior assuming common infection risk across individuals in a given year
#' @param infection_history the infection history matrix
#' @param age_mask the age mask, giving the first index of the infection_history matrix that each individual can be exposed to. One entry per individual
#' @param alpha alpha parameter for beta distribution
#' @param beta beta parameter for beta distribution
#' @return a prior probability for each individual
#' @family priors
#' @export
inf_mat_prior <- function(infection_history, age_mask, alpha, beta1) {
  n_alive <- sapply(1:ncol(infection_history), function(x) length(age_mask[age_mask <= x]))
  lk <- 0
  for (i in 1:length(n_alive)) {
    lk <- lk + log(beta(sum(infection_history[, i]) + alpha, n_alive[i] - sum(infection_history[, i]) + beta1) / beta(alpha, beta1))
  }
  return(lk)
}

#' Fit beta distribution to MCMC output
#'
#' Attempts to fit a beta distribution to a data frame of MCMC output from a previous run.
#' @param chain_samples the MCMC chain data frame to be fit to
#' @param par_name the column label to fit to
#' @param error_tol = 9999999999999, what's the error tolerance on the fit? Might take some tweaking
#' @param try_attempts = 10 how many fitting attempts to try before giving up
#' @param plot_fit = TRUE, if TRUE, plots the fit to the MCMC chain
#' @return the model fit object as returned by optim
#' @seealso \code{\link{fit_normal_prior}}
fit_beta_prior <- function(chain_samples, par_name = "", error_tol = 999999999, try_attempts = 10, plot_fit = TRUE) {
  data <- density(chain_samples)
  data <- data.frame(x = data$x, y = data$y)
  x <- data$x

  f <- function(pars) {
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dbeta(x, shape1, shape2)
    return(sum((out - data$y)^2))
  }
  res <- list(value = error_tol + 1)
  try_no <- 1
  while (res$value > error_tol & try_no < try_attempts) {
    message(cat("Try no: ", try_no))
    res <- optim(runif(2, 1, 100), f, method = "Nelder-Mead", control = list(abstol = 1e-8, reltol = 1e-8))
    try_no <- try_no + 1
  }
  if (try_no == try_attempts) {
    message("Fitting did not work - does this distribution make sense?")
    return(list(par = c(NA, NA)))
  }
  if (plot_fit) {
    plot(dbeta(x, res$par[1], res$par[2]) ~ x,
      type = "l", col = "blue",
      main = paste0("Beta dist fit to posterior of ", par_name),
      xlab = "Value",
      ylab = "Density"
    )
    legend(x = x[1], y = max(data$y), box.lty = 0, lty = 1, col = c("blue", "red"), legend = c("Model fit", "Posterior density"))
    lines(data, col = "red")
  }
  return(res)
}
#' Fit normal distribution to MCMC output
#'
#' Attempts to fit a normal distribution to a data frame of MCMC output from a previous run.
#' @param chain_samples the MCMC chain data frame to be fit to
#' @param par_name the column label to fit to
#' @param error_tol = 9999999999999, what's the error tolerance on the fit? Might take some tweaking
#' @param try_attempts = 10 how many fitting attempts to try before giving up
#' @param plot_fit = TRUE, if TRUE, plots the fit to the MCMC chain
#' @return the model fit object as returned by optim
#' @seealso \code{\link{fit_normal_prior}}
fit_normal_prior <- function(chain_samples, par_name = "", error_tol = 999999999, try_attempts = 10, plot_fit = TRUE) {
  data <- density(chain_samples)
  data <- data.frame(x = data$x, y = data$y)
  x <- data$x

  f <- function(pars) {
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dnorm(x, mean = shape1, sd = shape2)
    return(sum((out - data$y)^2))
  }
  start_mean <- mean(chain_samples)
  start_sd <- sd(chain_samples)
  res <- list(value = error_tol + 1)
  try_no <- 1
  while (res$value > error_tol & try_no < try_attempts) {
    message(cat("Try no: ", try_no))
    res <- optim(c(start_mean, start_sd), f, method = "Nelder-Mead", control = list(abstol = 1e-8, reltol = 1e-8))
    try_no <- try_no + 1
  }
  if (try_no == try_attempts) {
    message("Fitting did not work - does this distribution make sense?")
    return(list(par = c(NA, NA)))
  }
  if (plot_fit) {
    plot(dnorm(x, res$par[1], res$par[2]) ~ x,
      type = "l", col = "blue",
      main = paste0("Normal dist fit to posterior of ", par_name),
      xlab = "Value",
      ylab = "Density"
    )
    legend(x = x[1], y = max(data$y), box.lty = 0, lty = 1, col = c("blue", "red"), legend = c("Model fit", "Posterior density"))
    lines(data, col = "red")
  }
  return(res)
}

#' Find beta parameters for mean and variance
#'
#' Finds the alpha and beta parameters for a beta distribution that gives the desired mean and variance
#' @param mean the mean of the beta distribution (between 0 and 1)
#' @param var the variance of the beta distribution (this is likely going to be somewhere less than 0.25)
#' @param make_plot if TRUE, plots the resulting beta distibution to the R device
#' @return a list with alpha and beta
#' @export
find_beta_parameters <- function(mean, var, make_plot = FALSE) {
  if(mean < 0 | mean > 1) stop("Mean is outside of bounds (0, 1)")
  if(var < 0 | var > 0.25) stop("Var is outside of bounds (0, 0.25^2)")
  alpha <- ((1 - mean) / (var) - (1 / mean)) * mean^2
  beta <- alpha * (1 / mean - 1)
  x <- seq(0, 1, by = 0.001)
  y <- dbeta(x, alpha, beta)
  if (make_plot) plot(y ~ x, type = "l")
  return(list(alpha = alpha, beta = beta))
}

#' Find beta parameters with maximum variance
#'
#' Given a desired annual mean attack rate and the number of epochs to consider per year, gives the alpha and beta parameters with the desired mean but maximum variance
#' @param desired_annual_mean the desired ANNUAL mean attack rate
#' @param the number of buckets to split each year into
#' @return a list with alpha and beta
#' @export
generate_beta_prior_pars <- function(desired_annual_mean, buckets) {
  mean1 <- desired_annual_mean / buckets
  max_var <- mean1 * (1 - mean1) - 0.000001
  pars <- find_beta_parameters(mean1, max_var)
  alpha <- pars$alpha
  beta <- pars$beta
  y <- dbeta(seq(0.01, 1 - 0.01, by = 0.01), alpha, beta, log = FALSE)
  while (!(all(y == cummin(y)))) {
    max_var <- max_var - 0.001
    pars <- find_beta_parameters(mean1, max_var)
    alpha <- pars$alpha
    beta <- pars$beta
    y <- dbeta(seq(0.01, 1 - 0.01, by = 0.01), alpha, beta, log = FALSE)
  }
  return(pars)
}
