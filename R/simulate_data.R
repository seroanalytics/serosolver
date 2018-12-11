#' Simulate group data
#'
#' Simulates a full set of titre data for n_indiv individuals with known theta and infection_histories. Each individual gets nsamps random samples from sampleTimes, and infections can occur at any of strain_isolation_times
#' @param n_indiv the number of individuals to simulate
#' @param theta the named parameter vector
#' @param infection_histories the matrix of 1s and 0s giving presence/absence of infections for each individual
#' @param strain_isolation_times the vector of strain circulation times (ie. possible infection times)
#' @param sample_times the vector of times at which samples could be taken
#' @param nsamps the number of samples per individual
#' @param antigenic_map_long the collapsed antigenic map for long term cross reactivity, after multiplying by sigma1
#' @param antigenic_map_short the collapsed antigenic map for short term cross reactivity, after multiplying by sigma2
#' @param repeats number of repeated samples for each year
#' @param mus default NULL, optional vector of boosting parameters for each strain
#' @param mu_indices default NULL, optional vector giving the index of `mus` that each strain uses the boosting parameter from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @param measurement_bias default NULL, optional vector of measurement bias shift parameters for each strain
#' @param measurement_indices default NULL, optional vector giving the index of `measurement_bias` that each strain uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @param add_noise defaults TRUE, adds observation error to simulated titres
#' @return a data frame with columns individual, samples, virus and titre of simulated data
#' @export
#' @seealso \code{\link{simulate_individual}}
simulate_group <- function(n_indiv,
                           theta,
                           infection_histories,
                           strain_isolation_times,
                           sample_times,
                           nsamps,
                           antigenic_map_long,
                           antigenic_map_short,
                           repeats = 1,
                           mus = NULL,
                           mu_indices = NULL,
                           measurement_bias = NULL,
                           measurement_indices = NULL,
                           add_noise = TRUE) {
  dat <- NULL
  ## For each individual
  for (i in 1:n_indiv) {
    ## Choose random sampling times
    ## If there is one sampling time, then repeat the same sampling time
    if (length(sample_times) == 1) {
      samps <- rep(sample_times, nsamps)
    } else {
      samps <- sample(sample_times, nsamps)
      samps <- samps[order(samps)]
    }

    virus_samples <- rep(strain_isolation_times, length(samps))
    data_indices <- rep(length(strain_isolation_times), length(samps))
    virus_indices <- match(virus_samples, strain_isolation_times) - 1
    y <- as.data.frame(simulate_individual(
      theta, infection_histories[i, ],
      samps, data_indices, virus_samples,
      virus_indices,
      antigenic_map_long,
      antigenic_map_short,
      strain_isolation_times,
      repeats, mus, mu_indices,
      measurement_bias, measurement_indices,
      add_noise
    ))
    ## Record individual ID
    y$indiv <- i
    colnames(y) <- c("samples", "virus", "titre", "individual")
    ## Combine data
    dat <- rbind(dat, y[, c("individual", "samples", "virus", "titre")])
  }
  return(dat)
}

#' Simulate individual data
#'
#' Simulates a full set of titre data for an individual with known theta and infection_history.
#' @inheritParams simulate_group
#' @param infection_history the vector of 1s and 0s giving presence/absence of infections
#' @param data_indices see \code{\link{create_posterior_func}}
#' @param virus_indices see \code{\link{create_posterior_func}}
#' @return a data frame with columns samples, virus and titre of simulated data
#' @export
#' @seealso \code{\link{infection_model_indiv}}
simulate_individual <- function(theta,
                                infection_history,
                                sampling_times,
                                data_indices,
                                strain_isolation_times,
                                virus_indices,
                                antigenic_map_long,
                                antigenic_map_short,
                                strains,
                                repeats = 1,
                                mus = NULL,
                                mu_indices = NULL,
                                measurement_bias = NULL,
                                measurement_indices = NULL,
                                add_noise = TRUE) {
    numberStrains <- length(strains)
    dat <- matrix(ncol = 3, nrow = length(strain_isolation_times) * repeats)

    additional_arguments <- NULL
    if (!is.null(mus)) additional_args <- list("mus"=mus,"boosting_vec_indices"=mu_indices-1)

    titres <- titre_data_individual(
        theta, infection_history, strains, seq_along(strains) - 1, sampling_times,
        data_indices, match(strain_isolation_times, strains) - 1,
        antigenic_map_long, antigenic_map_short, numberStrains, 0,
        additional_arguments
    )
    titres <- rep(titres, each = repeats)
    sampling_times <- rep(sampling_times, data_indices)
    sampling_times <- rep(sampling_times, each = repeats)
    dat[, 1] <- sampling_times
    dat[, 2] <- rep(strain_isolation_times, each = repeats)
    if (add_noise) {
        if (!is.null(measurement_indices)) {
            dat[, 3] <- add_noise(titres, theta, measurement_bias, measurement_indices[match(dat[, 2], strains)])
        } else {
            dat[, 3] <- add_noise(titres, theta, NULL, NULL)
        }
    } else {
        dat[, 3] <- titres
    }
    return(dat)
}

#' Add noise
#'
#' Adds truncated noise to titre data
#' @param y the titre
#' @param theta a vector with MAX_TITRE and error parameters
#' @return a noisy titre
#' @export
add_noise <- function(y, theta, measurement_bias = NULL, indices = NULL) {
  ## Draw from normal
  if (!is.null(measurement_bias)) {
    noise_y <- floor(rnorm(length(y), mean = y + measurement_bias[indices], sd = theta["error"]))
  } else {
    noise_y <- floor(rnorm(length(y), mean = y, sd = theta["error"]))
  }

  ## If outside of bounds, truncate
  noise_y[noise_y < 0] <- 0
  noise_y[noise_y > theta["MAX_TITRE"]] <- theta["MAX_TITRE"]
  return(noise_y)
}

#' Simulate attack rates
#'
#' Given a number of possible infection years, simulates attack rates from a log normal distribution with specified mean and standard deviation.
#' @param infection_years the number of infection years
#' @param mean_par the mean of the log normal
#' @param sd_par the sd of the log normal
#' @param large_first_year simulate an extra large attach rate in the first year?
#' @param big_year_mean if large first year, what mean to use?
#' @return a vector of attack rates
#' @export
simulate_attack_rates <- function(infection_years, mean_par = 0.15, sd_par = 0.5,
                                  large_first_year = FALSE, big_year_mean = 0.5) {
  attack_year <- rlnorm(infection_years, meanlog = log(mean_par) - sd_par^2 / 2, sdlog = sd_par)
  if (large_first_year) attack_year[1] <- rlnorm(1, meanlog = log(big_year_mean) - (sd_par / 2)^2 / 2, sdlog = sd_par / 2)
  return(attack_year)
}

#' Simulate infection histories
#'
#' Given a vector of infection probabilities and potential infection times, simulates infections for each element of ages (ie. each element is an individual age. Only adds infections for alive individuals)
#' @param p_inf a vector of attack rates (infection probabilities) for each year
#' @param strain_isolation_times the vector of possible infection times
#' @param sampling_times vector of potential sampling times
#' @param ages a vector of ages for each individual
#' @return a list with a matrix of infection histories for each individual in ages and the true attack rate for each epoch
#' @export
simulate_infection_histories <- function(p_inf, strain_isolation_times, sampling_times, ages) {
  n_strains <- length(p_inf) # How many strains
  n_indiv <- length(ages) # How many individuals
  indivs <- 1:n_indiv
  ## Empty matrix
  infection_histories <- matrix(0, ncol = n_strains, nrow = n_indiv)

  ## Simulate attack rates
  attack_rates <- p_inf

  ## Should this be necessary?
  attack_rates[attack_rates > 1] <- 1
  ARs <- numeric(n_strains)
  ## For each strain (ie. each infection year)
  for (i in 1:n_strains) {
    # If there are strains circulating beyond the max sampling times, then alive==0
    if (max(sampling_times) >= strain_isolation_times[i]) {
      ## Find who was alive (all we need sampling_times for is its max value)
      alive <- (max(sampling_times) - ages) <= strain_isolation_times[i]
    } else {
      alive <- rep(0, n_indiv)
    }

    ## Sample a number of infections for the alive individuals, and set these entries to 1
    y <- round(length(indivs[alive]) * attack_rates[i])
    # y <- rbinom(1, length(indivs[alive]),attack_rates[i])
    ARs[i] <- y / length(indivs[alive])
    x <- sample(indivs[alive], y)
    infection_histories[x, i] <- 1
  }
  return(list(infection_histories, ARs))
}

#' Generates attack rates from an SIR model with fixed beta/gamma, specified final attack rate and the number of time "buckets" to solve over ie. buckets=12 returns attack rates for 12 time periods
#' @export
generate_ar_annual <- function(AR, buckets) {
  SIR_odes <- function(t, x, params) {
    S <- x[1]
    I <- x[2]
    R <- x[3]
    inc <- x[4]

    beta <- params[1]
    gamma <- params[2]
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    dinc <- beta * S * I
    list(c(dS, dI, dR, dinc))
  }
  R0 <- 1.2
  gamma <- 1 / 5
  beta <- R0 * gamma
  t <- seq(0, 360, by = 0.1)
  results <- as.data.frame(deSolve::ode(
    y = c(S = 1, I = 0.0001, R = 0, inc = 0),
    times = t, func = SIR_odes,
    parms = c(beta, gamma)
  ))
  incidence <- diff(results$inc)
  incidence <- incidence * AR / sum(incidence)
  group <- 360 * 10 / buckets
  monthly_risk <- colSums(matrix(incidence, nrow = group))
  return(monthly_risk)
}

#' @export
simulate_ars_buckets <- function(infection_years, buckets, mean_par = 0.15, sd_par = 0.5,
                                 large_first_year = FALSE, big_year_mean = 0.5) {
  n <- ceiling(length(infection_years) / buckets)
  attack_year <- rlnorm(n, meanlog = log(mean_par) - sd_par^2 / 2, sdlog = sd_par)
  if (large_first_year) attack_year[1] <- rlnorm(1, meanlog = log(big_year_mean) - (sd_par / 2)^2 / 2, sdlog = sd_par / 2)
  ars <- NULL

  for (i in seq_along(attack_year)) {
    ars <- c(ars, generate_ar_annual(attack_year[i], buckets))
  }

  ars <- ars[1:length(infection_years)]
  return(ars)
}

#' @export
simulate_ars_spline <- function(infection_years, buckets, mean_par = 0.15, sd_par = 0.5, large_first_year = FALSE, big_year_mean = 0.5, knots, theta) {
  infection_years <- infection_years[seq(1, length(infection_years), by = buckets)] / buckets
  n <- length(infection_years)
  attack_year <- rlnorm(n, meanlog = log(mean_par) - sd_par^2 / 2, sdlog = sd_par)
  if (large_first_year) attack_year[1] <- rlnorm(1, meanlog = log(big_year_mean) - (sd_par / 2)^2 / 2, sdlog = sd_par / 2)
  ars <- generate_lambdas(attack_year, knots, theta, n, buckets)
  return(ars)
}



#' Simulate full data set
#'
#' Simulates a full data set for a given set of parameters etc.
#' @param par_tab the full parameter table controlling parameter ranges and values
#' @param group which group index to give this simulated data
#' @param n_indiv number of individuals to simulate
#' @param strain_isolation_times vector of strain ciruclation times
#' @param sampling_times possible sampling times for the individuals
#' @param nsamps the number of samples each individual has
#' @param antigenic_map the raw antigenic map with colnames x_coord, y_coord and inf_years
#' @param titre_sensoring what proportion of titres are randomly missing?
#' @param age_min minimum age to simulate
#' @param age_max maximum age to simulate
#' @param attack_rates a vector of attack_rates to be used in the simulation
#' @param repeats number of repeats for each year
#' @param mu_indices default NULL, optional vector giving the index of `mus` that each strain uses the boosting parameter from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @return a list with: 1) the data frame of titre data as returned by \code{\link{simulate_group}}; 2) a matrix of infection histories as returned by \code{\link{simulate_infection_histories}}; 3) a vector of ages
#' @export
simulate_data <- function(par_tab, group = 1, n_indiv, buckets = 12,
                          strain_isolation_times, sampling_times, nsamps = 2,
                          antigenic_map,
                          titre_sensoring = 0,
                          age_min = 5, age_max = 80,
                          attack_rates,
                          repeats = 1,
                          mu_indices = NULL,
                          measurement_indices = NULL,
                          add_noise = TRUE) {

  ## Check attack_rates entry
    check_attack_rates(attack_rates, strain_isolation_times)

    ## Extract parameter type indices from par_tab, to split up
    ## similar parameters in model solving functions
    option_indices <- which(par_tab$type == 0)
    theta_indices <- which(par_tab$type %in% c(0, 1))
    measurement_indices_par_tab <- which(par_tab$type == 3)
    mu_indices_par_tab <- which(par_tab$type == 6)
    
    ## Extract parameters
    par_names_theta <- par_tab[theta_indices, "names"]
    pars <- par_tab$values
    theta <- pars[theta_indices]
    names(theta) <- par_names_theta
    
    mus <- NULL
    if (!is.null(mu_indices)) {
        mus <- pars[mu_indices_par_tab]
    }
    
    measurement_bias <- NULL
    if (!is.null(measurement_indices)) {
        measurement_bias <- pars[measurement_indices_par_tab]
    }

  ## Check the inputs of par_tab
  check_par_tab(par_tab)

  ## Create antigenic map for short and long term boosting
  antigenic_map1 <- melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")])

  antigenic_map_long <- 1 - theta["sigma1"] * c(antigenic_map1)
  antigenic_map_short <- 1 - theta["sigma2"] * c(antigenic_map1)

  antigenic_map_long[antigenic_map_long < 0] <- 0
  antigenic_map_short[antigenic_map_short < 0] <- 0

  ## Simulate ages
  ages <- floor(runif(n_indiv, age_min, age_max))

  ## Simulate infection histories
  tmp <- simulate_infection_histories(
    attack_rates, strain_isolation_times,
    sampling_times, ages
  )

  infection_history <- tmp[[1]]
  ARs <- tmp[[2]]

  ## Simulate titre data
  y <- simulate_group(
    n_indiv, theta, infection_history, strain_isolation_times, sampling_times,
    nsamps, antigenic_map_long, antigenic_map_short, repeats,
    mus, mu_indices, measurement_bias, measurement_indices, add_noise
  )

  ## Randomly censor titre values
  y$titre <- y$titre * sample(c(NA, 1), nrow(y), prob = c(titre_sensoring, 1 - titre_sensoring), replace = TRUE)
  y$run <- 1
  y$group <- 1

  DOB <- max(sampling_times) - ages
  ages <- data.frame("individual" = 1:n_indiv, "DOB" = DOB)
  attack_rates <- data.frame("year" = strain_isolation_times, "AR" = ARs)
  return(list(data = y, infection_histories = infection_history, ages = ages, attack_rates = attack_rates, lambdas = attack_rates))
}

#' Create useable antigenic map
#'
#' Creates an antigenic map from an input data frame that can be used to calculate cross reactivity. This will end up being an NxN matrix, where there are N strains circulating.
#' @param anti.map.in can either be a 1D antigenic line to calculate distance from, or a two dimensional matrix with x and y coordinates on an antigenic map
#' @return the euclidean antigenic distance between each pair of viruses in anti.map.in
#' @export
melt_antigenic_coords <- function(anti.map.in) { # anti.map.in can be vector or matrix - rows give inf_years, columns give location
  # Calculate antigenic distances
  if (is.null(dim(anti.map.in))) { # check if input map is one or 2 dimensions
    # If 1D antigenic 'line' defined, calculate distances directory from input
    (dmatrix <- sapply(anti.map.in, function(x) {
      y <- abs(anti.map.in - x)
      y
    }))
  } else { # If 2D antigenic map defined, calculate distances directory from input
    (dmatrix <- apply(
      anti.map.in, 1,
      function(x) {
        y <- sqrt(colSums(apply(anti.map.in, 1, function(y) {
          (y - x)^2
        })))
        y
      }
    ))
  }
}

#' Generate antigenic map
#'
#' Firts a smoothing spline through a set of antigenic coordinates, and uses this to predict antigenic coordinates for all potential infection time points
#' @param antigenic_distances a data frame of antigenic coordinates, with columns labelled X, Y and Strain for x coord, y coord and Strain label respectively. Strain labels should be in the virus_key vector which is at the top of this source code.
#' @param buckets the number of epochs per year. 1 means that each year has 1 strain; 12 means that each year has 12 strains (monthly resolution)
#' @param spar to be passed to smooth.spline
#' @return a fitted antigenic map
#' @seealso \code{\link{generate_antigenic_map_flexible}}
#' @export
generate_antigenic_map <- function(antigenic_distances, buckets = 1, spar = 0.3) {
  ## Following assumptions:
  ## 1. X31 == 1969
  ## 2. PE2009 is like the strain circulating in 2010
  virus_key <- c(
    "HK68" = 1968, "EN72" = 1972, "VI75" = 1975, "TX77" = 1977, "BK79" = 1979, "SI87" = 1987, "BE89" = 1989, "BJ89" = 1989,
    "BE92" = 1992, "WU95" = 1995, "SY97" = 1997, "FU02" = 2002, "CA04" = 2004, "WI05" = 2005, "PE06" = 2006
  ) * buckets
  antigenic_distances$Strain <- virus_key[antigenic_distances$Strain]
  fit <- smooth.spline(antigenic_distances$X, antigenic_distances$Y, spar = spar)
  x_line <- lm(data = antigenic_distances, X ~ Strain)
  Strain <- seq(1968 * buckets, 2016 * buckets - 1, by = 1)
  x_predict <- predict(x_line, data.frame(Strain))
  y_predict <- predict(fit, x = x_predict)
  fit_dat <- data.frame(x = y_predict$x, y = y_predict$y)
  fit_dat$strain <- Strain
  colnames(fit_dat) <- c("x_coord", "y_coord", "inf_years")
  return(fit_dat)
}
#' Generate antigenic map, flexible
#'
#' Firts a smoothing spline through a set of antigenic coordinates, and uses this to predict antigenic coordinates for all potential infection time points. This version is more flexible than \code{\link{generate_antigenic_map}}, and allows the user to specify "clusters" to assume that strains circulating in a given period are all identical, rather than on a continuous path through space as a function of time.
#' @param antigenic_distances a data frame of antigenic coordinates, with columns labelled X, Y and Strain for x coord, y coord and Strain label respectively. "Strain" should be a single number giving the year of circulation of that strain.
#' @param buckets = 1 the number of epochs per year. 1 means that each year has 1 strain; 12 means that each year has 12 strains (monthly resolution)
#' @param clusters = NULL a data frame of cluster labels, indicating which cluster each circulation year belongs to. Note that each row (year) gets repeated by the number of buckets. Column names "year" and "cluster_used"
#' @param use_clusters = FALSE if TRUE, uses the clusters data frame above, otherwise just returns as normal
#' @param spar = 0.3 to be passed to smooth.spline
#' @param year_min = 1968 first year in the antigenic map (usually 1968)
#' @param year_max = 2016 last year in the antigenic map
#' @return a fitted antigenic map
#' @seealso \code{\link{generate_antigenic_map}}
#' @export
generate_antigenic_map_flexible <- function(antigenic_distances, buckets = 1, clusters = NULL,
                                            use_clusters = FALSE, spar = 0.3, year_min = 1968, year_max = 2016) {
  ## Convert strains to correct time dimensions
  antigenic_distances$Strain <- antigenic_distances$Strain * buckets
  ## Fit spline through antigenic coordinates
  fit <- smooth.spline(antigenic_distances$X, antigenic_distances$Y, spar = spar)

  ## Work out relationship between strain circulation time and x coordinate
  x_line <- lm(data = antigenic_distances, X ~ Strain)

  ## Enumerate all strains that could circulate
  Strain <- seq(year_min * buckets, year_max * buckets - 1, by = 1)

  ## Predict x and y coordinates for each possible strain from this spline
  x_predict <- predict(x_line, data.frame(Strain))
  y_predict <- predict(fit, x = x_predict)

  fit_dat <- data.frame(x = y_predict$x, y = y_predict$y)
  fit_dat$strain <- Strain
  colnames(fit_dat) <- c("x_coord", "y_coord", "inf_years")

  ## If using clusters
  if (use_clusters) {
    ## Repeat each row by the number of buckets per year
    clusters <- clusters[rep(seq_len(nrow(clusters)), each = buckets), ]

    ## Enumerate out such that each row has a unique time
    clusters$year <- seq(year_min * buckets, length.out = nrow(clusters))

    ## Which time point does each cluster start?
    cluster_starts <- ddply(clusters, ~cluster_used, function(x) x$year[1])
    colnames(cluster_starts)[2] <- "first_cluster_year"
    clusters1 <- merge(clusters, cluster_starts, by = c("cluster_used"))
    clusters1 <- clusters1[, c("cluster_used", "first_cluster_year", "year")]
    colnames(fit_dat)[3] <- "first_cluster_year"

    ## Merge on "inf_years" with first_cluster_year, such that all viruses
    ## in a cluster have the same location as the first virus in that cluster
    fit_dat <- fit_dat[fit_dat$first_cluster_year %in% clusters1$first_cluster_year, ]
    fit_dat <- merge(clusters1, fit_dat, by = "first_cluster_year")
    fit_dat <- fit_dat[, c("x_coord", "y_coord", "year")]
    colnames(fit_dat)[3] <- "inf_years"
    fit_dat <- fit_dat[order(fit_dat$inf_years), ]
  }

  return(fit_dat)
}

#' @export
euc_distance <- function(i1, i2, fit_dat) {
  return(sqrt((fit_dat[i1, "x_coord"] - fit_dat[i2, "x_coord"])^2 + (fit_dat[i1, "y_coord"] - fit_dat[i2, "y_coord"])^2))
}
