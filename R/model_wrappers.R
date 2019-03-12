#' Solves the antibody kinetics model for a single individual
#'
#' @param par_tab the parameter table controlling information such as bounds, initial values etc
#' @param infection_history vector of 1s and 0s giving the infection history
#' @param sample_time the time at which a blood sample is taken
#' @param antigenic_map a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param measurement_indices_by_time if not NULL, then use these indices to specify which measurement bias parameter index corresponds to which time
#' @param mu_indices if not NULL, then use these indices to specify which boosting parameter index corresponds to which time
#' @param DOB the date of birth of this individual, matching the time resolution of the model
#' @return a vector of predicted titres
#' @export
solve_model_individual <- function(par_tab, infection_history, sample_time, antigenic_map,
                                   measurement_indices_by_time = NULL, mu_indices = NULL, DOB = 0) {
  ## Unique strains that an individual could be exposed to
  strain_isolation_times <- unique(antigenic_map$inf_years)
  ## The entry of each strain in the antigenic map table
  strain_indices <- match(strain_isolation_times, strain_isolation_times) - 1
  antigenic_map_melted <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))

  ## Extract parameter type indices from par_tab, to split up
  ## similar parameters in model solving functions
  theta_indices <- which(par_tab$type %in% c(0, 1))
  measurement_indices_par_tab <- which(par_tab$type == 3)
  mu_indices_par_tab <- which(par_tab$type == 6)
  #########################################################

  par_names_theta <- par_tab[theta_indices, "names"]
  pars <- par_tab$values

  ## Find which options are being used in advance for speed
  use_measurement_bias <- (length(measurement_indices_par_tab) > 0) & !is.null(measurement_indices_by_time)
  titre_shifts <- NULL
  expected_indices <- NULL
  measurement_bias <- NULL
  use_strain_dependent <- (length(mu_indices) > 0) & !is.null(mu_indices)
  additional_arguments <- NULL


  if (use_measurement_bias) {
    expected_indices <- measurement_indices_by_time[match(titre_dat$virus, strain_isolation_times)]
  }

  if (use_strain_dependent) {
    additional_arguments <- list(
      "boosting_vec_indices" = mu_indices - 1,
      "mus" = rep(2, length(strain_isolation_times))
    )
  }

  theta <- pars[theta_indices]
  mus <- pars[mu_indices_par_tab]

  if (use_measurement_bias) {
    measurement_bias <- pars[measurement_indices_par_tab]
  }

  if (use_strain_dependent) {
    additional_arguments[["mus"]] <- mus
  }
  names(theta) <- par_names_theta

  ## Work out short and long term boosting cross reactivity - C++ function
  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

  y <- infection_model_indiv(
    theta, infection_history, strain_isolation_times, strain_indices, sample_time,
    strain_indices, antigenic_map_long, antigenic_map_short, length(infection_history),
    DOB, additional_arguments
  )

  if (use_measurement_bias) {
    measurement_bias <- pars[measurement_indices_par_tab]
    titre_shifts <- measurement_bias[expected_indices]
    y <- y + titre_shifts
  }
  y
}
