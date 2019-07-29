#' Get number alive
#'
#' Given the titre_dat data frame, calculates the number that are alive (alive to be infected, that is) for each time in times
#' @param titre_dat the data frame of titre data. See \code{\link{example_titre_dat}}
#' @param times the vector of times to calculate number alive for
#' @return a vector giving the number alive in each time point
#' @family get_summary
#' @examples
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- unique(example_antigenic_map$inf_years)
#' get_n_alive(example_titre_dat, times)
#' @export
get_n_alive <- function(titre_dat, times) {
  DOBs <- unique(titre_dat[, c("individual", "DOB")])[, 2]
  age_mask <- create_age_mask(DOBs, times)
  strain_mask <- create_strain_mask(titre_dat, times)
  masks <- data.frame(cbind(age_mask, strain_mask))
  n_alive <- sapply(seq(1, length(times)), function(x)
    nrow(masks[masks$age_mask <= x & masks$strain_mask >= x, ]))
}

#' Get number alive bylocation
#'
#' Given the titre_dat data frame with entries for location, calculates the number that are alive (alive to be infected, that is) for each time in times by location
#' @param titre_dat the data frame of titre data. See \code{\link{example_titre_dat}}
#' @param times the vector of times to calculate number alive for
#' @param melt_dat if TRUE, returns a melted data frame. Returns a wide matrix otherwise.
#' @return a matrix giving the number alive in each time point in each location
#' @family get_summary
#' @examples
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- unique(example_antigenic_map$inf_years)
#' get_n_alive_group(example_titre_dat, times)
#' @export
get_n_alive_group <- function(titre_dat, times, melt_dat = FALSE) {
  DOBs <- unique(titre_dat[, c("individual", "group", "DOB")])
  age_mask <- create_age_mask(DOBs[, "DOB"], times)
  strain_mask <- create_strain_mask(titre_dat, times)
  masks <- data.frame(cbind(age_mask, strain_mask))
  DOBs <- cbind(DOBs, masks)
  n_alive <- plyr::ddply(DOBs, ~group, function(y) sapply(seq(1, length(times)), function(x)
      nrow(y[y$age_mask <= x & y$strain_mask >= x, ])))
  n_alive <- as.matrix(n_alive[, 2:ncol(n_alive)])
  colnames(n_alive) <- times
  if (melt_dat) {
    n_alive <- data.frame(n_alive)
    n_alive$group <- 1:nrow(n_alive)
    n_alive <- melt(n_alive, id.vars = c("group"))
    colnames(n_alive)[2] <- "j"
    n_alive$j <- as.numeric(as.factor(n_alive$j))
    colnames(n_alive)[3] <- "n_alive"
  }
  n_alive
}



#' Create age mask
#'
#' Converts a data frame of individual ages to give the index of the first infection that individual could have had
#' @param DOBs the vector of dates of birth, same time units as strain_isolation_times
#' @param strain_isolation_times the vector of times that individuals can be infected
#' @return a vector giving the first index of strain_isolation_times that an individual can be infected
#' @family create_masks
#' @examples
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_years
#' DOBs <- unique(example_titre_dat[,c("individual","DOB")])
#' age_mask <- create_age_mask(DOBs$DOB, times)
#' @export
create_age_mask <- function(DOBs, strain_isolation_times) {
  age_mask <- sapply(DOBs, function(x) {
    if (is.na(x)) {
      1
    } else {
      which(as.numeric(x <= strain_isolation_times) > 0)[1]
    }
  })
  return(age_mask)
}
#' Create strain mask
#'
#' Converts a data frame of individual sampling times to give the index of the last infection that individual could have had
#' @param titre_dat the data frame of titre data. See \code{\link{example_titre_dat}}
#' @param strain_isolation_times the vector of times that individuals can be infected
#' @return a vector giving the last index of strain_isolation_times that an individual can be infected
#' @family create_masks
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_years
#' strain_mask <- create_strain_mask(example_titre_dat, times)
#' @export
create_strain_mask <- function(titre_dat, strain_isolation_times) {
  ids <- unique(titre_dat$individual)
  strain_mask <- sapply(ids, function(x) {
    sample_times <- titre_dat$samples[titre_dat$individual == x]
    max(which(max(sample_times) >= strain_isolation_times))
  })
  return(strain_mask)
}


#' Expands default MCMC saved inf_chain
#'
#' The MCMC function saves sparse matrix summaries of the infection history chain to
#' save space and run time. This function returns the expanded infection history chain,
#' as in the original version of the code. Returned value is a data table with leftmost columns
#' giving sample number and individual, with columns expanded to the right for each infection
#' period.
#' @param inf_chain a data table with the MCMC saved infection history chain
#' @return the MCMC saved infection history expanded with infection times as columns
#' @export
expand_summary_inf_chain <- function(inf_chain) {
  full_inf_chain <- data.table::CJ(i = 1:max(inf_chain$i), j = 1:max(inf_chain$j), sampno = min(inf_chain$sampno):max(inf_chain$sampno))
  inf_chain <- data.table::data.table(apply(inf_chain, 2, as.numeric))
  summary_with_non_infections <- merge(inf_chain, full_inf_chain, by = c("sampno", "j", "i"), all = TRUE)
  summary_with_non_infections[is.na(summary_with_non_infections$x), "x"] <- 0
  colnames(summary_with_non_infections) <- c("sampno", "j", "individual", "x")
  expanded_chain <- data.table::dcast(summary_with_non_infections, sampno + individual ~ j, value.var = "x")
  return(expanded_chain)
}


#' Best pars
#'
#' Given an MCMC chain, returns the set of best fitting parameters (MLE)
#' @param chain the MCMC chain
#' @return a name vector of the best parameters
#' @family mcmc_diagnostics
#' @examples
#' \dontrun{
#' mcmc_chains <- load_theta_chains()
#' best_pars <- get_best_pars(mcmc_chains$chain)
#' }
#' @export
#' @useDynLib serosolver
get_best_pars <- function(chain) {
  tmp_names <- colnames(chain)[2:(ncol(chain) - 1)]
  best_pars <- as.numeric(chain[which.max(chain[, "lnlike"]), 2:(ncol(chain) - 1)])
  names(best_pars) <- tmp_names
  return(best_pars)
}

#' Index pars
#'
#' Given an MCMC chain, returns the parameters at the specified index
#' @param chain the MCMC chain
#' @param index the index
#' @return a named vector of the best parameters
#' @family mcmc_diagnostics
#' @examples
#' \dontrun{
#' mcmc_chains <- load_theta_chains()
#' pars <- get_index_pars(mcmc_chains$chain, 1000)
#' }
#' @export
get_index_pars <- function(chain, index) {
  tmp_names <- colnames(chain)[2:(ncol(chain) - 1)]
  pars <- as.numeric(chain[chain$sampno == index, 2:(ncol(chain) - 1)])
  names(pars) <- tmp_names
  return(pars)
}

#' PDF - Rich's function to print to device without potential for bad errors
#'
#' Prints to pdf, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @family safe_plot_saving
#' @export
to.pdf <- function(expr, filename, ..., verbose = TRUE) {
  if (verbose) {
    cat(sprintf("Creating %s\n", filename))
  }
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

#' PNG - Rich's function to print to device without potential for bad errors
#'
#' Prints to png, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @family safe_plot_saving
#' @export
to.png <- function(expr, filename, ..., verbose = TRUE) {
  if (verbose) {
    cat(sprintf("Creating %s\n", filename))
  }
  png(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

#' SVG - Rich's function to print to device without potential for bad errors
#'
#' Prints to SVG, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @family safe_plot_saving
#' @export
to.svg <- function(expr, filename, ..., verbose = TRUE) {
  if (verbose) {
    cat(sprintf("Creating %s\n", filename))
  }
  svg(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}


#' Protect function
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
#' @param f the function to be protected
#' @return the protected function
#' @export
#' @useDynLib serosolver
protect <- function(f) {
  function(...) {
    tryCatch(f(...), error = function(e) {
      message("caught error: ", e$message)
      -10000000
    })
  }
}

#' Convert to unit scale
toUnitScale <- function(x, min, max) {
  return((x - min) / (max - min))
}

#' Convert from unit scale to original scale
fromUnitScale <- function(x, min, max) {
  return(min + (max - min) * x)
}

#' @export
describe_proposals <- function() {
  message("Which version to use in run_MCMC? The following text describes the proposal step for updating infection histories.")
  message("Version 1: Beta prior on per time attack rates. Explicit FOI on each epoch using probability of infection term. Proposal performs N `flip` proposals at random locations in an individual's infection history, switching 1->0 or 0->1. Otherwise, swaps the contents of two random locations")
  message("Version 2: Beta prior on per time attack rates. Gibbs sampling of infection histories as in Indian Buffet Process papers, integrating out each probability of infection term.")
  message("Version 3: Beta prior on probability of infection for an individual, assuming independence between individuals. Samples from a beta binomial with alpha and beta specified by the par_tab input. Proposes nInfs moves at a time for add/remove, or when swapping, swaps locations up to moveSize time steps away")
  message("Version 4: Beta prior on probability of any infection. Gibbs sampling of infection histories using total number of infections across all times and all individuals as the prior")
}

#' @export
logistic_transform <- function(x, maxX) {
  return(maxX / (1 + exp(-x)))
}
#' @export
logit_transform <- function(p, maxX) {
  return(log(p / (maxX - p)))
}


#' Pad par_tab with alpha and betas
#'
#' Pads par_tab with a new row for each infection epoch, such that each epoch has its own alpha and beta
#' @param par_tab as per usual
#' @param n_times the number of additional rows to add for each alpha and beta
#' @examples
#' n_times <- 40
#' data(example_par_tab)
#' new_par_tab <- pad_alphas_and_betas(example_par_tab, n_times)
#' @export
pad_alphas_and_betas <- function(par_tab, n_times) {
  alpha_row <- par_tab[par_tab$names == "alpha", ]
  beta_row <- par_tab[par_tab$names == "beta", ]

  for (i in 1:(n_times - 1)) {
    par_tab <- rbind(par_tab, alpha_row, beta_row)
  }
  par_tab
}

## From prodlim package - finds matching rows between two data frames. "Thus the function returns a vector with the row numbers of (first) matches of its first argument in its second.", https://www.rdocumentation.org/packages/prodlim/versions/2018.04.18/topics/row.match
row.match <- function(x, table, nomatch = NA) {
  if (class(table) == "matrix") {
    table <- as.data.frame(table)
  }
  if (is.null(dim(x))) {
    x <- as.data.frame(matrix(x, nrow = 1))
  }
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}

#' Setup titre data indices
#'
#' Sets up a large list of pre-indexing and pre-processing to speed up the model solving during MCMC fitting.
#' @inheritParams create_posterior_func
#' @return a very long list. See source code directly.
#' @seealso \code{\link{create_posterior_func}}
#' @export
setup_titredat_for_posterior_func <- function(titre_dat, antigenic_map, age_mask = NULL, n_alive = NULL) {
  essential_colnames <- c("individual", "samples", "titre", "virus", "run", "group")

  strain_isolation_times <- antigenic_map$inf_years
  antigenic_map_melted <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))

  measured_strain_indices <- match(titre_dat$virus, antigenic_map$inf_years) - 1 ## For each virus tested, what is its index in the antigenic map?
  infection_strain_indices <- match(strain_isolation_times, strain_isolation_times) - 1 ## For each virus that circulated, what is its index in the antigenic map?

  ## Get unique measurement sets for each individual at
  ## each sampling time for each repeat
  ## ie. each row of this is a unique blood sample taken
  samples <- unique(titre_dat[, c("individual", "samples", "run")])
  sample_times <- samples$samples ## What were the times that these samples were taken?
  individuals <- samples$individual ## Who are the individuals that these samples correspond to?
  n_indiv <- length(unique(individuals))

  groups <- unique(titre_dat$group)
  group_table <- unique(titre_dat[, c("individual", "group")])
  group_id_vec <- group_table$group - 1

  ## Firstly, how many rows in the titre data correspond to each unique individual, sample and titre repeat?
  ## ie. each element of this vector corresponds to one set of titres that need to be predicted
  # nrows_per_blood_sample <- NULL
  nrows_per_blood_sample <- plyr::ddply(titre_dat, .(individual, samples, run), nrow)$V1

  ## Which indices in the sampling times vector correspond to each individual?
  ## ie. each contiguous pair of entries in this vector corresponds to the
  ## first and last entry in the samples matrix that correspond to each individual
  rows_per_indiv_in_samples <- c(0)
  for (individual in unique(individuals)) {
    rows_per_indiv_in_samples <- c(rows_per_indiv_in_samples, length(individuals[individuals == individual]))
  }
  rows_per_indiv_in_samples <- cumsum(rows_per_indiv_in_samples)

  ## Which indices in the titre data matrix correspond to each individual?
  ## And, how many rows match each individual?
  nrows_per_individual_in_data <- plyr::ddply(titre_dat, .(individual), nrow)$V1
  cum_nrows_per_individual_in_data <- cumsum(c(0, nrows_per_individual_in_data))

  if (!is.null(titre_dat$DOB)) {
    DOBs <- unique(titre_dat[, c("individual", "DOB")])[, 2]
  } else {
    DOBs <- rep(min(strain_isolation_times), n_indiv)
  }
  age_mask <- create_age_mask(DOBs, strain_isolation_times)
  strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
  masks <- data.frame(cbind(age_mask, strain_mask))

  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
  }

  return(list(
    "individuals" = individuals,
    "antigenic_map_melted" = antigenic_map_melted,
    "strain_isolation_times" = strain_isolation_times,
    "infection_strain_indices" = infection_strain_indices,
    "sample_times" = sample_times,
    "rows_per_indiv_in_samples" = rows_per_indiv_in_samples,
    "nrows_per_individual_in_data" = nrows_per_individual_in_data,
    "cum_nrows_per_individual_in_data" = cum_nrows_per_individual_in_data,
    "group_id_vec" = group_id_vec,
    "nrows_per_blood_sample" = nrows_per_blood_sample,
    "measured_strain_indices" = measured_strain_indices,
    "n_indiv" = n_indiv,
    "age_mask" = age_mask,
    "strain_mask" = strain_mask,
    "n_alive" = n_alive,
    "DOBs" = DOBs
  ))
}


#' @export
euc_distance <- function(i1, i2, fit_dat) {
  return(sqrt((fit_dat[i1, "x_coord"] - fit_dat[i2, "x_coord"])^2 + (fit_dat[i1, "y_coord"] - fit_dat[i2, "y_coord"])^2))
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
#' Fits a smoothing spline through a set of antigenic coordinates, and uses this to predict antigenic coordinates for all potential infection time points
#' @param antigenic_distances a data frame of antigenic coordinates, with columns labelled X, Y and Strain for x coord, y coord and Strain label respectively. Strain labels should be in the virus_key vector which is at the top of this source code. See \code{\link{example_antigenic_map}}
#' @param buckets the number of epochs per year. 1 means that each year has 1 strain; 12 means that each year has 12 strains (monthly resolution)
#' @param spar to be passed to smooth.spline
#' @return a fitted antigenic map
#' @family antigenic_maps
#' @examples
#' \dontrun{
#' antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
#' antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)
#' antigenic_map <- generate_antigenic_map(antigenic_coords, 1)
#' }
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
#' Fits a smoothing spline through a set of antigenic coordinates, and uses this to predict antigenic coordinates for all potential infection time points. This version is more flexible than \code{\link{generate_antigenic_map}}, and allows the user to specify "clusters" to assume that strains circulating in a given period are all identical, rather than on a continuous path through space as a function of time.
#' @param antigenic_distances a data frame of antigenic coordinates, with columns labelled X, Y and Strain for x coord, y coord and Strain label respectively. "Strain" should be a single number giving the year of circulation of that strain. See \code{\link{example_antigenic_map}}
#' @param buckets = 1 the number of epochs per year. 1 means that each year has 1 strain; 12 means that each year has 12 strains (monthly resolution)
#' @param clusters = NULL a data frame of cluster labels, indicating which cluster each circulation year belongs to. Note that each row (year) gets repeated by the number of buckets. Column names "year" and "cluster_used"
#' @param use_clusters = FALSE if TRUE, uses the clusters data frame above, otherwise just returns as normal
#' @param spar = 0.3 to be passed to smooth.spline
#' @param year_min = 1968 first year in the antigenic map (usually 1968)
#' @param year_max = 2016 last year in the antigenic map
#' @return a fitted antigenic map
#' @family antigenic_maps
#' @examples
#' \dontrun{
#' antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
#' antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)
#' antigenic_coords$Strain <- c(68,72,75,77,79,87,89,92,95,97,102,104,105,106) + 1900
#' antigenic_map <- generate_antigenic_map_flexible(antigenic_coords, buckets=1, year_min=1968, year_max=2015,spar=0.3)
#' 
#' times <- 1968:2010
#' n_times <- length(times)
#' clusters <- rep(1:5, each=10)
#' clusters <- clusters[1:n_times]
#' clusters <- data.frame(year=times, cluster_used=clusters)
#' antigenic_map <- generate_antigenic_map_flexible(antigenic_coords, buckets=1, 
#'                                                 clusters=clusters,use_clusters=TRUE,
#'                                                 year_min=1968, year_max=2010,spar=0.5)
#' }
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
    cluster_starts <- plyr::ddply(clusters, ~cluster_used, function(x) x$year[1])
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
