#' Get number alive
#'
#' Given the titre_dat data frame, calculates the number that are alive (alive to be infected, that is) for each time in times
#' @param titre_dat the data frame of titre data
#' @param times the vector of times to calculate number alive for
#' @return a vector giving the number alive in each time point
#' @export
get_n_alive <- function(titre_dat, times) {
  DOBs <- unique(titre_dat[, c("individual", "DOB")])[, 2]
  age_mask <- create_age_mask(DOBs, times)
  strain_mask <- create_strain_mask(titre_dat, times)
  masks <- data.frame(cbind(age_mask, strain_mask))
  n_alive <- sapply(seq(1, length(times)), function(x)
    nrow(masks[masks$age_mask <= x & masks$strain_mask >= x, ]))
}

#' Create age mask
#'
#' Converts a data frame of individual ages to give the index of the first infection that individual could have had
#' @param DOBs the vector of dates of birth, same time units as strain_isolation_times
#' @param strain_isolation_times the vector of times that individuals can be infected
#' @return a vector giving the first index of strain_isolation_times that an individual can be infected
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
#' @param titre_dat the data frame of titre data
#' @param strain_isolation_times the vector of times that individuals can be infected
#' @return a vector giving the last index of strain_isolation_times that an individual can be infected
#' @export
create_strain_mask <- function(titre_dat, strain_isolation_times) {
  ids <- unique(titre_dat$individual)
  strain_mask <- sapply(ids, function(x) {
    sample_times <- titre_dat$samples[titre_dat$individual == x]
    max(which(max(sample_times) >= strain_isolation_times))
  })
  return(strain_mask)
}


#' Expands default MCMC saved infChain
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
#' @export
get_index_pars <- function(chain, index) {
  tmp_names <- colnames(chain)[2:(ncol(chain) - 1)]
  pars <- as.numeric(chain[index, 2:(ncol(chain) - 1)])
  names(pars) <- tmp_names
  return(pars)
}

#' PDF - Rich's function to print to device without potential for bad errors
#'
#' Prints to pdf, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
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
  print("Which version to use in run_MCMC? The following text describes the proposal step for updating infection histories.")
  print("Version 1: Explicit FOI on each epoch using lambda term. Proposal performs N `flip` proposals at random locations in an individual's infectoin history, switching 1->0 or 0->1. Otherwise, swaps the contents of two random locations")
  print("Version 2: gibbs sampling of infection histories as in Indian Buffet Process papers")
  print("Version 3: samples from a beta binomial with alpha and beta specified by the par_tab input. Proposes nInfs moves at a time for add/remove, or when swapping, swaps locations up to moveSize time steps away")
  print("Version 4: gibbs sampling of infection histories using total number of infections across all times and all individuals as the prior")
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

  ## Firstly, how many rows in the titre data correspond to each unique individual, sample and titre repeat?
  ## ie. each element of this vector corresponds to one set of titres that need to be predicted
  nrows_per_blood_sample <- NULL
  for (i in 1:nrow(samples)) {
    nrows_per_blood_sample <- c(nrows_per_blood_sample, nrow(samples[titre_dat$individual == samples[i, "individual"] &
      titre_dat$samples == samples[i, "samples"] &
      titre_dat$run == samples[i, "run"], ]))
  }

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
  nrows_per_individual_in_data <- NULL
  for (individual in unique(individuals)) {
    nrows_per_individual_in_data <- c(nrows_per_individual_in_data, nrow(titre_dat[titre_dat$individual == individual, ]))
  }
  cum_nrows_per_individual_in_data <- cumsum(c(0, nrows_per_individual_in_data))

  if (!is.null(titre_dat$DOB)) {
    DOBs <- unique(titre_dat[, c("individual", "DOB")])[, 2]
  } else {
    DOBs <- rep(min(strain_isolation_times), n_indiv)
  }
  if (is.null(age_mask)) {
    if (!is.null(titre_dat$DOB)) {
      age_mask <- create_age_mask(DOBs, strain_isolation_times)
    } else {
      age_mask <- rep(1, n_indiv)
    }
  }
  strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
  masks <- data.frame(cbind(age_mask, strain_mask))
  if (is.null(n_alive)) {
    n_alive <- sapply(seq(1, length(strain_isolation_times)), function(x)
      nrow(masks[masks$age_mask <= x & masks$strain_mask >= x, ]))
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
    "nrows_per_blood_sample" = nrows_per_blood_sample,
    "measured_strain_indices" = measured_strain_indices,
    "n_indiv" = n_indiv,
    "age_mask" = age_mask,
    "strain_mask" = strain_mask,
    "n_alive" = n_alive,
    "DOBs" = DOBs
  ))
}
