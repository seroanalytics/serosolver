#' Get number alive
#'
#' Given the antibody_data data frame, calculates the number that are alive (alive to be infected, that is) for each time in times
#' @param antibody_data the data frame of titre data. See \code{\link{example_antibody_data}}
#' @param times the vector of times to calculate number alive for
#' @return a vector giving the number alive in each time point
#' @family get_summary
#' @examples
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' times <- unique(example_antigenic_map$inf_times)
#' get_n_alive(example_antibody_data, times)
#' @export
get_n_alive <- function(antibody_data, times) {
  DOBs <- unique(antibody_data[, c("individual", "birth")])[, 2]
  age_mask <- create_age_mask(DOBs, times)
  sample_mask <- create_sample_mask(antibody_data, times)
  masks <- data.frame(cbind(age_mask, sample_mask))
  n_alive <- sapply(seq(1, length(times)), function(x)
    nrow(masks[masks$age_mask <= x & masks$sample_mask >= x, ]))
}

#' Get DOBs
#'
#' Gets the dates of birth for each individual in antibody_data
#' @param antibody_data the data frame of antibody data. See \code{\link{example_antibody_data}}
#' @return a data frame with individual and birth
#' @family get_summary
#' @examples
#' data(example_antibody_data)
#' get_DOBs(antibody_data)
#' @export
get_DOBs <- function(antibody_data){
    DOBs <- unique(antibody_data[,c("individual","birth")])
}

#' Get number alive by location
#'
#' Given the antibody_data data frame with entries for location, calculates the number that are alive (alive to be infected, that is) for each time in times by location
#' @param antibody_data the data frame of antibody data. See \code{\link{example_antibody_data}}
#' @param times the vector of times to calculate number alive for
#' @param melt_data if TRUE, returns a melted data frame. Returns a wide matrix otherwise.
#' @return a matrix giving the number alive in each time point in each location
#' @family get_summary
#' @examples
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' times <- unique(example_antigenic_map$inf_times)
#' get_n_alive_group(example_antibody_data, times)
#' @export
get_n_alive_group <- function(antibody_data, times, melt_data = FALSE) {
  DOBs <- unique(antibody_data[, c("individual", "population_group", "birth")])
  age_mask <- create_age_mask(DOBs[, "birth"], times)
  sample_mask <- create_sample_mask(antibody_data, times)
  masks <- data.frame(cbind(age_mask, sample_mask))
  DOBs <- cbind(DOBs, masks)
  n_alive <- plyr::ddply(DOBs, ~population_group, function(y) sapply(seq(1, length(times)), function(x)
      nrow(y[y$age_mask <= x & y$sample_mask >= x, ])))
  n_alive <- as.matrix(n_alive[, 2:ncol(n_alive)])
  colnames(n_alive) <- times
  if (melt_data) {
    n_alive <- data.frame(n_alive)
    n_alive$population_group <- 1:nrow(n_alive)
    n_alive <- reshape2::melt(n_alive, id.vars = c("population_group"))
    colnames(n_alive)[2] <- "j"
    n_alive$j <- times[n_alive$j]
    colnames(n_alive)[3] <- "n_alive"
  }
  n_alive
}

#' @export
create_prior_lookup <- function(antibody_data, possible_exposure_times, infection_model_prior_shape1, beta1, n_alive=NULL){
    if(is.null(n_alive)){
        n_alive <- get_n_alive(antibody_data, possible_exposure_times)
    }
    lookup_tab <- matrix(nrow=max(n_alive)+1,ncol=length(possible_exposure_times))
    max_alive <- max(n_alive)
    for(i in seq_along(possible_exposure_times)){
        results <- rep(-100000, max_alive+1)
        m <- seq_len(n_alive[i]+1)-1
        results[1:(n_alive[i]+1)] <- lbeta(infection_model_prior_shape1 + m, n_alive[i] - m + beta1) + lbeta(infection_model_prior_shape1, beta1)
        lookup_tab[,i] <- results        
    }
    lookup_tab[!is.finite(lookup_tab)] <- -100000
    lookup_tab[is.nan(lookup_tab)] <- -100000
    lookup_tab
}

#' @export
create_prior_lookup_groups <- function(antibody_data, possible_exposure_times, infection_model_prior_shape1, beta1, n_alive=NULL){
    if(is.null(n_alive)){
        n_alive <- get_n_alive_group(antibody_data, possible_exposure_times)
    }
    lookup_tab <- array(-100000, dim=c(max(n_alive)+1,length(possible_exposure_times),nrow(n_alive)))
    max_alive <- max(n_alive)
    for(g in 1:nrow(n_alive)){
        for(i in seq_along(possible_exposure_times)){
            results <- rep(-100000, max_alive+1)
            m <- seq_len(n_alive[g,i]+1)-1
            results[1:(n_alive[g,i]+1)] <- lbeta(infection_model_prior_shape1 + m, n_alive[g,i] - m + beta1) + lbeta(infection_model_prior_shape1, beta1)
            lookup_tab[,i,g] <- results        
        }
    }
    
    lookup_tab[!is.finite(lookup_tab)] <- -100000
    lookup_tab[is.nan(lookup_tab)] <- -100000
    lookup_tab
}


#' Create age mask
#'
#' Converts a data frame of individual ages to give the index of the first infection that individual could have had
#' @param DOBs the vector of dates of birth, same time units as possible_exposure_times
#' @param possible_exposure_times the vector of times that individuals can be infected
#' @return a vector giving the first index of possible_exposure_times that an individual can be infected
#' @family create_masks
#' @examples
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' DOBs <- unique(example_antibody_data[,c("individual","birth")])
#' age_mask <- create_age_mask(DOBs$birth, times)
#' @export
create_age_mask <- function(DOBs, possible_exposure_times) {
  age_mask <- sapply(DOBs, function(x) {
    if (is.na(x)) {
      1
    } else {
      which(as.numeric(x <= possible_exposure_times) > 0)[1]
    }
  })
  return(age_mask)
}
#' Create sample mask
#'
#' Converts a data frame of individual sampling times to give the index of the last infection that individual could have had (as we can't observe anything after the last observation)
#' @param antibody_data the data frame of antibody data. See \code{\link{example_antibody_data}}
#' @param possible_exposure_times the vector of times that individuals can be infected
#' @return a vector giving the last index of possible_exposure_times that an individual can be infected
#' @family create_masks
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' sample_mask <- create_sample_mask(example_antibody_data, times)
#' @export
create_sample_mask <- function(antibody_data, possible_exposure_times) {
  ids <- unique(antibody_data$individual)
  sample_mask <- sapply(ids, function(x) {
    sample_times <- antibody_data$sample_time[antibody_data$individual == x]
    max(which(max(sample_times) >= possible_exposure_times))
  })
  return(sample_mask)
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
  full_inf_chain <- data.table::CJ(i = 1:max(inf_chain$i), j = 1:max(inf_chain$j), samp_no = min(inf_chain$samp_no):max(inf_chain$samp_no))
  inf_chain <- data.table::data.table(apply(inf_chain, 2, as.numeric))
  summary_with_non_infections <- merge(inf_chain, full_inf_chain, by = c("samp_no", "j", "i"), all = TRUE)
  summary_with_non_infections[is.na(summary_with_non_infections$x), "x"] <- 0
  colnames(summary_with_non_infections) <- c("samp_no", "j", "individual", "x")
  expanded_chain <- data.table::dcast(summary_with_non_infections, samp_no + individual ~ j, value.var = "x")
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
  tmp_names <- colnames(chain)
  tmp_names <- tmp_names[!(tmp_names) %in% c("samp_no","chain_no","likelihood","posterior_prob","prior_prob")]
  best_pars <- as.numeric(as.data.frame(chain)[which.max(as.data.frame(chain)[, "posterior_prob"]), tmp_names])
  names(best_pars) <- tmp_names
  return(best_pars)
}

#' Index pars
#'
#' Given an MCMC chain, returns the parameters at the specified index
#' @param chain the MCMC chain
#' @param samp_no the sample number
#' @param index if not using `samp_no`, `index returns the desired row number`
#' @param chain_no the chain number to subset
#' @return a named vector of the best parameters
#' @family mcmc_diagnostics
#' @examples
#' \dontrun{
#' mcmc_chains <- load_theta_chains()
#' pars <- get_index_pars(mcmc_chains$chain, index=100)
#' }
#' @export
get_index_pars <- function(chain, samp_no=NULL,index=NULL,chain_no=NULL) {
  tmp_names <- colnames(chain)
  tmp_names <- tmp_names[!(tmp_names) %in% c("samp_no","chain_no","likelihood","posterior_prob","prior_prob")]

  if(is.null(samp_no) & is.null(index)){
    stop("Error in get_index_pars - must specify one of samp_no or index")
  }
  
  if(!is.null(index)){
    pars <- as.numeric(chain[index, tmp_names])
  } else if(!is.null(chain_no) & "chain_no" %in% colnames(chain)) {
    pars <- as.numeric(chain[chain$samp_no == samp_no & chain$chain_no==chain_no, tmp_names])
  } else {
    pars <- as.numeric(chain[chain$samp_no == samp_no, tmp_names])
  }
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
#' Wrapper function to protect calls to a function. If the function does not compute correctly, returns -100000.
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

#' Protect function (posterior function)
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
#' @param f the function to be protected
#' @return the protected function
#' @export
#' @useDynLib serosolver
protect_posterior <- function(f) {
  function(...) {
    tryCatch(f(...), error = function(e) {
      message("Error thrown in posterior solving function: ", e$message)
      message("If this error seems fairly cryptic (e.g., subscript out of bounds), the error likely comes from the C++ code, indicating an error with the antibody model function. You may have misspecified some parameters in par_tab, or there may be a misalignment between par_tab, antibody_data and antigenic_map.")
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

#' Describe infection history priors
#' @export
describe_priors <- function() {
  message("Which version to use in serosolver? The following text describes the proposal step for updating infection histories.")
  message("Version 1: Beta prior on per time attack rates. Explicit FOI on each epoch using probability of infection term. Proposal performs N `flip` proposals at random locations in an individual's infection history, switching 1->0 or 0->1. Otherwise, swaps the contents of two random locations")
  message("Version 2: Beta prior on per time attack rates. Gibbs sampling of infection histories as in Indian Buffet Process papers, integrating out each probability of infection term.")
  message("Version 3: Beta prior on probability of infection for an individual, assuming independence between individuals. Samples from a beta binomial with shape1 and shape2 specified by the par_tab input. Proposes nInfs moves at a time for add/remove, or when swapping, swaps locations up to moveSize time steps away")
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


#' Pad par_tab with Beta distribution shape parameters
#'
#' Pads par_tab with a new row for each infection epoch, such that each epoch has its own shape1 and shape2
#' @param par_tab as per usual
#' @param n_times the number of additional rows to add for each alpha and beta
#' @examples
#' n_times <- 40
#' data(example_par_tab)
#' new_par_tab <- pad_infection_model_prior_parameters(example_par_tab, n_times)
#' @export
pad_infection_model_prior_parameters <- function(par_tab, n_times) {
  shape1_row <- par_tab[par_tab$names == "infection_model_prior_shape1", ]
  shape2_row <- par_tab[par_tab$names == "infection_model_prior_shape2", ]

  for (i in 1:(n_times - 1)) {
    par_tab <- rbind(par_tab, shape1_row, shape2_row)
  }
  par_tab
}

## From prodlim package - finds matching rows between two data frames. "Thus the function returns a vector with the row numbers of (first) matches of its first argument in its second.", https://www.rdocumentation.org/packages/prodlim/versions/2018.04.18/topics/row.match
#' @export
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

#' Setup antibody data indices
#'
#' Sets up a large list of pre-indexing and pre-processing to speed up the model solving during MCMC fitting.
#' Note that this should be `antibody_data` after subsetting to only `run==1`, as we will figure out elsewhere which solves to use as repeats
#' @inheritParams create_posterior_func
#' @param verbose if TRUE, brings warning messages
#' @return a very long list. See source code directly.
#' @seealso \code{\link{create_posterior_func}}
#' @export
setup_antibody_data_for_posterior_func <- function(antibody_data, antigenic_map=NULL, possible_exposure_times=NULL,
                                              age_mask = NULL, n_alive = NULL,verbose=FALSE) {
  essential_colnames <- c("individual", "sample_time", "measurement", "biomarker_id", "biomarker_group","population_group")
  ## How many observation types are there?
  n_indiv <- length(unique(antibody_data$individual))
  unique_biomarker_groups <- unique(antibody_data$biomarker_group)
  n_biomarker_groups <- length(unique_biomarker_groups)
  
  ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
  if (!is.null(antigenic_map)) {
      possible_exposure_times_tmp <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
      if(!is.null(possible_exposure_times) & !identical(possible_exposure_times, possible_exposure_times_tmp)){
          if(verbose) message(cat("Warning: provided possible_exposure_times argument does not match entries in the antigenic map. Please make sure that there is an entry in the antigenic map for each possible circulation time. Using the antigenic map times.\n"))
      }
      possible_exposure_times <- possible_exposure_times_tmp
      
      ## If no observation types assumed, set all to 1.
      if (!("biomarker_group" %in% colnames(antigenic_map))) {
        if(verbose) message(cat("Note: no biomarker_group detection in antigenic_map. Aligning antigenic map with par_tab.\n"))
          antigenic_map_tmp <- replicate(n_biomarker_groups,antigenic_map,simplify=FALSE)
          for(biomarker_group in unique_biomarker_groups){
              antigenic_map_tmp[[biomarker_group]]$biomarker_group <- biomarker_group
          }
          antigenic_map <- do.call(rbind,antigenic_map_tmp)
      }
      
  } else {
      ## Create a dummy map with entries for each observation type
      antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,
                                  "inf_times"=rep(possible_exposure_times, n_biomarker_groups), 
                                  "biomarker_group"=rep(unique_biomarker_groups,each=length(possible_exposure_times)))
  }
  possible_exposure_times <- unique(antigenic_map$inf_times)
  
  ## Create a melted antigenic map for each observation type
  antigenic_maps_melted <- lapply(unique_biomarker_groups, function(b){
      tmp <- antigenic_map[antigenic_map$biomarker_group == b,]
      c(melt_antigenic_coords(tmp[,c("x_coord","y_coord")]))
    })
  
  biomarker_id_indices <- match(antibody_data$biomarker_id, possible_exposure_times) - 1 ## For each biomarker_id tested, what is its index in the antigenic map?
  exposure_id_indices <- match(possible_exposure_times, possible_exposure_times) - 1 ## For each biomarker_id that circulated, what is its index in the antigenic map?

  ## Get unique measurement sets for each individual at
  ## each sampling time, for each observation type, for each repeat
  ## ie. each row of this is a unique blood sample and observation type taken
  sample_data <- unique(antibody_data[, c("individual", "sample_time", "biomarker_group")])
  sample_data_start <- cumsum(c(0,plyr::ddply(sample_data, .(individual, biomarker_group), nrow)$V1))
  sample_times <- sample_data$sample_time ## What were the times that these samples were taken?
  
  type_data <- unique(antibody_data[,c("individual","biomarker_group")])
  type_data_start <- cumsum(c(0, plyr::ddply(type_data, .(individual), nrow)$V1))
  biomarker_groups <- type_data$biomarker_group
  
  ## How many rows in the antibody data correspond to each unique individual, sample, observation type?
  ## ie. each element of this vector corresponds to one set of antibodies that need to be predicted
  nrows_per_sample <- plyr::ddply(antibody_data, .(individual,biomarker_group, sample_time), nrow)$V1
  antibody_data_start <- cumsum(c(0,nrows_per_sample))

    ## Get unique groups
  groups <- unique(antibody_data$population_group)
  group_table <- unique(antibody_data[, c("individual", "population_group")])
  group_id_vec <- group_table$population_group - 1

  if (!is.null(antibody_data$birth)) {
    DOBs <- unique(antibody_data[, c("individual", "birth")])[, 2]
  } else {
    DOBs <- rep(min(possible_exposure_times), n_indiv)
    antibody_data$DOB <- min(possible_exposure_times)
  }
  age_mask <- create_age_mask(DOBs, possible_exposure_times)
  sample_mask <- create_sample_mask(antibody_data, possible_exposure_times)
  masks <- data.frame(cbind(age_mask, sample_mask))

  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(antibody_data, possible_exposure_times)
  }

  return(list(
    "biomarker_groups"=biomarker_groups,
    "antigenic_map_melted" = antigenic_maps_melted,
    "possible_exposure_times" = possible_exposure_times,
    "exposure_id_indices" = exposure_id_indices,
    
    "sample_times" = sample_times,
    
    "antibody_data_start"=antibody_data_start,
    "nrows_per_sample"=nrows_per_sample,
    "sample_data_start"=sample_data_start,
    "type_data_start"=type_data_start,
    
    "group_id_vec" = group_id_vec,
    
    ## This one I need to figure out -- does it need to have a different set of indices for each observation type? Probably not
    "biomarker_id_indices" = biomarker_id_indices,
    "n_indiv" = n_indiv,
    "age_mask" = age_mask,
    "sample_mask" = sample_mask,
    "n_alive" = n_alive,
    "DOBs" = DOBs
  ))
}


#' @export
euc_distance <- function(i1, i2, fit_data) {
  return(sqrt((fit_data[i1, "x_coord"] - fit_data[i2, "x_coord"])^2 + (fit_data[i1, "y_coord"] - fit_data[i2, "y_coord"])^2))
}


#' Create useable antigenic map
#'
#' Creates an antigenic map from an input data frame that can be used to calculate cross reactivity. This will end up being an NxN matrix, where there are N strains circulating.
#' @param anti.map.in can either be a 1D antigenic line to calculate distance from, or a two dimensional matrix with x and y coordinates on an antigenic map
#' @return the euclidean antigenic distance between each pair of viruses in anti.map.in
#' @export
melt_antigenic_coords <- function(anti.map.in) { # anti.map.in can be vector or matrix - rows give inf_times, columns give location
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

  fit_data <- data.frame(x = y_predict$x, y = y_predict$y)
  fit_data$strain <- Strain
  colnames(fit_data) <- c("x_coord", "y_coord", "inf_times")

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
    colnames(fit_data)[3] <- "first_cluster_year"

    ## Merge on "inf_times" with first_cluster_year, such that all viruses
    ## in a cluster have the same location as the first virus in that cluster
    fit_data <- fit_data[fit_data$first_cluster_year %in% clusters1$first_cluster_year, ]
    fit_data <- merge(clusters1, fit_data, by = "first_cluster_year")
    fit_data <- fit_data[, c("x_coord", "y_coord", "year")]
    colnames(fit_data)[3] <- "inf_times"
    fit_data <- fit_data[order(fit_data$inf_times), ]
  }

  return(fit_data)
}


#' Pad infection history chain
#'
#' Given that the MCMC sampler only stores present infections (ie. there are no entries for 0s from the infection history matrix), for some summaries we need to add these 0s back in to avoid bias.
#' @param inf_chain the data table with infection history samples from \code{\link{serosolver}}
#' @return the same inf_chain that was passed in, but with 0s for missing i/j/samp_no combinations
#' @export
pad_inf_chain <- function(inf_chain) {
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (is.null(inf_chain$population_group)) {
    inf_chain$population_group <- 1
  }
  
  is <- unique(inf_chain$i)
  js <- unique(inf_chain$j)
  
  samp_nos <- unique(inf_chain$samp_no)
  chain_nos <- unique(inf_chain$chain_no)
  groups <- unique(inf_chain$population_group)
  expanded_values <- data.table::CJ(
    i = is,
    j = js,
    samp_no = samp_nos,
    chain_no = chain_nos,
    population_group = groups
  )
  diff_infs <- fsetdiff(expanded_values, inf_chain[, c("i", "j", "samp_no", "chain_no", "population_group")])
  diff_infs$x <- 0
  inf_chain <- rbind(inf_chain, diff_infs)
  return(inf_chain)
}

#' Unregister parallel backgrounds forcefully
#'
#' When the serosolver function uses a parallel backend but is not closed correctly, this can confuse `dopar` if the function is called again. This function tidies the parallel backend and corrects the error.
#' @return NULL
#' @export
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#' Generate starting antibody levels
#'
#' Generates either random or data-driven starting antibody levels for each measured biomarker group/id combination per individual. This is mostly used elsewhere in the serosolver model
#' @param antibody_data the antibody data, see \code{\link{example_antibody_data}}
#' @param start_level_summary string telling the function how to use the `antibody_data` object to create starting values. One of: min, max, mean, median, full_random.
#' @param randomize if TRUE and data is discretized, then sets the starting level to a random value between floor(x) and floor(x)+1
#' @return a list with two objects: 1) a tibble giving the starting antibody level for each individual, biomarker group and biomarker_id combinations; 2) a list of indices (starting at 0) of length matching `nrow(antibody_data)` giving the index of the antibody starting level to use for each measurement
#' @family antigenic_maps
#' @examples
#' \dontrun{
#' create_start_level_data(example_antibody_data,"min",FALSE)
#' create_start_level_data(example_antibody_data,"min",TRUE)
#' create_start_level_data(example_antibody_data,"max",FALSE)
#' create_start_level_data(example_antibody_data,"max",TRUE)
#' create_start_level_data(example_antibody_data,"mean",FALSE)
#' create_start_level_data(example_antibody_data,"mean",TRUE)
#' create_start_level_data(example_antibody_data,"median",FALSE)
#' create_start_level_data(example_antibody_data,"median",TRUE)
#' create_start_level_data(example_antibody_data,"other",FALSE)
#' create_start_level_data(example_antibody_data,"other",TRUE)
#' create_start_level_data(example_antibody_data,"full_random",FALSE)
#' create_start_level_data(example_antibody_data,"full_random",TRUE)
#' }
#' @export
create_start_level_data <- function(antibody_data, start_level_summary = "min", randomize=FALSE){
  ## Get earliest measurement per biomarker ID as starting level. Need to decide if using min, max, mean or median.
  starting_levels <- antibody_data %>% 
    dplyr::select(individual,biomarker_id,biomarker_group, measurement,sample_time) %>% 
    dplyr::group_by(individual, biomarker_id, biomarker_group) %>% 
    dplyr::filter(sample_time == min(sample_time)) %>% 
    dplyr::group_by(individual,biomarker_id, biomarker_group)
  
  ## Bounds of data and whether it's continuous or discrete
  data_controls <- antibody_data %>% 
    dplyr::mutate(diff_from_self = measurement - floor(measurement)) %>% 
    dplyr::group_by(biomarker_group) %>% 
    dplyr::summarize(min_measurement=min(measurement,na.rm=TRUE),max_measurement=max(measurement,na.rm=TRUE),type=if_else(any(diff_from_self != 0),"continuous","discrete")) 
  
  if(start_level_summary == "min"){
    starting_levels <- starting_levels %>% dplyr::summarize(starting_level = min(measurement,na.rm=TRUE)) 
  }else if(start_level_summary == "max"){
    starting_levels <- starting_levels %>% dplyr::summarize(starting_level = max(measurement,na.rm=TRUE)) 
  }else if(start_level_summary == "mean"){
    starting_levels <- starting_levels %>% dplyr::summarize(starting_level = mean(measurement,na.rm=TRUE)) 
  }else if(start_level_summary == "median"){
    starting_levels <- starting_levels %>% dplyr::summarize(starting_level = median(measurement,na.rm=TRUE)) 
  } else if(start_level_summary == "full_random"){
    starting_levels <- starting_levels %>% dplyr::select(individual, biomarker_group, biomarker_id) %>% dplyr::distinct() %>% dplyr::left_join(data_controls,by="biomarker_group") %>% dplyr::ungroup() %>% dplyr::mutate(starting_level = runif(n(), min_measurement,max_measurement))
  } else {
    starting_levels <- starting_levels %>% dplyr::select(individual, biomarker_group, biomarker_id) %>% dplyr::distinct() %>% dplyr::mutate(starting_level = 0)
  }
  starting_levels <- starting_levels %>% dplyr::select(individual, biomarker_id, biomarker_group, starting_level)
  starting_levels <- starting_levels %>% dplyr::arrange(individual, biomarker_group, biomarker_id) %>% dplyr::ungroup() %>% dplyr::mutate(start_index = 1:n())
  if(randomize){
    starting_levels <- starting_levels %>% dplyr::left_join(data_controls,by="biomarker_group") %>% ungroup() %>% dplyr::mutate(starting_level = if_else(type=="discrete",runif(n(),floor(starting_level), floor(starting_level) + 1),starting_level))
  }
  
  start_indices <- left_join(antibody_data, starting_levels, by=c("individual", "biomarker_id", "biomarker_group"))
  start_indices
}



