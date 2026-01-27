#' Load a starting parameter table from file
#'
#' Searches the specified working directory for a file matching "start_tab.csv" and reads in the first matching file as a data frame. This is typically used to initialize MCMC parameter values.
#' @param location A character string specifying the path to the directory where the function should search for the file. Defaults to the current working directory.
#' @return A data frame containing the contents of the first matching start tab file found, or `NULL` if none is found.
#' @family load_data_functions
#' @seealso [load_antibody_data()], [load_mcmc_chains()]
#' @examples
#' \dontrun{
#'   # Load the start tab from the current working directory
#'   start_tab <- load_start_tab()
#'
#'   # Load from a specific subdirectory
#'   start_tab <- load_start_tab("mcmc_chains/")
#' }
#' @export
load_start_tab <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_start_tab.csv"))
  files_old <- Sys.glob(file.path(location, "*_startTab.csv"))
  files <- c(files, files_old)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  start_tab <- read.csv(files[1], stringsAsFactors = FALSE)
  return(start_tab)
}

#' Load serological antibody data from CSV file
#'
#' Searches the specified working directory for a csv file with "antibody_data.csv" or "titer_dat.csv" in the name and reads in the first file found as a data frame.
#' @inheritParams load_start_tab
#' @return A data frame with antibody data, or `NULL` if no matching file is found.
#' @family load_data_functions
#' @seealso [load_start_tab()], [load_mcmc_chains()]
#' @examples
#' \dontrun{
#'   antibody_data <- load_antibody_data()
#'   antibody_data <- load_antibody_data("data/")
#' }
#'
#' @export
load_antibody_data <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_antibody_data.csv"))
  files_antibody_data <- Sys.glob(file.path(location, "*_titre_dat.csv"))
  files_titreDat <- Sys.glob(file.path(location, "*_titreDat.csv"))
  files <- c(files, files_antibody_data, files_titreDat)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  antibody_data <- read.csv(files[1], stringsAsFactors = FALSE)
  return(antibody_data)
}


#' Load in antigenic map file from CSV
#'
#' Reads in an antigenic map from the first CSV file found matching the pattern `*_antigenic_map.csv`. 
#' @inheritParams load_start_tab
#' @return A data frame containing the antigenic map. Returns `NULL` if no matching file is found.
#' @family load_data_functions
#' @seealso [load_antibody_data()], [load_start_tab()]
#' @examples
#' \dontrun{
#'   antigenic_map <- load_antigenic_map_file()
#'   antigenic_map <- load_antigenic_map_file("data/")
#' }
#' @export
load_antigenic_map <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_antigenic_map.csv"))
  files_old <- Sys.glob(file.path(location, "*_antigenicMap.csv"))
  files <- c(files, files_old)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  antigenic_map <- read.csv(files[1], stringsAsFactors = FALSE)
  return(antigenic_map)
}


#' Load MCMC chains from CSV files for antibody kinetics parameters (theta)
#'
#' Searches the given working directory for MCMC outputs from \code{\link{serosolver}} matching "_chain.csv", reads these in, subsets for burn in and thinning, and formats as both lists and a combined data frame.
#' @param location Character string path to the directory to search. Defaults to `getwd()`.
#' @param par_tab Data frame; optional parameter table for filtering. Defaults to `NULL`.
#' @param estimated_only if TRUE, only returns free model parameters (par_tab$fixed == 0) if par_tab specified
#' @param thin Integer; thinning interval to apply to the MCMC samples. Defaults to `1` (no thinning).
#' @param burnin Integer; number of initial samples to discard. Defaults to `0`.
#' @param convert_mcmc Logical; whether to convert output to `coda::mcmc.list`. Defaults to `TRUE`.
#' @param verbose Logical; whether to print progress messages
#' @return a list with a) a list of each chain separately; b) a combined data frame, indexing each iteration by which chain it comes from
#' @family load_data_functions
#' @examples
#' \dontrun{load_theta_chains(location="mcmc_chains", par_tab=par_tab, estimated_only=TRUE,thin=10,burnin=5000,convert_mcmc=TRUE)}
#' @export
load_theta_chains <- function(location = getwd(), par_tab = NULL, estimated_only = TRUE, thin = 1, burnin = 0, convert_mcmc = TRUE,verbose=TRUE) {
  chains <- Sys.glob(file.path(location, "*_chain.csv"))
  if(verbose) message("Chains detected: ", length(chains), sep = "\t")
  if (length(chains) < 1) {
    if(verbose) message("Error - no chains found")
    return(NULL)
  }
  
  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains, read.csv)
  
  if(verbose) message("Highest MCMC sample interations: \n")
  if(verbose) lapply(read_chains, function(x) message(max(x$samp_no)))
  
  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[seq(1, nrow(x), by = thin), ])
  read_chains <- lapply(read_chains, function(x) x[x$samp_no > burnin, ])
  max_samp_no <- min(as.numeric(lapply(read_chains, function(x) max(x$samp_no))))
  read_chains <- lapply(read_chains, function(x) x[x$samp_no <= max_samp_no, ])
  unique_samp_nos <- lapply(read_chains, function(x) unique(x[, "samp_no"]))
  unique_samp_nos <- Reduce(intersect, unique_samp_nos)
  read_chains <- lapply(read_chains, function(x) x[x$samp_no %in% unique_samp_nos, ])
  
  for (i in 1:length(read_chains)) read_chains[[i]]$chain_no <- i
  
  ## Get the estimated parameters only
  if (estimated_only & !is.null(par_tab)) {
    fixed <- par_tab$fixed
    fixed_names <- par_tab$names[which(fixed == 0)]
    ## Go through the vector of strings called fixed_names and append a number to each non-unique name
    if (length(fixed_names) > 0) {
      fixed_names <- make.unique(fixed_names)
    }
    
    use_colnames <- intersect(c("samp_no", fixed_names, "posterior_prob", "likelihood", "prior_prob", "chain_no"), colnames(read_chains[[1]]))
    read_chains <- lapply(read_chains, function(x) x[, use_colnames])
  }
  
  ## Try to create an MCMC list. This might not work, which is why we have a try catch
  list_chains <- read_chains
  if (convert_mcmc) {
    list_chains <- tryCatch({
      tmp_list <- lapply(read_chains, coda::as.mcmc)
      tmp_list <- coda::as.mcmc.list(tmp_list)
    },
    warning = function(w) {
      if(verbose) message(w)
      NULL
    }, error = function(e) {
      if(verbose) message(e)
      NULL
    },
    finally = {
      tmp_list
    }
    )
  }
  
  chain <- do.call("rbind", read_chains)
  if (convert_mcmc) chain <- coda::as.mcmc(chain)
  return(list("list" = list_chains, "chain" = chain))
}

#' Load MCMC chains from CSV files for infection histories
#'
#' Searches the given working directory for MCMC outputs from \code{\link{serosolver}} ending "infection_histories.csv", loads these in, subsets for burn in and thinning, and formats as both lists and a combined data table.
#' @inheritParams load_theta_chains
#' @param chain_subset if not NULL, a vector of indices to only load and store a subset of the chains detected. eg. chain_subset = 1:3 means that only the first 3 detected files will be processed.
#' @return a list with a) a list of each chain as a data table separately; b) a combined data table, indexing each iteration by which chain it comes from
#' @family load_data_functions
#' @seealso [load_mcmc_chains()], [load_start_tab()], [plot_infection_histories()]
#' @examples
#' \dontrun{load_infection_chains(thin=10,burnin=5000,chain_subset=1:3)}
#' @export
load_infection_chains <- function(location = getwd(), thin = 1, burnin = 0, chain_subset = NULL, verbose=TRUE) {
  chains <- Sys.glob(file.path(location, "*_infection_histories.csv"))
  chains_old <- Sys.glob(file.path(location, "*_infectionHistories.csv"))
  chains <- c(chains, chains_old)
  if (!is.null(chain_subset)) chains <- chains[chain_subset]
  if(verbose) message("Chains detected: ", chains, sep = "\n")
  if (length(chains) < 1) {
    if(verbose) message("Error - no chains found")
    return(NULL)
  }
  
  if(verbose) message("Reading in infection history chains. May take a while.")
  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains, data.table::fread)
  
  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[samp_no > burnin, ])
  read_chains <- lapply(read_chains, function(x) {
    samp_nos <- unique(x$samp_no)
    samp_nos <- samp_nos[seq(1, length(samp_nos), by = thin)]
    x <- x[samp_no %in% samp_nos, ]
  })
  max_samp_no <- min(as.numeric(lapply(read_chains, function(x) max(x$samp_no))))
  read_chains <- lapply(read_chains, function(x) x[samp_no <= max_samp_no, ])
  
  if(verbose) message("Number of rows: ")
  if(verbose) lapply(read_chains, function(x) message(nrow(x)))
  
  for (i in 1:length(read_chains)) read_chains[[i]]$chain_no <- i
  chain <- do.call("rbind", read_chains)
  return(list("list" = read_chains, "chain" = chain))
}


#' Load MCMC chains for the antibody kinetics parameters and infection histories from CSV outputs
#'
#' Reads in all MCMC chains for theta and infection histories from the specified directory matcing files ending "_chain.csv", adding in the total number of infections
#' @inheritParams load_theta_chains
#' @return a list of the concatenated and individual chains (4 elements, either data frames of coda::mcmc objects)
#' @family load_data_functions
#' @examples
#' \dontrun{load_mcmc_chains(par_tab=par_tab, estimated_only=TRUE,thin=10,burnin=5000,convert_mcmc=TRUE)}
#' @export
load_mcmc_chains <- function(location = getwd(), par_tab = NULL, estimated_only = FALSE, thin = 1, burnin = 0, convert_mcmc = FALSE, verbose=TRUE) {
  ## Load in theta chains
  theta_chains <- load_theta_chains(location, par_tab, estimated_only, thin, burnin, verbose=verbose)
  ## Load in infection history chains
  inf_chains <- load_infection_chains(location, thin, burnin, verbose=verbose)
  
  chain <- theta_chains$chain
  inf_chain <- inf_chains$chain
  theta_list_chains <- theta_chains$list
  inf_list_chains <- inf_chains$list
  
  ## Concatenate total number of infections per MCMC sample
  total_inf_chain <- get_total_number_infections(inf_chain, FALSE)
  chain <- merge(chain, total_inf_chain)
  
  ## Combine total number of infections with theta chain
  list_total_inf_chains <- lapply(inf_list_chains, function(x) get_total_number_infections(x, FALSE))
  theta_list_chains <- Map(
    function(x, y) merge(data.table(x), data.table(y), by = c("samp_no", "chain_no")),
    theta_list_chains, list_total_inf_chains
  )
  theta_list_chains <- lapply(theta_list_chains, function(x) x[order(x[, "samp_no"]), ])
  unique_samp_nos <- lapply(theta_list_chains, function(x) unique(x[, "samp_no"])$samp_no)
  unique_samp_nos <- Reduce(intersect, unique_samp_nos)
  theta_list_chains <- lapply(theta_list_chains, function(x) x[x$samp_no %in% unique_samp_nos, ])
  
  ## Convert to MCMC objects if desired
  if (convert_mcmc) {
    theta_list_chains <- lapply(theta_list_chains, coda::as.mcmc)
    chain <- as.mcmc(chain)
  }
  
  return(list("theta_chain" = chain, "inf_chain" = inf_chain, "theta_list_chains" = theta_list_chains, "inf_list_chains" = inf_list_chains))
}