#' Read in start_tab
#'
#' Searches the working directory for a file called "start_tab.csv" and reads in the first one found
#' @param location string giving the relative file path to be searched
#' @return data frame, the first found starting parameter tables in the current working directory
#' @family load_data_functions
#' @examples
#' \dontrun{load_start_tab()}
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

#' Read in titre_dat
#'
#' Searches the working directory for a file with "titre_dat.csv" in the name and reads in the first one found
#' @inheritParams load_start_tab
#' @return data frame, the first found starting parameter tables in the current working directory
#' @family load_data_functions
#' @examples
#' \dontrun{load_titre_dat()}
#' @export
load_titre_dat <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_titre_dat.csv"))
  files_old <- Sys.glob(file.path(location, "*_titreDat.csv"))
  files <- c(files, files_old)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  titre_dat <- read.csv(files[1], stringsAsFactors = FALSE)
  return(titre_dat)
}


#' Read in antigenic_map
#'
#' Searches the working directory for a file with "antigenic_map.csv" in the name and reads in the first one found
#' @inheritParams load_start_tab
#' @return data frame, the first found antigenic map tables in the current working directory
#' @family load_data_functions
#' @examples
#' \dontrun{load_antigenic_map_file()}
#' @export
load_antigenic_map_file <- function(location = getwd()) {
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




#' Combine theta and infection history chains
#'
#' Reads in all MCMC chains for theta and infection histories from the specified directory, adding in the total number of infections
#' @inheritParams load_theta_chains
#' @return a list of the concatenated and individual chains (4 elements, either data frames of coda::mcmc objects)
#' @family load_data_functions
#' @examples
#' \dontrun{load_mcmc_chains(par_tab=par_tab, unfixed=TRUE,thin=10,burnin=5000,convert_mcmc=TRUE)}
#' @export
load_mcmc_chains <- function(location = getwd(), par_tab = NULL, unfixed = TRUE, thin = 1, burnin = 0, convert_mcmc = FALSE) {
    ## Load in theta chains
    theta_chains <- load_theta_chains(location, par_tab, unfixed, thin, burnin)
    ## Load in infection history chains
  inf_chains <- load_infection_chains(location, thin, burnin)

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
    function(x, y) merge(data.table(x), data.table(y), by = c("sampno", "chain_no")),
    theta_list_chains, list_total_inf_chains
  )
  theta_list_chains <- lapply(theta_list_chains, function(x) x[order(x[, "sampno"]), ])
  unique_sampnos <- lapply(theta_list_chains, function(x) unique(x[, "sampno"])$sampno)
  unique_sampnos <- Reduce(intersect, unique_sampnos)
  theta_list_chains <- lapply(theta_list_chains, function(x) x[x$sampno %in% unique_sampnos, ])

    ## Convert to MCMC objects if desired
  if (convert_mcmc) {
    theta_list_chains <- lapply(theta_list_chains, as.mcmc)
    chain <- as.mcmc(chain)
  }

  return(list("theta_chain" = chain, "inf_chain" = inf_chain, "theta_list_chains" = theta_list_chains, "inf_list_chains" = inf_list_chains))
}

#' Load MCMC chains for theta
#'
#' Searches the given working directory for MCMC outputs from \code{\link{run_MCMC}}, loads these in, subsets for burn in and thinning, and formats as both lists and a combined data frame.
#' @param location defaults to current working directory. Gives relative file path to look for files ending in "_chain.csv"
#' @param par_tab if not NULL, can use this to only extract free model parameters
#' @param unfixed if TRUE, only returns free model parameters (par_tab$fixed == 0) if par_tab specified
#' @param thin thin the chains by every thin'th sample
#' @param burnin discard the first burnin samples from the MCMC chain
#' @param convert_mcmc if TRUE, converts everything to MCMC objects (from the `coda` R package)
#' @return a list with a) a list of each chain separately; b) a combined data frame, indexing each iteration by which chain it comes from
#' @family load_data_functions
#' @examples
#' \dontrun{load_theta_chains(par_tab=par_tab, unfixed=TRUE,thin=10,burnin=5000,convert_mcmc=TRUE)}
#' @export
load_theta_chains <- function(location = getwd(), par_tab = NULL, unfixed = TRUE, thin = 1, burnin = 0, convert_mcmc = TRUE) {
  chains <- Sys.glob(file.path(location, "*_chain.csv"))
  message(cat("Chains detected: ", length(chains), sep = "\t"))
  if (length(chains) < 1) {
    message("Error - no chains found")
    return(NULL)
  }

  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains, read.csv)

  message(cat("Highest MCMC sample interation: "))
  lapply(read_chains, function(x) message(max(x$sampno)))
  
  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[seq(1, nrow(x), by = thin), ])
  read_chains <- lapply(read_chains, function(x) x[x$sampno > burnin, ])
  max_sampno <- min(as.numeric(lapply(read_chains, function(x) max(x$sampno))))
  read_chains <- lapply(read_chains, function(x) x[x$sampno <= max_sampno, ])
  unique_sampnos <- lapply(read_chains, function(x) unique(x[, "sampno"]))
  unique_sampnos <- Reduce(intersect, unique_sampnos)
  read_chains <- lapply(read_chains, function(x) x[x$sampno %in% unique_sampnos, ])

  message("Number of rows: ")
  message(lapply(read_chains, nrow))

  for (i in 1:length(read_chains)) read_chains[[i]]$chain_no <- i

  ## Get the estimated parameters only
  if (unfixed & !is.null(par_tab)) {
    fixed <- par_tab$fixed
    use_colnames <- intersect(c("sampno", par_tab$names[which(fixed == 0)], "lnlike", "likelihood", "prior_prob", "chain_no"), colnames(read_chains[[1]]))
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
      message(w)
      NULL
    }, error = function(e) {
      message(e)
      NULL
    },
    finally = {
      tmp_list
    }
    )
  }

  chain <- do.call("rbind", read_chains)
  if (convert_mcmc) chain <- as.mcmc(chain)
  return(list("list" = list_chains, "chain" = chain))
}

#' Load MCMC chains for infection histories
#'
#' Searches the given working directory for MCMC outputs from \code{\link{run_MCMC}}, loads these in, subsets for burn in and thinning, and formats as both lists and a combined data table.
#' @param location defaults to current working directory. Where to look for MCMC chains? These are files ending in "_infection_histories.csv"
#' @inheritParams load_theta_chains
#' @param chain_subset if not NULL, a vector of indices to only load and store a subset of the chains detected. eg. chain_subset = 1:3 means that only the first 3 detected files will be processed.
#' @return a list with a) a list of each chain as a data table separately; b) a combined data table, indexing each iteration by which chain it comes from
#' @family load_data_functions
#' @examples
#' \dontrun{load_infection_chains(thin=10,burnin=5000,chain_subset=1:3)}
#' @export
load_infection_chains <- function(location = getwd(), thin = 1, burnin = 0, chain_subset = NULL) {
  chains <- Sys.glob(file.path(location, "*_infection_histories.csv"))
  chains_old <- Sys.glob(file.path(location, "*_infectionHistories.csv"))
  chains <- c(chains, chains_old)
  if (!is.null(chain_subset)) chains <- chains[chain_subset]
  message(cat("Chains detected: ", chains, sep = "\n"))
  if (length(chains) < 1) {
    message("Error - no chains found")
    return(NULL)
  }

  message("Reading in infection history chains. May take a while.")
  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains, data.table::fread)

  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[sampno > burnin, ])
  read_chains <- lapply(read_chains, function(x) {
    sampnos <- unique(x$sampno)
    sampnos <- sampnos[seq(1, length(sampnos), by = thin)]
    x <- x[sampno %in% sampnos, ]
  })
  max_sampno <- min(as.numeric(lapply(read_chains, function(x) max(x$sampno))))
  read_chains <- lapply(read_chains, function(x) x[sampno <= max_sampno, ])

  message("Number of rows: ")
  print(lapply(read_chains, nrow))

  for (i in 1:length(read_chains)) read_chains[[i]]$chain_no <- i
  chain <- do.call("rbind", read_chains)
  return(list("list" = read_chains, "chain" = chain))
}


#' Get total number of infections
#'
#' Finds the total number of infections for each iteration of an MCMC chain
#' @inheritParams plot_infection_history_chains_time
#' @return a data table
#' @examples
#' \dontrun{
#' inf_chain <- load_infection_chains(thin=10,burnin=5000,chain_subset=1:3)
#' n_infs <- get_total_number_infections(inf_chain$chain, pad_chain=FALSE)
#' }
#' @export
get_total_number_infections <- function(inf_chain, pad_chain = TRUE) {
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  n_strain <- max(inf_chain$j)
  data.table::setkey(inf_chain, "sampno", "chain_no")
  ## For each individual, how many infections did they have in each sample in total?
  n_inf_chain <- inf_chain[, list(total_infections = sum(x)), by = key(inf_chain)]
}

#' Estimate vector mode
#'
#' @param x the vector to be estimated
#' @return the estimated mode of the given vector of values
#' @examples
#' x <- runif(1000)
#' y <- estimate_mode(x)
#' @export
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
