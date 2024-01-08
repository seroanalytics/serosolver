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

  message(cat("Highest MCMC sample interations: \n"))
  lapply(read_chains, function(x) message(max(x$sampno)))
  
  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[seq(1, nrow(x), by = thin), ])
  read_chains <- lapply(read_chains, function(x) x[x$sampno > burnin, ])
  max_sampno <- min(as.numeric(lapply(read_chains, function(x) max(x$sampno))))
  read_chains <- lapply(read_chains, function(x) x[x$sampno <= max_sampno, ])
  unique_sampnos <- lapply(read_chains, function(x) unique(x[, "sampno"]))
  unique_sampnos <- Reduce(intersect, unique_sampnos)
  read_chains <- lapply(read_chains, function(x) x[x$sampno %in% unique_sampnos, ])

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

#' Estimate attack rates
#'
#' Extracts attack rates from the MCMC output on infection histories
#' @param infection_histories the MCMC chain for infection histories
#' @param titre_dat the data frame of titre data
#' @param strain_isolation_times vector of the epochs of potential circulation
#' @param by_group if TRUE, stratifies by group ID
#' @param n_alive vector with the number of people alive in each year of circulation. Can be left as NULL, and ages will be used to infer this
#' @return a tibble with the estimated attack rates for each potential epoch of circulation
#' @export
get_attack_rates <- function(infection_histories, strain_isolation_times,titre_dat,n_alive=NULL, by_group=FALSE){
  if (!by_group) {
    titre_dat$group <- 1
    infection_histories$group <- 1
  }
  ## Find inferred total number of infections from the MCMC output
  ## Scale by number of individuals that were alive in each epoch
  ## and generate quantiles
  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
  }
  
  n_alive <- as.data.frame(n_alive)
  n_alive$group <- 1:nrow(n_alive)
  
  n_groups <- length(unique(titre_dat$group))
  n_alive_tot <- get_n_alive(titre_dat, strain_isolation_times)
  colnames(infection_histories)[1] <- "individual"
  infection_histories <- merge(infection_histories, data.table(unique(titre_dat[, c("individual", "group")])), by = c("individual","group"))
  years <- c(strain_isolation_times, max(strain_isolation_times) + 3)
  data.table::setkey(infection_histories, "sampno", "j", "chain_no", "group")
  tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
  tmp$taken <- years[tmp$j] %in% unique(titre_dat$samples)
  tmp$taken <- ifelse(tmp$taken, "Yes", "No")
  
  n_alive_tmp <- reshape2::melt(n_alive, id.vars = "group")
  n_alive_tmp$variable <- as.numeric(n_alive_tmp$variable)
  colnames(n_alive_tmp) <- c("group", "j", "n_alive")
  tmp <- merge(tmp, data.table(n_alive_tmp), by = c("group", "j"))
  tmp$V1 <- tmp$V1 / tmp$n_alive
  tmp$V1 <- pmin(tmp$V1, 1)
  
  quantiles <- ddply(tmp, .(j, group), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
  colnames(quantiles) <- c("j", "group", "lower", "median", "upper")
  # quantiles[c("lower", "median", "upper")] <- quantiles[c("lower", "median", "upper")]# / n_alive1[quantiles$j]
  quantiles$j <- years[quantiles$j]
  quantiles$taken <- quantiles$j %in% unique(titre_dat$samples)
  quantiles$taken <- ifelse(quantiles$taken, "Yes", "No")
  
  quantiles$tested <- quantiles$j %in% unique(titre_dat$virus)
  quantiles$tested <- ifelse(quantiles$tested, "Yes", "No")
  quantiles
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

#' Identify runs of infections from posterior
#'
#' For each individual and MCMC iteration, uses the infection history MCMC chain and detects runs of consecutive infections.
#' @param inf_chain data table of the infection histories posterior
#' @return a tibble giving the consecutive infection run length, the start and end time of each run, which index the run is (ie., which distinct infection), and the time from the end of the previous run, for each i and sampno
#' @examples
#' \dontrun{
#' identify_run_lengths(inf_chain)
#' }
#' @export
identify_run_lengths <- function(inf_chain) {
  inf_chain %>% 
    arrange(sampno, i, j) %>%
    as_tibble() %>%
    group_by(i, sampno) %>%
    mutate(run_group = cumsum(c(0, diff(x != 1) != 0))) %>%
    filter(x == 1) %>%
    group_by(i, sampno, run_group) %>%
    summarise(run_length = n(),
              start_time = first(j),
              end_time = last(j))%>%
    group_by(i, sampno) %>%
    mutate(infection_index = row_number(),
           time_diff = start_time - lag(end_time,1,default=NA)) %>%
    ungroup()%>%
    select(-run_group)
}

#' Summarize runs of infections from posterior
#'
#' Either takes the output of \code{\link{identify_run_lengths}} and produces summary statistics giving the posterior median and 95% quantiles for the length of each infection run
#' @param inf_chain data table of the infection histories posterior or tibble from \code{\link{identify_run_lengths}}
#' @return a tibble summarizing the infection run lengths for each individual and distinct infection event
#' @examples
#' \dontrun{
#' summarize_run_lengths(inf_chain)
#' summarize_run_lengths(identify_run_lengths(inf_chain))
#' }
#' @export
summarize_run_lengths <- function(infections){
  if(!("run_length" %in% colnames(infections))){
    summary <- identify_run_lengths(infections) 
  } else {
    summary <- infections
  }
  summary %>%
    group_by(i, infection_index) %>% 
    dplyr::summarize(median_run_length=median(run_length),
                     lower95_run_length=quantile(run_length,0.025),
                     upper95_run_length=quantile(run_length,0.975))
}
