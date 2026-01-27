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

#' Formatted quantiles
#'
#' Given a vector of MCMC samples, generates and formats the desired quantile estimates
#' @param x the vector to summarise
#' @param sig_f how many significant figures to print
#' @param qs the vector of quantiles
#' @param as_text if TRUE, formats nicely as text rather than a vector of numbers
#' @return the formatted quantiles
#' @examples
#' data(example_theta_chain)
#' x <- example_theta_chain$boost_long
#' generate_quantiles(x)
#' @export
generate_quantiles <- function(x, sig_f = 3, qs = c(0.025, 0.5, 0.975), as_text = TRUE) {
  res <- signif(quantile(x, qs), sig_f)
  if (as_text) {
    res <- paste(res[2], " (", res[1], "-", res[3], ")", sep = "")
  }
  return(res)
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
  data.table::setkey(inf_chain, "samp_no", "chain_no")
  ## For each individual, how many infections did they have in each sample in total?
  n_inf_chain <- inf_chain[, list(total_infections = sum(x)), by = key(inf_chain)]
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
summarize_run_lengths <- function(inf_chain){
  if(!("run_length" %in% colnames(inf_chain))){
    summary <- identify_run_lengths(inf_chain) 
  } else {
    summary <- inf_chain
  }
  summary %>%
    group_by(i, infection_index) %>% 
    dplyr::summarize(median_run_length=median(run_length),
                     lower95_run_length=quantile(run_length,0.025),
                     upper95_run_length=quantile(run_length,0.975))
}

#' Identify runs of infections from posterior
#'
#' For each individual and MCMC iteration, uses the infection history MCMC chain and detects runs of consecutive infections.
#' @param inf_chain data table of the infection histories posterior
#' @return a tibble giving the consecutive infection run length, the start and end time of each run, which index the run is (ie., which distinct infection), and the time from the end of the previous run, for each i and samp_no
#' @examples
#' \dontrun{
#' identify_run_lengths(inf_chain)
#' }
#' @export
identify_run_lengths <- function(inf_chain) {
  inf_chain %>% 
    arrange(samp_no, i, j) %>%
    as_tibble() %>%
    group_by(i, samp_no) %>%
    mutate(run_group = cumsum(c(0, diff(x != 1) != 0))) %>%
    filter(x == 1) %>%
    group_by(i, samp_no, run_group) %>%
    summarise(run_length = n(),
              start_time = first(j),
              end_time = last(j))%>%
    group_by(i, samp_no) %>%
    mutate(infection_index = row_number(),
           time_diff = start_time - lag(end_time,1,default=NA)) %>%
    ungroup()%>%
    select(-run_group)
}