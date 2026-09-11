# Modified by an AI assistant on 2026-09-11 using GPT-5. Added package-level
# roxygen imports for functions used from imported dependencies; no functional
# code was changed.
.datatable.aware <- TRUE
#' @useDynLib serosolver
#' @importFrom Rcpp evalCpp
#' @importFrom data.table `%like%` data.table fsetdiff key
#' @importFrom dplyr `%>%` across all_of arrange bind_cols bind_rows
#'   case_when cur_group_id distinct everything filter group_by if_else
#'   left_join mutate n pull rename row_number sample_n select summarise
#'   summarize tally tibble ungroup
#' @importFrom tidyr complete drop_na expand_grid pivot_longer pivot_wider
#'   unnest
#' @importFrom coda as.mcmc as.mcmc.list effectiveSize gelman.diag
#' @importFrom foreach foreach `%do%`
NULL
