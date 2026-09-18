# Modified by an AI assistant on 2026-09-18 using GPT-5. Added targeted roxygen
# imports for base-package functions reported by R CMD check; no functional
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
#' @importFrom doParallel registerDoParallel
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 `%+replace%` aes coord_cartesian element_blank element_line element_rect
#'   element_text facet_grid facet_wrap geom_density geom_histogram geom_hline
#'   geom_line geom_point geom_pointrange geom_rect geom_ribbon geom_violin
#'   geom_vline ggplot ggtitle guide_colourbar guides scale_alpha_continuous
#'   scale_color_manual scale_color_viridis_d scale_fill_gradient scale_fill_manual
#'   scale_fill_viridis_d scale_linetype_manual scale_shape_manual
#'   scale_x_continuous scale_y_continuous stat_density_2d theme_bw theme_classic
#'   theme_minimal theme xlab xlim ylab ylim margin
#' @importFrom grid unit
#' @importFrom graphics legend lines
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom stats cov dbeta density dnorm lag lm median optim pnorm predict
#'   qnorm quantile rbeta rbinom rlnorm rnorm rpois runif sd smooth.spline time
#' @importFrom utils head read.csv
NULL
