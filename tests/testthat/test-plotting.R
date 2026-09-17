# Modified by an AI assistant on 2026-09-17 using GPT-5. Replaced obsolete
# plotting tests with checks against the current plotting and prediction APIs.

context("Plotting functions")

test_that("infection-history plotting functions return their current outputs", {
  data(example_antibody_data)
  data(example_antigenic_map)
  data(example_inf_chain)
  data(example_inf_hist)

  antibody_data <- example_antibody_data
  antibody_data$population_group <- 1
  possible_exposure_times <- example_antigenic_map$inf_times
  n_alive <- get_n_alive(antibody_data, possible_exposure_times)
  n_alive_group <- get_n_alive_group(
    antibody_data, possible_exposure_times, melt_data = TRUE
  )
  n_alive_group$j <- possible_exposure_times[n_alive_group$j]
  known_ar <- data.frame(
    j = possible_exposure_times,
    AR = colSums(example_inf_hist) / n_alive,
    population_group = 1
  )

  posterior_plots <- plot_infection_history_posteriors(
    example_inf_chain,
    possible_exposure_times,
    n_alive_group,
    known_ar = known_ar,
    known_infection_history = example_inf_hist,
    samples = 5,
    pad_chain = FALSE
  )
  cumulative_plots <- plot_cumulative_infection_histories(
    example_inf_chain,
    indivs = 1:4,
    real_inf_hist = example_inf_hist,
    possible_exposure_times = possible_exposure_times,
    nsamp = 5,
    pad_chain = FALSE
  )

  expect_named(
    posterior_plots,
    c("by_time_trace", "by_indiv_trace", "indiv_infections", "estimates")
  )
  expect_length(cumulative_plots, 2)
  expect_s3_class(posterior_plots$indiv_infections, "ggplot")
  expect_s3_class(cumulative_plots[[1]], "ggplot")
})

test_that("antibody predictions use the current prediction helper", {
  data(example_theta_chain)
  data(example_inf_chain)
  data(example_antibody_data)
  data(example_antigenic_map)
  data(example_par_tab)

  predictions <- get_antibody_level_predictions(
    example_theta_chain,
    example_inf_chain,
    example_antibody_data,
    individuals = 1:3,
    antigenic_map = example_antigenic_map,
    par_tab = example_par_tab,
    nsamp = 3
  )

  expect_named(
    predictions,
    c("predictions", "histories", "best_infhist", "predicted_observations")
  )
  expect_true(all(c("median", "lower", "upper") %in%
                    names(predictions$predicted_observations)))
  expect_gt(nrow(predictions$predicted_observations), 0)
})

test_that("theta posterior plots use the current arguments", {
  data(example_par_tab)
  data(example_theta_chain)

  plots <- plot_posteriors_theta(example_theta_chain, example_par_tab)

  expect_type(plots, "list")
  expect_length(plots, 8)
  expect_true(any(vapply(plots, inherits, logical(1), what = "ggplot")))
})
