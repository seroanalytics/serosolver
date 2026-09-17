# Modified by an AI assistant on 2026-09-17 using GPT-5. Added fast checks for
# advanced simulation inputs, parameter helpers, and fixed/free phi loading.

context("Feature interfaces")

test_that("simulate_data handles multiple numeric observation types", {
  data(example_par_tab)

  par_tab <- extend_par_tab_biomarker_groups(example_par_tab, 2)
  set.seed(2)
  simulated <- simulate_data(
    par_tab = par_tab,
    n_indiv = 6,
    possible_exposure_times = 1:5,
    measured_biomarker_ids = c(4, 5),
    sampling_times = 1:5,
    nsamps = 2,
    attack_rates = simulate_attack_rates(1:5, mean_par = 0.1, sd_par = 0),
    data_type = c("discrete", "continuous")
  )

  discrete <- simulated$antibody_data$biomarker_group == 1
  continuous <- simulated$antibody_data$biomarker_group == 2
  expect_true(all(simulated$antibody_data$measurement[discrete] ==
                    floor(simulated$antibody_data$measurement[discrete])))
  expect_true(any(simulated$antibody_data$measurement[continuous] !=
                    floor(simulated$antibody_data$measurement[continuous])))
})

test_that("observation-type labels map to the existing numeric codes", {
  expect_equal(
    serosolver:::normalize_data_type(
      c("discrete", "continuous", "false_positive")
    ),
    c(1L, 2L, 3L)
  )
  expect_equal(serosolver:::normalize_data_type("continuous", 2), c(2L, 2L))
  expect_equal(serosolver:::normalize_data_type(2, 2), c(2L, 2L))
  expect_error(
    serosolver:::normalize_data_type(c("continuous", "unknown")),
    "data_type"
  )
  expect_error(
    serosolver:::normalize_data_type(c("continuous", "discrete"), 3),
    "one value per biomarker_group"
  )
})

test_that("coefficient_values sets generated stratification coefficients", {
  data(example_par_tab)

  par_tab <- example_par_tab
  par_tab$stratification <- NA_character_
  par_tab[par_tab$names == "boost_short", "stratification"] <- "urban"
  demographics <- data.frame(
    individual = 1:10,
    birth = 0,
    urban = rep(0:1, each = 5)
  )
  coefficient_values <- data.frame(
    parameter = "boost_short",
    stratification = "urban",
    stratification_level = 1,
    biomarker_group = 1,
    value = 0.75
  )

  set.seed(1)
  simulated <- simulate_data(
    par_tab = par_tab,
    n_indiv = 10,
    possible_exposure_times = 1:5,
    measured_biomarker_ids = 5,
    sampling_times = 1:5,
    nsamps = 2,
    demographics = demographics,
    attack_rates = simulate_attack_rates(1:5, mean_par = 0.1, sd_par = 0),
    data_type = 2,
    coefficient_values = coefficient_values
  )

  coefficient <- simulated$par_tab[
    simulated$par_tab$names == "boost_short_biomarker_1_coef_urban_1", ,
    drop = FALSE
  ]
  expect_equal(nrow(coefficient), 1)
  expect_equal(coefficient$values, 0.75)
  expect_equal(coefficient$par_type, 4)
})

test_that("coefficient_values supports additive covariates and biomarker groups", {
  data(example_par_tab)

  par_tab <- extend_par_tab_biomarker_groups(example_par_tab, 2)
  par_tab$stratification <- NA_character_
  par_tab[par_tab$names == "boost_short", "stratification"] <- "urban, location"
  demographics <- data.frame(
    individual = 1:12,
    birth = 0,
    urban = rep(0:1, each = 6),
    location = rep(0:1, 6)
  )
  coefficient_values <- expand.grid(
    parameter = "boost_short",
    stratification = c("urban", "location"),
    stratification_level = 1,
    biomarker_group = 1:2,
    value = c(0.4, -0.2, 0.8, -0.6)
  )
  coefficient_values <- coefficient_values[
    c(1, 3, 2, 4),
    c("parameter", "stratification", "stratification_level",
      "biomarker_group", "value")
  ]

  set.seed(3)
  simulated <- simulate_data(
    par_tab = par_tab,
    n_indiv = 12,
    possible_exposure_times = 1:5,
    measured_biomarker_ids = c(4, 5),
    sampling_times = 1:5,
    nsamps = 2,
    demographics = demographics,
    attack_rates = simulate_attack_rates(1:5, mean_par = 0.1, sd_par = 0),
    data_type = c(1, 2),
    coefficient_values = coefficient_values
  )

  expected_names <- c(
    "boost_short_biomarker_1_coef_urban_1",
    "boost_short_biomarker_1_coef_location_1",
    "boost_short_biomarker_2_coef_urban_1",
    "boost_short_biomarker_2_coef_location_1"
  )
  coefficient_rows <- simulated$par_tab[
    simulated$par_tab$names %in% expected_names,
    c("names", "values"), drop = FALSE
  ]
  expected_values <- setNames(coefficient_values$value, expected_names)

  expect_equal(nrow(coefficient_rows), 4)
  expect_equal(
    setNames(coefficient_rows$values, coefficient_rows$names),
    expected_values
  )
  expect_false(any(grepl("_coef_(urban|location)_0$", simulated$par_tab$names)))
  expect_equal(
    simulated$par_tab$values[simulated$par_tab$names == "boost_short"],
    c(2, 2)
  )
})

test_that("coefficient_values rejects duplicate and unmatched specifications", {
  data(example_par_tab)
  par_tab <- example_par_tab
  par_tab$stratification <- NA_character_
  par_tab[par_tab$names == "boost_short", "stratification"] <- "urban"
  demographics <- data.frame(individual = 1:4, birth = 0, urban = 0:1)
  antibody_data <- data.frame(
    individual = rep(1:4, each = 2),
    biomarker_id = 5,
    biomarker_group = 1,
    sample_time = 1,
    birth = 0,
    measurement = 1,
    repeat_number = 1,
    urban = rep(0:1, each = 2)
  )
  generated <- add_scale_pars(par_tab, antibody_data, demographics)
  specification <- data.frame(
    parameter = "boost_short",
    stratification = "urban",
    stratification_level = 1,
    biomarker_group = 1,
    value = 0.5
  )

  expect_error(
    serosolver:::apply_coefficient_values(
      generated, rbind(specification, specification)
    ),
    "duplicate"
  )
  specification$parameter <- "not_a_parameter"
  expect_error(
    serosolver:::apply_coefficient_values(generated, specification),
    "do not match"
  )
})

test_that("starting-level and measurement-offset helpers return aligned tables", {
  data(example_antibody_data)
  data(example_par_tab)

  starts <- create_start_level_data(example_antibody_data, "median")
  rho_setup <- add_rhos_par_tab(example_par_tab, c(1968, 1970), n_obs_types = 2)

  expect_true(all(c("individual", "biomarker_id", "biomarker_group",
                    "starting_level") %in% names(starts)))
  expect_equal(nrow(rho_setup[[2]]), 4)
  expect_equal(nrow(rho_setup[[1]][rho_setup[[1]]$names == "rho", ]), 4)
  expect_equal(sort(unique(rho_setup[[2]]$biomarker_group)), c(1, 2))
})

test_that("the antibody model exposes distinct linear and exponential options", {
  pars <- c(
    boost_long = 2, boost_short = 3, boost_delay = 1,
    wane_short = 0.2, wane_long = 0.01, antigenic_seniority = 0,
    cr_long = 0.1, cr_short = 0.03
  )
  antigenic_map <- data.frame(
    inf_times = 1:5,
    x_coord = 1:5,
    y_coord = 1:5
  )

  linear <- simulate_antibody_model(
    pars, infection_history = 1, antigenic_map = antigenic_map,
    times = 1:5, exponential_waning = FALSE
  )
  exponential <- simulate_antibody_model(
    pars, infection_history = 1, antigenic_map = antigenic_map,
    times = 1:5, exponential_waning = TRUE
  )

  expect_equal(nrow(linear), nrow(exponential))
  expect_false(isTRUE(all.equal(linear$antibody_level,
                                exponential$antibody_level)))
})

test_that("fixed and free duplicated phi columns are loaded correctly", {
  location <- tempfile()
  dir.create(location)
  write.csv(
    data.frame(
      samp_no = 1:2,
      phi = c(0.1, 0.2),
      phi.1 = c(0.3, 0.4),
      phi.2 = c(0.5, 0.6),
      posterior_prob = 1:2,
      likelihood = 1:2,
      prior_prob = 1:2
    ),
    file.path(location, "example_chain.csv"),
    row.names = FALSE
  )
  par_tab <- data.frame(
    names = c("phi", "phi", "phi"),
    fixed = c(0, 1, 0)
  )

  chains <- load_theta_chains(
    location, par_tab = par_tab, estimated_only = TRUE,
    convert_mcmc = FALSE, verbose = FALSE
  )

  expect_true(all(c("phi", "phi.2") %in% names(chains$chain)))
  expect_false("phi.1" %in% names(chains$chain))
})
