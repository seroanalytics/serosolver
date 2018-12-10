context("Infection model")

library(serosolver)

test_that("Infection model returns expected values for base cases", {
  theta <- c(
    "mu" = 2, "mu_short" = 2.7, "tau" = 0.05, "wane" = 0.8, "sigma1" = 0.1, "sigma2" = 0.03, "boost_limit" = 6, "gradient" = 0.2,
    "error" = 1, "wane_type" = 0, "titre_dependent" = 0, "kappa" = 0.9, "t_change" = 6
  )
  sample_time <- 2015
  strain_isolation_times <- 1968:2015
  number_strains <- length(strain_isolation_times)
  antigenic_map <- data.frame(
    inf_years = strain_isolation_times,
    x_coord = seq(1, number_strains),
    y_coord = seq(1, number_strains)
  )
  antigenic_map_melted <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))

  infection_history <- rep(0, number_strains)

  measured_strains <- strain_isolation_times

  measurement_map_indices <- match(measured_strains, strain_isolation_times) - 1
  infection_map_indices <- seq_along(strain_isolation_times) - 1

  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])


  ## Test that returns all 0 when no infections
  y <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    sample_time,
    measurement_map_indices,
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  expect_equal(y, rep(0, number_strains))

  ## Test that returns known values when some infections occur, default parameters
  infection_history[c(10, 20, 30, 40)] <- 1

  correct_titre <- c(
    0, 0, 0.0201010126776668, 0.302943725152286, 0.585786437626905,
    0.868629150101524, 1.15147186257614, 1.43431457505076, 1.71715728752538,
    2, 1.71715728752538, 1.43431457505076, 1.17056782461993, 1.1564256889962,
    1.14228355337246, 1.12814141774873, 1.113999282125, 1.36259884629822,
    1.63129942314911, 1.9, 1.63129942314911, 1.36259884629822, 1.11198918085724,
    1.0978470452335, 1.08370490960977, 1.06956277398604, 1.05542063836231,
    1.29088311754569, 1.54544155877284, 1.8, 1.54544155877284, 1.29088311754569,
    1.05341053709455, 1.03926840147081, 1.02512626584708, 1.01098413022335,
    0.996841994599622, 1.21916738879315, 1.45958369439657, 1.7, 1.45958369439657,
    1.21916738879315, 0.978751083189721, 0.738334777586295, 0.497918471982869,
    0.257502166379443, 0.0170858607760168, 0
  )

  y_with_infections <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    sample_time,
    measurement_map_indices,
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )

  expect_equal(y_with_infections, correct_titre)

  ## Test that we get the same result when only passing actual infection years
  inf_years <- which(infection_history == 1)

  y_indexed <- infection_model_indiv(
    theta,
    infection_history[inf_years],
    strain_isolation_times[inf_years],
    infection_map_indices[inf_years],
    sample_time,
    measurement_map_indices,
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )

  expect_equal(y_indexed, correct_titre)

  ## Equivalent when only asking for a subset of measured strains
  subset_measured <- seq(1, length(measurement_map_indices), by = 2)
  y_measured_subset <- infection_model_indiv(
    theta,
    infection_history[inf_years],
    strain_isolation_times[inf_years],
    infection_map_indices[inf_years],
    sample_time,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  expect_equal(y_measured_subset, correct_titre[subset_measured])

  ## Try with no waning, cross reactivity and seniority
  theta["wane"] <- 1
  theta["sigma1"] <- 1
  theta["sigma2"] <- 1
  theta["tau"] <- 0
  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

  y_base <- infection_model_indiv(
    theta,
    infection_history[inf_years],
    strain_isolation_times[inf_years],
    infection_map_indices[inf_years],
    sample_time,
    measurement_map_indices,
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  correct_base_titres <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,
    0, 0, 0, 0, 0, 0, 0
  )
  expect_equal(y_base, correct_base_titres)
})
test_that("Infection model returns expected values for waning mechanisms", {
  theta <- c(
    "mu" = 2, "mu_short" = 2.7, "tau" = 0.05, "wane" = 0.8, "sigma1" = 0.1, "sigma2" = 0.03, "boost_limit" = 6, "gradient" = 0.2,
    "error" = 1, "wane_type" = 0, "titre_dependent" = 0, "kappa" = 0.9, "t_change" = 6
  )
  sample_time <- 2015
  strain_isolation_times <- 1968:2015
  number_strains <- length(strain_isolation_times)
  antigenic_map <- data.frame(
    inf_years = strain_isolation_times,
    x_coord = seq(1, number_strains),
    y_coord = seq(1, number_strains)
  )
  antigenic_map_melted <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))

  infection_history <- rep(0, number_strains)
  infection_history[c(10, 20, 30, 40)] <- 1
  inf_years <- which(infection_history == 1)

  measured_strains <- strain_isolation_times

  measurement_map_indices <- match(measured_strains, strain_isolation_times) - 1
  infection_map_indices <- seq_along(strain_isolation_times) - 1

  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])
  ## Test that waning is working correctly
  sample_time1 <- strain_isolation_times[40]
  sample_time2 <- strain_isolation_times[41]
  subset_measured <- seq(1, length(measurement_map_indices), by = 2)

  ## First time point
  y_sample_time1 <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    sample_time1,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  correct_sample_time1 <- c(
    0, 0.0201010126776668, 0.585786437626905, 1.15147186257614,
    1.71715728752538, 1.71715728752538, 1.17056782461993, 1.14228355337246,
    1.16952139542909, 1.88155874399197, 2.07629595153075, 1.75172291677765,
    1.91817585306896, 2.08462878936027, 2.76938691730958, 2.96412412484835,
    2.66683031070883, 2.83328324700015, 2.99973618329146, 3.65721509062719,
    3.65721509062719, 2.98164527188156, 2.30607545313593, 1.6305056343903
  )
  expect_equal(y_sample_time1, correct_sample_time1)

  ## Second time point
  y_sample_time2 <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    sample_time2,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  correct_sample_time2 <- c(
    0, 0.0201010126776668, 0.585786437626905, 1.15147186257614,
    1.71715728752538, 1.71715728752538, 1.17056782461993, 1.14228355337246,
    1.12510370478582, 1.68135128731768, 1.72029872882544, 1.23993592804132,
    1.25059909830161, 1.2612622685619, 1.79023063048019, 1.82917807198795,
    1.3760944918174, 1.3867576620777, 1.39742083233799, 1.8991099736427,
    1.8991099736427, 1.37932992092809, 0.859549868213481, 0.339769815498874
  )
  expect_equal(y_sample_time2, correct_sample_time2)

  ## Waning type 1
  theta["wane_type"] <- 1
  theta["wane"] <- 0.1
  y_waning <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    2010,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  correct_titre_waning_2010 <- c(
    0.331365264754357, 0.431308532522922, 1.07683621256306, 1.74176397219297,
    2.45533512357313, 2.57343751238047, 2.11516710065878, 2.17520188059508,
    2.27410213984424, 3.01603737732925, 3.19061355632283, 2.79582057555722,
    2.89205356583602, 2.98828655611483, 3.60282473805162, 3.65929852823787,
    3.1556978254058, 3.13713559329928, 3.12706497852162, 3.60802033478837,
    3.48991794598103, 2.7546669006897, 2.05520232373987, 1.37001019591586
  )


  expect_equal(y_waning, correct_titre_waning_2010)

  theta["wane_type"] <- 0
  y_waning_old <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    2010,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  expect_false(isTRUE(all.equal(y_waning_old, y_waning)))
})

test_that("Infection model returns expected values for titre dependent boosting", {
  theta <- c(
    "mu" = 2, "mu_short" = 2.7, "tau" = 0.05, "wane" = 0.8, "sigma1" = 0.1, "sigma2" = 0.03, "boost_limit" = 6, "gradient" = 0.2,
    "error" = 1, "wane_type" = 0, "titre_dependent" = 0, "kappa" = 0.9, "t_change" = 6
  )
  sample_time <- 2015
  strain_isolation_times <- 1968:2015
  number_strains <- length(strain_isolation_times)
  antigenic_map <- data.frame(
    inf_years = strain_isolation_times,
    x_coord = seq(1, number_strains),
    y_coord = seq(1, number_strains)
  )
  antigenic_map_melted <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))

  infection_history <- rep(0, number_strains)
  infection_history[c(10, 20, 30, 40)] <- 1
  inf_years <- which(infection_history == 1)

  measured_strains <- strain_isolation_times

  measurement_map_indices <- match(measured_strains, strain_isolation_times) - 1
  infection_map_indices <- seq_along(strain_isolation_times) - 1

  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])
  ## Test titre dependent boosting
  infection_history[15] <- 1
  subset_measured <- seq(1, length(measurement_map_indices), by = 2)
  y_no_titre_dependent <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    sample_time,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  theta["titre_dependent"] <- 1
  y_titre_dependent <- infection_model_indiv(
    theta,
    infection_history,
    strain_isolation_times,
    infection_map_indices,
    sample_time,
    measurement_map_indices[subset_measured],
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  correct_titre_dependent <- c(
    0, 0.0201010126776668, 0.585786437626905, 1.15147186257614,
    1.97123636456396, 2.44567705678503, 2.36862614248912, 2.67682577310048,
    2.03614401926971, 1.94104667771051, 1.46660598548945, 0.826487818075144,
    0.804296464556281, 0.782105111037419, 1.14516430242582, 1.14516430242582,
    0.7805280201005, 0.758336666581637, 0.736145313062774, 1.07780169640077,
    1.07780169640077, 0.722740040098964, 0.367678383797157, 0.0126167274953504
  )
  expect_equal(y_titre_dependent, correct_titre_dependent)
  expect_false(isTRUE(all.equal(y_titre_dependent, y_no_titre_dependent)))
  infection_history[15] <- 0
  theta["titre_dependent"] <- 0
})

test_that("Infection model returns expected values for antigenic seniority", {
  ## Try with no waning, cross reactivity and seniority
  theta <- c(
    "mu" = 2, "mu_short" = 2.7, "tau" = 0.05, "wane" = 0.8, "sigma1" = 0.1, "sigma2" = 0.03, "boost_limit" = 6, "gradient" = 0.2,
    "error" = 1, "wane_type" = 0, "titre_dependent" = 0, "kappa" = 0.9, "t_change" = 6
  )
  sample_time <- 2015
  strain_isolation_times <- 1968:2015
  number_strains <- length(strain_isolation_times)
  antigenic_map <- data.frame(
    inf_years = strain_isolation_times,
    x_coord = seq(1, number_strains),
    y_coord = seq(1, number_strains)
  )
  antigenic_map_melted <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))

  infection_history <- rep(0, number_strains)
  infection_history[c(10, 20, 30, 40)] <- 1
  inf_years <- which(infection_history == 1)

  measured_strains <- strain_isolation_times

  measurement_map_indices <- match(measured_strains, strain_isolation_times) - 1
  infection_map_indices <- seq_along(strain_isolation_times) - 1

  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])
  theta["wane"] <- 1
  theta["sigma1"] <- 1
  theta["sigma2"] <- 1
  theta["tau"] <- 0
  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])


  ## Try with seniority
  theta["tau"] <- 0.1
  y_seniority <- infection_model_indiv(
    theta,
    infection_history[inf_years],
    strain_isolation_times[inf_years],
    infection_map_indices[inf_years],
    sample_time,
    measurement_map_indices,
    antigenic_map_long,
    antigenic_map_short,
    number_strains
  )
  correct_seniority_titres <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.8,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.4,
    0, 0, 0, 0, 0, 0, 0, 0
  )
  expect_equal(y_seniority, correct_seniority_titres)
})
