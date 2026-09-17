# Modified by an AI assistant on 2026-09-17 using GPT-5. Replaced the obsolete
# exploratory simulation script with focused checks of the current `simulate_data()` output.

context("Simulation code")

test_that("simulate_data returns a usable longitudinal data set", {
  data(example_par_tab)

  set.seed(1)
  simulated <- simulate_data(
    par_tab = example_par_tab,
    n_indiv = 8,
    possible_exposure_times = 1:5,
    measured_biomarker_ids = 5,
    sampling_times = 1:5,
    nsamps = 3,
    age_min = 5,
    age_max = 5,
    attack_rates = simulate_attack_rates(1:5, mean_par = 0.2, sd_par = 0),
    data_type = 2
  )

  expect_true(all(c("antibody_data", "infection_histories", "attack_rates") %in%
                    names(simulated)))
  expect_equal(nrow(simulated$infection_histories), 8)
  expect_equal(ncol(simulated$infection_histories), 5)
  expect_true(all(c("individual", "sample_time", "biomarker_id", "measurement") %in%
                    names(simulated$antibody_data)))
})
