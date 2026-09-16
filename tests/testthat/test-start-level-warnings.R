# Modified by an AI assistant on 2026-09-16 using GPT-5. Added regression
# checks that automatic starting levels remain silent while incomplete
# user-supplied starting-level tables still warn.

context("Starting-level warnings")

test_that("automatic starting levels do not warn when biomarker IDs are expanded", {
  data(example_theta_chain)
  data(example_inf_chain)
  data(example_antibody_data)
  data(example_antigenic_map)
  data(example_par_tab)

  expect_silent(
    plot_model_fits(
      example_theta_chain,
      example_inf_chain,
      example_antibody_data,
      individuals = 1:3,
      par_tab = example_par_tab,
      antigenic_map = example_antigenic_map,
      expand_to_all_biomarker_ids = TRUE,
      nsamp = 2,
      start_level = "none"
    )
  )
})

test_that("incomplete user-supplied starting levels still warn", {
  data(example_theta_chain)
  data(example_inf_chain)
  data(example_antibody_data)
  data(example_antigenic_map)
  data(example_par_tab)

  antibody_data <- example_antibody_data[example_antibody_data$individual %in% 1:3, ]
  start_levels <- create_start_level_data(antibody_data, "none", FALSE)[1:2, ]

  expect_warning(
    get_antibody_level_predictions(
      example_theta_chain,
      example_inf_chain,
      antibody_data,
      individuals = 1:3,
      antigenic_map = example_antigenic_map,
      par_tab = example_par_tab,
      expand_to_all_biomarker_ids = TRUE,
      nsamp = 2,
      start_level = start_levels
    ),
    "No starting levels were supplied"
  )
})
