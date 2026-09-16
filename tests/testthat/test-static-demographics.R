# Modified by an AI assistant on 2026-09-16 using GPT-5. Added a regression
# test for fixed demographic tables supplied without a `time` column.

context("Static demographic inputs")

test_that("static demographics are not treated as time-varying", {
  data(example_antibody_data)
  data(example_par_tab)

  static_demographics <- unique(example_antibody_data[, c("individual", "birth")])
  static_demographics$urban <- as.integer(
    static_demographics$individual > median(static_demographics$individual)
  )

  par_tab <- example_par_tab
  par_tab$stratification <- NA_character_
  par_tab[par_tab$names == "boost_short", "stratification"] <- "urban"

  par_tab <- add_scale_pars(par_tab, example_antibody_data, static_demographics)
  demographic_groups <- get_demographic_groups(
    par_tab,
    example_antibody_data,
    static_demographics
  )
  aligned <- add_stratifying_variables(
    example_antibody_data,
    static_demographics,
    par_tab,
    demographic_groups$use_demographic_groups
  )

  expect_false(demographic_groups$timevarying_demographics)
  expect_true("urban" %in% names(aligned$antibody_data))
  expect_true("demographic_group" %in% names(aligned$antibody_data))
  expect_false(is.null(aligned$demographics))
})
