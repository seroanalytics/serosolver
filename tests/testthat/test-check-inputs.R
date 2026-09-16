# Modified by an AI assistant on 2026-09-16 using GPT-5. Added focused tests
# for the user-facing defaults and required-column checks in `check_par_tab()`.

context("Input checks")

test_that("check_par_tab has a usable default for MCMC checks", {
  data(example_par_tab)

  checked <- check_par_tab(example_par_tab, mcmc = TRUE)

  expect_true("steps" %in% names(checked))
  expect_false(any(checked$names == "phi"))
})

test_that("check_par_tab retains the version 1 option", {
  data(example_par_tab)

  checked <- check_par_tab(
    example_par_tab,
    mcmc = TRUE,
    version = 1,
    possible_exposure_times = 1:2
  )

  expect_equal(sum(checked$names == "phi"), 2)
})

test_that("check_par_tab reports missing required columns", {
  data(example_par_tab)
  incomplete <- example_par_tab
  incomplete$values <- NULL

  expect_error(check_par_tab(incomplete), "missing required columns")
})

test_that("fixed parameters do not trigger starting-bound warnings", {
  data(example_par_tab)
  par_tab <- example_par_tab
  par_tab[par_tab$names == "max_measurement", "fixed"] <- 1
  par_tab[par_tab$names == "max_measurement", c("lower_bound", "upper_bound")] <- 11

  expect_silent(check_par_tab(par_tab, mcmc = TRUE))
})

test_that("estimated parameters retain starting-bound warnings", {
  data(example_par_tab)
  par_tab <- example_par_tab
  par_tab[par_tab$names == "boost_long", "lower_bound"] <- 2

  expect_warning(
    check_par_tab(par_tab, mcmc = TRUE),
    "boost_long"
  )
})
