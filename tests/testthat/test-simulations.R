context("Simulation code")

library(serosolver)

test_that("Test that simulating a longitudinal cohort study works correctly", {
    plot_data(example_titre_dat, example_inf_hist, strain_isolation_times, 5)
})
