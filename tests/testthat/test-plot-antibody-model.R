# Modified by an AI assistant on 2026-09-16 using GPT-5. Added a regression
# test for the return structure of `plot_antibody_model()`.

context("Antibody model plotting")

test_that("plot_antibody_model returns only the longitudinal plot without a map", {
  data(example_par_tab)
  data(example_antigenic_map)

  longitudinal_plot <- plot_antibody_model(example_par_tab, times = 1:4)
  mapped_plots <- plot_antibody_model(
    example_par_tab,
    times = example_antigenic_map$inf_times[1:4],
    antigenic_map = example_antigenic_map
  )

  expect_s3_class(longitudinal_plot, "ggplot")
  expect_type(mapped_plots, "list")
  expect_length(mapped_plots, 2)
})
