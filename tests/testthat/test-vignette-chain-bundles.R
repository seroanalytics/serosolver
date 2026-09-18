# Modified by an AI assistant on 2026-09-18 using GPT-5. Added integrity checks
# for the losslessly compressed MCMC bundles shipped with the vignettes.

test_that("packaged vignette chain bundles have the expected structure", {
  source_bundle_dir <- testthat::test_path(
    "..", "..", "inst", "extdata", "vignette_chains"
  )
  installed_bundle_dir <- system.file(
    "extdata", "vignette_chains", package = "serosolver"
  )
  if (file.exists(file.path(source_bundle_dir, "guide.rds"))) {
    bundle_dir <- source_bundle_dir
  } else {
    bundle_dir <- installed_bundle_dir
  }

  expected <- list(
    advanced_features = c(
      "fixed", "starting", "offset", "exponential", "multiple_biomarker",
      "prior_version_1"
    ),
    cs1_hong_kong = c("prior", "main", "simulation"),
    cs2_fluscape = c("prior", "main", "simulation"),
    demographics_covariates = c("prior", "main", "variant"),
    guide = "guide"
  )

  for (bundle_name in names(expected)) {
    bundle <- readRDS(file.path(bundle_dir, paste0(bundle_name, ".rds")))
    expect_true(all(expected[[bundle_name]] %in% names(bundle)))

    for (fit_name in expected[[bundle_name]]) {
      fit <- bundle[[fit_name]]
      expect_true(all(c("settings", "full_theta", "full_inf") %in% names(fit)))
      expect_length(fit$full_theta, 3)
      expect_length(fit$full_inf, 3)
      for (chain_index in seq_along(fit$full_theta)) {
        theta <- fit$full_theta[[chain_index]]
        inf <- fit$full_inf[[chain_index]]
        expect_gt(nrow(theta), 0)
        expect_gt(nrow(inf), 0)
        expect_true("samp_no" %in% names(theta))
        expect_true(all(c("i", "j", "x", "samp_no") %in% names(inf)))
        expect_equal(length(unique(theta$samp_no)), nrow(theta))
        expect_gt(length(unique(inf$i)), 0)
        expect_gt(length(unique(inf$j)), 0)
        expect_gt(length(unique(inf$samp_no)), 0)
      }
    }
  }
})
