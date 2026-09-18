# Modified by an AI assistant on 2026-09-18 using GPT-5. Added a small
# vignette-only loader for packaged full-length MCMC chain bundles.

load_vignette_chain_bundle <- function(vignette_dir, bundle_file) {
  source_file <- file.path(
    vignette_dir, "..", "inst", "extdata", "vignette_chains", bundle_file
  )
  installed_file <- system.file(
    "extdata", "vignette_chains", bundle_file, package = "serosolver"
  )
  bundle_path <- if (file.exists(source_file)) source_file else installed_file
  if (!nzchar(bundle_path) || !file.exists(bundle_path)) {
    return(NULL)
  }
  readRDS(bundle_path)
}

materialize_vignette_chain_fit <- function(bundle_fit, cache_dir, file_stem) {
  if (is.null(bundle_fit)) {
    return(cache_dir)
  }
  if (!all(c("settings", "full_theta", "full_inf") %in% names(bundle_fit))) {
    stop("The vignette chain bundle is missing required chain data.")
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(bundle_fit$full_theta)) {
    utils::write.csv(
      bundle_fit$full_theta[[i]],
      file.path(cache_dir, paste0(file_stem, "_", i, "_chain.csv")),
      row.names = FALSE
    )
    data.table::fwrite(
      bundle_fit$full_inf[[i]],
      file.path(cache_dir, paste0(file_stem, "_", i, "_infection_histories.csv"))
    )
  }
  serosolver_settings <- bundle_fit$settings
  save(
    serosolver_settings,
    file = file.path(cache_dir, paste0(file_stem, "_serosolver_settings.RData"))
  )
  cache_dir
}
