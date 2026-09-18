# Modified by an AI assistant on 2026-09-18 using GPT-5. Added a script for
# rebuilding the losslessly compressed MCMC bundles used by the vignettes.

# Run the vignette MCMC fits first, if needed, then run this script from any
# location with: Rscript inst/extdata/create_vignette_chain_bundles.R

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- script_args[grepl("^--file=", script_args)]
if (length(script_file) > 0L) {
  repo_root <- normalizePath(
    file.path(
      dirname(sub("^--file=", "", script_file[1])), "..", ".."
    ),
    mustWork = TRUE
  )
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}

chain_root <- file.path(repo_root, "inst", "extdata")
bundle_root <- file.path(chain_root, "vignette_chains")
dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)

read_chain_bundle <- function(chain_dir, file_stem) {
  settings_file <- file.path(
    chain_dir, paste0(file_stem, "_serosolver_settings.RData")
  )
  chain_files <- list.files(
    chain_dir,
    pattern = paste0("^", file_stem, "_[0-9]+_chain\\.csv$"),
    full.names = TRUE
  )
  chain_numbers <- as.integer(sub(
    paste0("^", file_stem, "_([0-9]+)_chain\\.csv$"), "\\1",
    basename(chain_files)
  ))
  chain_files <- chain_files[order(chain_numbers)]
  history_files <- sub(
    "_chain\\.csv$", "_infection_histories.csv", chain_files
  )

  if (!file.exists(settings_file) || length(chain_files) != 3L ||
      !all(file.exists(history_files))) {
    stop(
      "Expected three complete chains and settings in ", chain_dir,
      " for file stem '", file_stem, "'.",
      call. = FALSE
    )
  }

  settings_environment <- new.env(parent = emptyenv())
  load(settings_file, envir = settings_environment)
  if (!exists("serosolver_settings", envir = settings_environment,
              inherits = FALSE)) {
    stop("No serosolver_settings object found in ", settings_file, ".",
         call. = FALSE)
  }

  full_theta <- lapply(chain_files, utils::read.csv, check.names = FALSE)
  full_inf <- lapply(history_files, data.table::fread)
  names(full_theta) <- basename(chain_files)
  names(full_inf) <- basename(history_files)

  list(
    settings = settings_environment$serosolver_settings,
    full_theta = full_theta,
    full_inf = full_inf
  )
}

write_vignette_bundle <- function(bundle_name, fits) {
  output_file <- file.path(bundle_root, paste0(bundle_name, ".rds"))
  saveRDS(fits, output_file, compress = "xz")
  message("Wrote ", output_file)
}

## These are the retained vignette fits. The README chains are intentionally
## excluded because the README is expected to run serosolver directly.
write_vignette_bundle("guide", list(
  guide = read_chain_bundle(file.path(chain_root, "guide"), "guide")
))

write_vignette_bundle("cs1_hong_kong", list(
  prior = read_chain_bundle(
    file.path(chain_root, "case_study_1", "prior_mcmc_chains"), "prior"
  ),
  main = read_chain_bundle(
    file.path(chain_root, "case_study_1", "hong_kong_mcmc_chains"), "hong_kong"
  ),
  simulation = read_chain_bundle(
    file.path(chain_root, "case_study_1", "simulation_mcmc_chains"), "simulated"
  )
))

write_vignette_bundle("cs2_fluscape", list(
  prior = read_chain_bundle(
    file.path(chain_root, "case_study_2", "prior_mcmc_chains"), "prior"
  ),
  main = read_chain_bundle(
    file.path(chain_root, "case_study_2", "fluscape_mcmc_chains"), "fluscape"
  ),
  simulation = read_chain_bundle(
    file.path(chain_root, "case_study_2", "simulation_mcmc_chains"), "simulated"
  )
))

write_vignette_bundle("demographics_covariates", list(
  prior = read_chain_bundle(
    file.path(chain_root, "demographics_covariates", "prior_mcmc_chains"),
    "prior_covariates"
  ),
  main = read_chain_bundle(
    file.path(chain_root, "demographics_covariates", "mcmc_chains"),
    "demographic_fit"
  ),
  variant = read_chain_bundle(
    file.path(chain_root, "demographics_covariates", "variant_mcmc_chains"),
    "variant_fit"
  )
))

write_vignette_bundle("advanced_features", list(
  fixed = read_chain_bundle(
    file.path(chain_root, "advanced_features", "fixed_infection_states"),
    "fixed_infection_states"
  ),
  starting = read_chain_bundle(
    file.path(chain_root, "advanced_features", "starting_levels"),
    "starting_levels_one_biomarker"
  ),
  offset = read_chain_bundle(
    file.path(chain_root, "advanced_features", "measurement_bias"),
    "measurement_bias"
  ),
  exponential = read_chain_bundle(
    file.path(chain_root, "advanced_features", "exponential_waning"),
    "exponential_waning"
  ),
  multiple_biomarker = read_chain_bundle(
    file.path(chain_root, "advanced_features", "multiple_biomarker_groups"),
    "multiple_biomarker_groups"
  ),
  prior_version_1 = read_chain_bundle(
    file.path(chain_root, "advanced_features", "prior_version_1"),
    "prior_version_1"
  )
))
