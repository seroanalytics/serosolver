
<!-- Modified by an AI assistant on 2026-09-17 using GPT-5. Updated the worked
data_type example to use character observation-model labels. -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# serosolver

<img src="figure/logo.png" align="right" width="250" alt="serosolver logo">

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`serosolver` uses a hierarchical model with a custom Markov chain Monte
Carlo sampler to simultaneously infer antibody kinetics and infection
histories from cross-sectional or longitudinal serological data.
`serosolver` is a [time-since-infection serodynamics
model](https://osf.io/preprints/osf/kqdsn), meaning that infection times
are back-calculated from one or more antibody measurements through an
antibody kinetics model. `serosolver` can be used to infer infection
timings during a study period using longitudinal measurements against a
single antigen, or lifetime infection histories using multi-antigen
serology panels. The package and model are described by Hay *et al.*
[here](https://doi.org/10.1371/journal.pcbi.1007840).

<br clear="right">

## New features

`serosolver` is in active development with new features and ongoing
improvements.
<details>

<summary>

List of new features:
</summary>

- Generalisation to multiple biomarker types per sample (e.g., antibody
  titre and avidity)
- Support for continuous and discrete observations (e.g., ELISA data and
  HAI titres)
- Stratification by demographic variables
- Fixing infection states during fitting
- Fixing or estimating starting titres
- Improved user interface

</details>

## Installation

`serosolver` is in the midst of an overhaul. Please use the `published`
branch to ensure continued compatibility with existing projects.

``` r
remotes::install_github("seroanalytics/serosolver", ref = "published")
```

Install the development version of serosolver from
[GitHub](https://github.com/seroanalytics/serosolver):

``` r
remotes::install_github("seroanalytics/serosolver")
library(serosolver)
```

## Dependencies

A working C++14 compiler is needed. The package uses `Rcpp`,
`RcppArmadillo`, and `RcppParallel`.

``` r
required_packages <- c(
  "data.table", "ggplot2", "dplyr", "tidyr", "Rcpp", "coda",
  "doParallel", "doRNG", "foreach", "Matrix", "MASS", "reshape2",
  "tibble", "RcppArmadillo", "RcppParallel"
)

additional_packages <- c(
  "remotes", "devtools", "plyr", "tidyverse"
)

install.packages(c(required_packages, additional_packages))
```

## Resources

Read the
[guide](https://seroanalytics.github.io/serosolver/articles/serosolver-guide.html)
to set up and run a simple implementation with a simulation model.

Additional vignettes:

- [Longitudinal
  data](https://seroanalytics.github.io/serosolver/articles/cs1_hong_kong.html):
  estimating infection timings using longitudinal data, example of
  influenza A/H1N1p in Hong Kong
- [Cross-sectional
  data](https://seroanalytics.github.io/serosolver/articles/cs2_vignette.html):
  estimating life-course infection histories from multi-strain serology,
  example of influenza A/H3N2 from the [Fluscape
  cohort](https://pubmed.ncbi.nlm.nih.gov/26875566/)
- [Optional
  features](https://seroanalytics.github.io/serosolver/articles/serosolver-guide.html):
  walkthrough of additional `serosolver` features and use cases, such as
  inclusion of biomarker-specific measurement offsets
- [Multiple
  measurements](https://seroanalytics.github.io/serosolver/articles/serosolver-guide.html):
  fitting `serosolver` to multiple biomarker types, example of binding
  avidity and ELISA measurements per sample
- [Group-level
  differences](https://seroanalytics.github.io/serosolver/articles/serosolver-guide.html):
  estimating demographic differences in antibody kinetics and attack
  rates
- [Naming
  conventions](https://seroanalytics.github.io/serosolver/articles/naming_convention.html):
  the current names for datasets, variables, and model inputs

## Example

This is a basic example of simulating some serological data and fitting
the model using the MCMC framework.

``` r
library(serosolver)

## Load in example parameter values and antigenic map
data(example_par_tab)
data(example_antigenic_map)
data(example_antibody_data)
data(example_inf_hist)

## Check the dataset and model control table for errors
example_antibody_data <- check_data(example_antibody_data)
example_par_tab <- check_par_tab(example_par_tab)
```

``` r
plot_antibody_data(example_antibody_data,example_antigenic_map$inf_times,n_indivs=1:5,infection_histories=example_inf_hist)
```

<img src="man/figures/README-example_plot-1.png" alt="" width="100%" />

``` r
## Run serosolver
readme_mcmc_dir <- file.path("inst", "extdata", "readme", "chains")
dir.create(readme_mcmc_dir, recursive = TRUE, showWarnings = FALSE)
output <- serosolver::serosolver(example_par_tab, example_antibody_data, antigenic_map=example_antigenic_map,
                filename=file.path(readme_mcmc_dir, "readme"), n_chains=3,parallel=TRUE,data_type="continuous",
                mcmc_pars=c(adaptive_iterations=10000, iterations=50000),verbose=TRUE)
#> ================================ Running serosolver ================================
#> Requested 3 chains in parallel, setting up parallel session using the parallel package
#> Progress messages will be piped to inst/extdata/readme/chains/readme_log.txt when `parallel` is set to true
#> Model fitting started
#> Model fitting done!
#> Generating MCMC diagnostics
#> Generating output plots
#> ================================ Finished ================================
```

``` r
output$plot_fits_cross_sectional
#> [[1]]
```

<img src="man/figures/README-example_model_fits-1.png" alt="" width="100%" />

## AI assistance

AI assistance was used during development of this package for
documentation drafting, repository audits, code review, and selected
implementation edits. Most of this assistance was provided using
OpenAI’s GPT-5 model. The package author is responsible for reviewing
and approving all changes. The core of `serosolver` remains the same as
the published version; new advanced features were developed and
implemented manually, with some AI assistance used to align complex
coding workflows.
