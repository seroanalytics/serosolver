# serosolver

`serosolver` is a modelling and inference package that uses a dynamic model to infer antibody dynamics and infection histories from cross-sectional or longitudinal serological data. The model infers individual-level infection histories, historical attack rates, and patterns of antibody dynamics by accounting for cross-reactive antibody responses and measurement error.

## Installation

1. Install [R][r-project]

1. Install development version:

    ```r
	#library(devtools)
	install_github("seroanalytics/serosolver")
	library(serosolver)
    ```
	
## Quick start

Read the [quick start vignette][vignette-doc] to set up and run a simple implementation with a simulation model.

[r-project]: http://cran.r-project.org
[vignette-doc]: https://github.com/seroanalytics/serosolver/blob/master/vignettes/serosolver-quick_start_guide.md

## Adding antibody kinetics mechanisms
To add a new antibody kinetics mechanism, the following functions need to be updated:
`titre_data_fast` in `src/infection_model.cpp`
`infection_history_proposal_gibbs_fast` in `src/proposal.cpp`
`create_posterior_func_fast` in `R/posteriors.R`
