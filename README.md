# serosolver

`serosolver` is a modelling and inference package that uses a dynamic model to infer antibody dynamics and infection histories from cross-sectional or longitudinal serological data. The model infers individual-level infection histories, historical attack rates, and patterns of antibody dynamics by accounting for cross-reactive antibody responses and measurement error.

## Installation

1. Install [R][r-project]

1. Install development version:

    ```r
	#library(devtools)
	install_github("adamkucharski/serosolver")
	library(serosolver)
    ```
	
## Quick start

Read the [quick start vignette][vignette-doc] to set up and run a simple implementation with a simulation model.

[r-project]: http://cran.r-project.org
[vignette-doc]: https://github.com/adamkucharski/serosolver/blob/master/vignettes/serosolver-quick_start_guide.md
