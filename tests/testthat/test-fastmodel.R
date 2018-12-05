context("Fast model compare")

library(serosolver)

test_that("Fast solver returns the same titres as the original version, base case", {
    set.seed(1)
    titreDat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
    antigenicMap <- read.csv("../testdata/fonville_annual_continuous.csv")
    strain_isolation_times <- antigenicMap$inf_years
    infection_history_mat <- setup_infection_histories_new(titreDat,strain_isolation_times,5,2)

    parTab <- read.csv("../testdata/parTab_base.csv",stringsAsFactors=FALSE)
    parTab <- parTab[parTab$names != "lambda",]

    f_slow <- create_posterior_func(parTab, titreDat, antigenicMap, 1, function_type=3)
    f_fast <-  create_posterior_func_fast(parTab, titreDat, antigenicMap, 1, function_type=3)
    
    par_names <- parTab$names
    pars <- parTab$values
    names(pars) <- par_names

    y_slow <- f_slow(pars, infection_history_mat)
    y_fast <- f_fast(pars, infection_history_mat)
    expect_equal(y_slow, y_fast)
})


test_that("Fast solver returns the same likelihood as the original version, base case", {
    set.seed(1)
    titreDat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
    antigenicMap <- read.csv("../testdata/fonville_annual_continuous.csv")
    strain_isolation_times <- antigenicMap$inf_years
    infection_history_mat <- setup_infection_histories_new(titreDat,strain_isolation_times,5,2)

    parTab <- read.csv("../testdata/parTab_base.csv",stringsAsFactors=FALSE)
    parTab <- parTab[parTab$names != "lambda",]

    f_slow <- create_posterior_func(parTab, titreDat, antigenicMap, 1, function_type=1)
    f_fast <- create_posterior_func_fast(parTab, titreDat, antigenicMap, 1, function_type=1)    

    par_names <- parTab$names
    pars <- parTab$values
    names(pars) <- par_names

    y_slow <- f_slow(pars, infection_history_mat)
    y_fast <- f_fast(pars, infection_history_mat)
    expect_equal(y_slow, y_fast)
})


