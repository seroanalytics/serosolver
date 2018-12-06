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

test_that("Inf hist proposal returns the correct post proposal probability after one run, base case", {
    titreDat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
    antigenicMap <- read.csv("../testdata/fonville_annual_continuous.csv")
    strain_isolation_times <- antigenicMap$inf_years
    infection_history_mat <- setup_infection_histories_new(titreDat,strain_isolation_times,5,2)

    parTab <- read.csv("../testdata/parTab_base.csv",stringsAsFactors=FALSE)
    parTab <- parTab[parTab$names != "lambda",]

    par_names <- parTab$names
    pars <- parTab$values
    names(pars) <- par_names
    
    f_slow <- create_posterior_func(parTab, titreDat, antigenicMap, 1, function_type=2)
    f_slow_prob <- create_posterior_func(parTab, titreDat, antigenicMap, 1, function_type=1)
    f_fast <- create_posterior_func_fast(parTab, titreDat, antigenicMap, 1, function_type=2) 
    f_fast_prob<- create_posterior_func_fast(parTab, titreDat, antigenicMap, 1, function_type=1)

    probs_fast <- f_fast_prob(pars, infection_history_mat)
    probs_slow <- f_slow_prob(pars, infection_history_mat)

    indivs <- unique(titreDat$individual)
    
    set.seed(1)
    res_fast <- f_fast(pars, infection_history_mat, probs_fast,
                       indivs, 1, 1, rep(10, length(indivs)),
                       0.5, 2, 1)

    new_probs_fast <- res_fast$old_probs
    new_inf_hist_fast <- res_fast$new_infection_history

    set.seed(1)
    new_inf_hist_slow <- f_slow(pars, infection_history_mat,
                                1, 1,
                                1, rep(10, length(indivs)),
                                0.5, 2, 1)
    new_probs_slow <- f_slow_prob(pars, new_inf_hist_fast)

    expect_equal(new_probs_fast, new_probs_slow)
})

