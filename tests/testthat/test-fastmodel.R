context("Fast model compare")

library(serosolver)

test_that("Fast solver returns the same titres as the original version, base case", {
  set.seed(1)
  titre_dat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
  antigenic_map <- read.csv("../testdata/fonville_annual_continuous.csv")
  strain_isolation_times <- antigenic_map$inf_years
  infection_history_mat <- setup_infection_histories_new(titre_dat, strain_isolation_times, 5, 2)

  par_tab <- read.csv("../testdata/par_tab_base.csv", stringsAsFactors = FALSE)
  par_tab <- par_tab[par_tab$names != "lambda", ]

  f_slow <- create_posterior_func(par_tab, titre_dat, antigenic_map, 1, function_type = 3)
  f_fast <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1, function_type = 3)

  par_names <- par_tab$names
  pars <- par_tab$values
  names(pars) <- par_names

  y_slow <- f_slow(pars, infection_history_mat)
  y_fast <- f_fast(pars, infection_history_mat)

  expect_equal(y_slow, y_fast)
})


test_that("Fast solver returns the same likelihood as the original version, base case", {
  set.seed(1)
  titre_dat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
  antigenic_map <- read.csv("../testdata/fonville_annual_continuous.csv")
  strain_isolation_times <- antigenic_map$inf_years
  infection_history_mat <- setup_infection_histories_new(titre_dat, strain_isolation_times, 5, 2)

  par_tab <- read.csv("../testdata/par_tab_base.csv", stringsAsFactors = FALSE)
  par_tab <- par_tab[par_tab$names != "lambda", ]

  f_slow <- create_posterior_func(par_tab, titre_dat, antigenic_map, 1, function_type = 1)
  f_fast <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1, function_type = 1)

  par_names <- par_tab$names
  pars <- par_tab$values
  names(pars) <- par_names

  y_slow <- f_slow(pars, infection_history_mat)
  y_fast <- f_fast(pars, infection_history_mat)
  expect_equal(y_slow, y_fast)
})

test_that("Inf hist proposal returns the correct post proposal probability after one run, base case", {
  titre_dat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
  antigenic_map <- read.csv("../testdata/fonville_annual_continuous.csv")
  strain_isolation_times <- antigenic_map$inf_years
  infection_history_mat <- setup_infection_histories_new(titre_dat, strain_isolation_times, 5, 2)

  par_tab <- read.csv("../testdata/par_tab_base.csv", stringsAsFactors = FALSE)
  par_tab <- par_tab[par_tab$names != "lambda", ]

  par_names <- par_tab$names
  pars <- par_tab$values
  names(pars) <- par_names

  f_slow <- create_posterior_func(par_tab, titre_dat, antigenic_map, 1, function_type = 2)
  f_slow_prob <- create_posterior_func(par_tab, titre_dat, antigenic_map, 1, function_type = 1)
  f_fast <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1, function_type = 2)
  f_fast_prob <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1, function_type = 1)

  probs_fast <- f_fast_prob(pars, infection_history_mat)
  probs_slow <- f_slow_prob(pars, infection_history_mat)

  indivs <- unique(titre_dat$individual)

  set.seed(1)
  res_fast <- f_fast(
    pars, infection_history_mat, probs_fast,
    indivs, 1, 1, rep(10, length(indivs)),
    0.5, 2, 1
  )

  new_probs_fast <- res_fast$old_probs
  new_inf_hist_fast <- res_fast$new_infection_history

  set.seed(1)
  new_inf_hist_slow <- f_slow(
    pars, infection_history_mat,
    1, 1,
    1, rep(10, length(indivs)),
    0.5, 2, 1
  )
  new_probs_slow <- f_slow_prob(pars, new_inf_hist_fast)

  expect_equal(new_probs_fast, new_probs_slow)
})


test_that("Fast solver returns the same titres as the original version, measurement bias", {
  set.seed(1)
  titre_dat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
  antigenic_map <- read.csv("../testdata/fonville_annual_continuous.csv")
  strain_isolation_times <- antigenic_map$inf_years
  infection_history_mat <- setup_infection_histories_new(titre_dat, strain_isolation_times, 5, 2)

  par_tab <- read.csv("../testdata/par_tab_base.csv", stringsAsFactors = FALSE)
  par_tab <- par_tab[par_tab$names != "lambda", ]
  measurement_bias <- rnorm(15,0,1)
                                        #measurement_bias[15] <- 0
  clusters <- read.csv("../testdata/fonville_clusters.csv")
  n_clusters <- length(unique(clusters$cluster1))
  measurement_indices <- clusters$cluster1
  measurement_indices <- rep(measurement_indices, each=1)

  for(i in 1:length(measurement_bias)){
      tmp <- data.frame(names="rho",values=1,fixed=0,steps=0.1,lower_bound=-10,
                        upper_bound=10,lower_start=-2,upper_start=2, type=3)
      par_tab <- rbind(par_tab, tmp)
  }
  par_tab[par_tab$names == "rho","values"] <- measurement_bias

  
  f_slow <- create_posterior_func(par_tab, titre_dat, antigenic_map, 1,
                                  measurement_indices = measurement_indices,function_type = 3)
  f_fast <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1,
                                       measurement_indices = measurement_indices, function_type = 3)
  f_fast_no_meas <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1,
                                       measurement_indices = NULL, function_type = 3)

  par_names <- par_tab$names
  pars <- par_tab$values
  names(pars) <- par_names

  y_slow <- f_slow(pars, infection_history_mat)
  y_fast <- f_fast(pars, infection_history_mat)
  y_no_meas <- f_fast_no_meas(pars, infection_history_mat)
  expect_false(isTRUE(all.equal(y_fast, y_no_meas)))
  expect_equal(y_slow, y_fast)
})


test_that("Fast solver returns the same titres as the original version, alternative waning", {
  set.seed(1)
  titre_dat <- read.csv("../testdata/fluscape_sim_annual_dat.csv")
  antigenic_map <- read.csv("../testdata/fonville_annual_continuous.csv")
  strain_isolation_times <- antigenic_map$inf_years
  infection_history_mat <- setup_infection_histories_new(titre_dat, strain_isolation_times, 5, 2)

  par_tab <- read.csv("../testdata/par_tab_base.csv", stringsAsFactors = FALSE)
  par_tab <- par_tab[par_tab$names != "lambda", ]
  par_tab[par_tab$names == "wane_type","values"] <- 1
  par_tab[par_tab$names == "wane","values"] <- 0.1
  par_tab[par_tab$names == "kappa","values"] <- 0.9
  par_tab[par_tab$names == "t_change","values"] <- 6
  
  f_slow <- create_posterior_func(par_tab, titre_dat, antigenic_map, 1,
                                  function_type = 3)
  f_fast <- create_posterior_func_fast(par_tab, titre_dat, antigenic_map, 1,
                                       function_type = 3)

  par_names <- par_tab$names
  pars <- par_tab$values
  names(pars) <- par_names

  y_slow <- f_slow(pars, infection_history_mat)
  y_fast <- f_fast(pars, infection_history_mat)
  expect_equal(y_slow, y_fast)
})

