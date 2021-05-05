context("Plotting functions")

library(serosolver)
test_that("Check that infection history plots run without error", {
    
    ## Load in exaple data
    data(example_inf_chain)
    data(example_antigenic_map)
    data(example_titre_dat)
    data(example_inf_hist)

    strain_isolation_times <- example_antigenic_map$inf_times

    ## Setup known attack rates
    n_alive <- get_n_alive(example_titre_dat, strain_isolation_times)
    n_infs <- colSums(example_inf_hist)
    known_ar <- n_infs/n_alive
    known_ar <- data.frame("j"=strain_isolation_times,"AR"=known_ar,"group"=1)

    ## Setup known infection histories
    known_inf_hist <- data.frame(example_inf_hist)
    colnames(known_inf_hist) <- strain_isolation_times

    n_alive_group <- get_n_alive_group(example_titre_dat, strain_isolation_times,melt_dat = TRUE)
    n_alive_group$j <- strain_isolation_times[n_alive_group$j]

    ## MCMC chain, time
    p1 <- plot_infection_history_chains_time(example_inf_chain, 0,
                                             sample(1:length(strain_isolation_times),10),
                                             n_alive,TRUE)
    ## MCMC chain, individual
    p2 <- plot_infection_history_chains_indiv(example_inf_chain, 0, 1:10, FALSE)
    ## Total number of infections for each individual
    p3 <- plot_number_infections(example_inf_chain, FALSE)
    ## Total number of infections overall
    p4 <- plot_total_number_infections(example_inf_chain)


    results <- calculate_infection_history_statistics(example_inf_chain, 0, strain_isolation_times,
                                                      n_alive=n_alive_group, known_ar=known_ar,
                                                      known_infection_history=known_inf_hist,
                                                      pad_chain = FALSE
                                                      )
    ## Names expected, by_time_trace, by_indiv_trace, indiv_infections, estimates
    
    all_plots <- plot_posteriors_infhist(example_inf_chain, strain_isolation_times,
                                         n_alive_group, known_ar=known_ar,
                                         known_infection_history =known_inf_hist,
                                         samples=100)
    
    
    ## Expect all AR inference was correct
    expect_equal(sum(results$by_year$correct), length(strain_isolation_times))
    ## I know that 48/50 individual infection histories were correct
    expect_equal(sum(results$by_indiv$correct), 48)

   
    
})


test_that("Check get_titre_predictions returns the correct values and data types when using the bas emodel", {
    data(example_theta_chain)
    data(example_inf_chain)
    data(example_titre_dat)
    data(example_antigenic_map)
    data(example_par_tab)


    ## Check baseline results
    y <- get_titre_predictions(example_theta_chain, example_inf_chain, example_titre_dat,
                               1:5, example_antigenic_map,add_residuals=TRUE,
                               example_par_tab,expand_titredat = FALSE)

    ## Check residual results
    residual_results <- get_titre_predictions(example_theta_chain, example_inf_chain,
                                              example_titre_dat,
                                              1:5,
                                              example_antigenic_map,
                                              for_res_plot=TRUE,
                                              example_par_tab,expand_titredat = FALSE)

    ## Check expanded version works
    expanded_results <- get_titre_predictions(example_theta_chain, example_inf_chain,
                                               example_titre_dat,
                                              1:5,
                                              example_antigenic_map,
                                              for_res_plot=FALSE,
                                              example_par_tab,
                                              expand_titredat = TRUE)
    
    ## Check that infection history fit plot works
    model_fit_plot <- plot_infection_histories(example_theta_chain, example_inf_chain,
                                               example_titre_dat, 1:10, example_antigenic_map,
                                               example_par_tab)
})


test_that("Check that theta diagnsotics work", {
    data(example_par_tab)
    data(example_theta_chain)

    res <- plot_posteriors_theta(example_theta_chain,example_par_tab,samples=100)

    
})


generate_cumulative_inf_plots(example_inf_chain, 0, indivs, example_inf_hist, NULL, times,
                              ages=ages, number_col=2,pad_chain=FALSE)
