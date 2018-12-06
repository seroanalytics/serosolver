
#' @export
create_posterior_func_fast <- function(par_tab,
                                       titre_dat,
                                       antigenic_map,
                                       version=1,
                                       solve_likelihood=TRUE,
                                       age_mask=NULL,
                                       measurement_indices_by_time=NULL,
                                       mu_indices=NULL,
                                       n_alive=NULL,
                                       function_type=1,
                                       ...){
   

    ## Seperate out initial readings and repeat readings - we only
    ## want to solve the model once for each unique indiv/sample/virus year tested
    titre_dat_unique <- titre_dat[titre_dat$run == 1,]
    titre_dat_repeats <- titre_dat[titre_dat$run != 1,]
    tmp <- row.match(titre_dat_repeats[,c("individual","samples","virus")], titre_dat_unique[,c("individual","samples","virus")])
    titre_dat_repeats$index <- tmp
    
    ## Setup data vectors and extract
    setup_dat <- setup_titredat_for_posterior_func(titre_dat_unique, antigenic_map, age_mask, n_alive)

    individuals <- setup_dat$individuals
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    strain_isolation_times <- setup_dat$strain_isolation_times
    infection_strain_indices <- setup_dat$infection_strain_indices
    sample_times <- setup_dat$sample_times
    rows_per_indiv_in_samples <- setup_dat$rows_per_indiv_in_samples
    nrows_per_individual_in_data <- setup_dat$nrows_per_individual_in_data
    cum_nrows_per_individual_in_data <- setup_dat$cum_nrows_per_individual_in_data
    nrows_per_blood_sample <- setup_dat$nrows_per_blood_sample
    measured_strain_indices <- setup_dat$measured_strain_indices
    n_alive <- setup_dat$n_alive
    age_mask <- setup_dat$age_mask
    strain_mask <- setup_dat$strain_mask
    n_indiv <- setup_dat$n_indiv
    DOBs <- setup_dat$DOBs

## Some additional setup for the repeat data    
    nrows_per_individual_in_data_repeats <- NULL
    for(individual in unique(individuals)){
        nrows_per_individual_in_data_repeats <- c(nrows_per_individual_in_data_repeats, nrow(titre_dat_repeats[titre_dat_repeats$individual == individual,]))
    }
    cum_nrows_per_individual_in_data_repeats <- cumsum(c(0,nrows_per_individual_in_data_repeats))
    
    titres_unique <- titre_dat_unique$titre
    titres_repeats <- titre_dat_repeats$titre
    repeat_indices <- titre_dat_repeats$index
    repeat_indices_cpp <- repeat_indices - 1

    all_titres <- c(titres_unique, titres_repeats)   
    par_names <- par_tab$names
    
    if(function_type == 1){
        f <- function(pars, infection_history_mat){
            names(pars) <- par_names
            antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma2"])
            y_new <- titre_data_fast(pars, infection_history_mat, strain_isolation_times, infection_strain_indices,
                                     sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
                                     nrows_per_blood_sample, measured_strain_indices, antigenic_map_long,
                                     antigenic_map_short, NULL)
            
            ## Calculate likelihood for unique titres and repeat data
            liks <- likelihood_func_fast(pars, titres_unique, y_new)
            liks_repeats <- likelihood_func_fast(pars, titres_repeats, y_new[repeat_indices])
            
            ## Sum these for each individual
            liks <- sum_buckets(liks, nrows_per_individual_in_data) +
                sum_buckets(liks_repeats, nrows_per_individual_in_data_repeats)
        }
    } else if(function_type == 2){
        f <- function(pars, infection_history_mat,
                      probs, sampled_indivs,
                      alpha, beta,
                      nInfs, swap_propn, swap_dist,
                      temp){
            names(pars) <- par_names

            ## Work out short and long term boosting cross reactivity - C++ function
            antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma2"])
            ## Now pass to the C++ function
            
            res <- infection_history_proposal_gibbs_fast(pars,
                                                         infection_history_mat,
                                                         probs,
                                                         sampled_indivs,
                                                         nInfs,
                                                         age_mask,
                                                         strain_mask,
                                                         n_alive,
                                                         swap_propn,
                                                         swap_dist,
                                                         alpha,
                                                         beta,
                                                         strain_isolation_times,
                                                         infection_strain_indices,
                                                         sample_times,
                                                         rows_per_indiv_in_samples,
                                                         cum_nrows_per_individual_in_data,
                                                         cum_nrows_per_individual_in_data_repeats,
                                                         nrows_per_blood_sample,
                                                         measured_strain_indices, 
                                                         antigenic_map_long,
                                                         antigenic_map_short,
                                                         titres_unique,
                                                         titres_repeats,
                                                         repeat_indices_cpp,
                                                         temp,
                                                         solve_likelihood
                                                         )
            return(res)
        }
    } else {
        f <- function(pars, infection_history_mat){
            names(pars) <- par_names
            antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma2"])
            y_new <- titre_data_fast(pars, infection_history_mat, strain_isolation_times, infection_strain_indices,
                                     sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
                                     nrows_per_blood_sample, measured_strain_indices, antigenic_map_long,
                                     antigenic_map_short, NULL)
        }
    }
    f
}
