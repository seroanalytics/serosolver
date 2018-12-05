row.match <- function (x, table, nomatch = NA) {
  if (class(table) == "matrix") 
    table <- as.data.frame(table)
  if (is.null(dim(x))) 
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}

#' @export
create_posterior_func_fast <- function(parTab,
                                  titreDat,
                                  antigenicMap,
                                  version=1,
                                  solve_likelihood=TRUE,
                                  ageMask=NULL,
                                  measurement_indices_by_time=NULL,
                                  mu_indices=NULL,
                                  n_alive=NULL,
                                  function_type=1,
                                  ...){
    strain_isolation_times <- antigenicMap$inf_years
    number_strains <- length(strain_isolation_times)
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))

    ## Seperate out initial readings and repeat readings - we only
    ## want to solve the model once for each unique indiv/sample/virus year tested
    titreDat_unique <- titreDat[titreDat$run == 1,]
    titreDat_repeats <- titreDat[titreDat$run != 1,]
    tmp <- row.match(titreDat_repeats[,c("individual","samples","virus")], titreDat_unique[,c("individual","samples","virus")])
    titreDat_repeats$index <- tmp

    measured_strain_indices <- match(titreDat_unique$virus, antigenicMap$inf_years) - 1 ## For each virus tested, what is its index in the antigenic map?
    infection_strain_indices <- match(strain_isolation_times, strain_isolation_times) -1 ## For each virus that circulated, what is its index in the antigenic map?

    ## Get unique measurement sets for each individual at
    ## each sampling time for each repeat
    ## ie. each row of this is a unique blood sample taken
    samples <- unique(titreDat_unique[,c("individual","samples","run")])
    sample_times <- samples$samples ## What were the times that these samples were taken?
    individuals <- samples$individual ## Who are the individuals that these samples correspond to?
    n_indiv <- length(unique(individuals))

    ## Firstly, how many rows in the titre data correspond to each unique individual, sample and titre repeat?
    ## ie. each element of this vector corresponds to one set of titres that need to be predicted
    nrows_per_blood_sample <- NULL
    for(i in 1:nrow(samples)){
        nrows_per_blood_sample <- c(nrows_per_blood_sample, nrow(samples[titreDat_unique$individual == samples[i,"individual"] &
                                                                         titreDat_unique$samples == samples[i,"samples"] &
                                                                         titreDat_unique$run == samples[i,"run"],]))
    }

    ## Which indices in the sampling times vector correspond to each individual?
    ## ie. each contiguous pair of entries in this vector corresponds to the 
    ## first and last entry in the samples matrix that correspond to each individual
    rows_per_indiv_in_samples <- c(0)
    for(individual in unique(individuals)){
        rows_per_indiv_in_samples <- c(rows_per_indiv_in_samples, length(individuals[individuals==individual]))
    }
    rows_per_indiv_in_samples <- cumsum(rows_per_indiv_in_samples)

    ## Which indices in the titre data matrix correspond to each individual?
    ## And, how many rows match each individual?
    nrows_per_individual_in_data <- NULL
    for(individual in unique(individuals)){
        nrows_per_individual_in_data <- c(nrows_per_individual_in_data, nrow(titreDat_unique[titreDat_unique$individual == individual,]))
    }
    cum_nrows_per_individual_in_data <- cumsum(c(0,nrows_per_individual_in_data))
    
    nrows_per_individual_in_data_repeats <- NULL
    for(individual in unique(individuals)){
        nrows_per_individual_in_data_repeats <- c(nrows_per_individual_in_data_repeats, nrow(titreDat_repeats[titreDat_repeats$individual == individual,]))
    }
    cum_nrows_per_individual_in_data_repeats <- cumsum(c(0,nrows_per_individual_in_data_repeats))
    
    titres_unique <- titreDat_unique$titre
    titres_repeats <- titreDat_repeats$titre
    repeat_indices <- titreDat_repeats$index
    repeat_indices_cpp <- repeat_indices - 1

    all_titres <- c(titres_unique, titres_repeats)


    if(!is.null(titreDat$DOB)){
        DOBs <- unique(titreDat[,c("individual","DOB")])[,2]
    } else {
        DOBs <- rep(min(strain_isolation_times), n_indiv)
    }
    if(is.null(ageMask)){
        if(!is.null(titreDat$DOB)){
            ageMask <- create_age_mask(DOBs, strain_isolation_times)
        } else {
            ageMask <- rep(1, n_indiv)
        }
    }
    strainMask <- create_strain_mask(titreDat,strain_isolation_times)
    masks <- data.frame(cbind(ageMask, strainMask))
    if (is.null(n_alive)) {
        n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x)
            nrow(masks[masks$ageMask <=x & masks$strainMask >= x,]))
    }    

    par_names <- parTab$names
    
    if(function_type == 1){
        f <- function(pars, infection_history_mat){
            names(pars) <- par_names
            antigenic_map_long <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma2"])
            y_new <- titre_data_fast(pars, infection_history_mat, strain_isolation_times, infection_strain_indices,
                                     sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
                                     nrows_per_blood_sample, measured_strain_indices, antigenic_map_long,
                                     antigenic_map_short)
            
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
            antigenic_map_long <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma2"])
            ## Now pass to the C++ function
            
            res <- infection_history_proposal_gibbs_fast(pars,
                                                         infection_history_mat,
                                                         probs,
                                                         sampled_indivs,
                                                         nInfs,
                                                         ageMask,
                                                         strainMask,
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
            antigenic_map_long <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma2"])
            y_new <- titre_data_fast(pars, infection_history_mat, strain_isolation_times, infection_strain_indices,
                                     sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
                                     nrows_per_blood_sample, measured_strain_indices, antigenic_map_long,
                                     antigenic_map_short)
        }
    }
    f
}
