create_proposal_bundle <- function(parTab, block, proposal_ver, posterior_ver, mcmcPars){
    if(proposal_ver == 1){
        bundle <- create_proposal_bundle_1(parTab, block, mcmcPars)
    } else if(proposal_ver == 2){
        bundle <- create_proposal_bundle_2(parTab, block, mcmcPars)
    } else if(proposal_ver == 3){
        bundle <- create_proposal_bundle_3(parTab, block, mcmcPars)
    } else {
        bundle <- create_proposal_bundle_4(parTab, block, mcmcPars)
    }    
    return(bundle)
}

create_proposal_bundle_1 <- function(parTab, block){
    pars <- new_pars <- parTab$values
    par_indices <- which(parTab$block == block & parTab$fixed == 0)
    par_names <- parTab[parTab$block == block, "names"]
    
    new_probab <- NULL
    new_probabs <- NULL
    steps <- parTab[parTab$block == block, "steps"]

    j <- 1
    n_par <- length(par_indices)

    tempiter <- tempaccepted <- reset <- integer(n_par)
    popt <- mcmcPars["popt"]
    pcur <- integer(n_par)
    lower_bounds <- parTab$lower_bound
    upper_bounds <- parTab$upper_bound

    bundle <- list(pars, new_pars, par_names, par_indices,
                   new_probab, new_probabs,
                   steps,
                   j, n_par,
                   tempiter, tempaccepted, reset,
                   popt, pcur,
                   lower_bounds, upper_bounds)
    return(bundle)
}


proposal_func_1 <- function(current_pars, infectionHistories, bundle){
    ## Find which parameter we're updating
    j <- bundle[["j"]]
    index <- bundle[["par_indices"]][j]

    ## Update this parameter
    bundle[["new_pars"]][index] <- univ_proposal(current_pars, bundle[["lower_bounds"]],
                                                 bundle[["upper_bounds"]],bundle[["steps"]], index)
    
    ## Update iterations for this parameter
    bundle[["tempiter"]][j] <- bundle[["tempiter"]][j] + 1
    j <- j + 1
    if(j > bundle[["n_par"]]) j <- 1
    bundle[["j"]] <- j
    return(bundle)    
}

acceptance_func_1 <- function(bundle, infectionHistories, probab, probabs){
    pars <- bundle[["new_pars"]]
    
    
}
