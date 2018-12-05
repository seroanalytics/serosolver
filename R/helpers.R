#' Best pars
#'
#' Given an MCMC chain, returns the set of best fitting parameters (MLE)
#' @param chain the MCMC chain
#' @return a name vector of the best parameters
#' @export
#' @useDynLib serosolver
get_best_pars <- function(chain){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    bestPars <- as.numeric(chain[which.max(chain[,"lnlike"]),2:(ncol(chain)-1)])
    names(bestPars) <- tmpNames
    return(bestPars)
}

#' Index pars
#'
#' Given an MCMC chain, returns the parameters at the specified index
#' @param chain the MCMC chain
#' @param index the index
#' @return a named vector of the best parameters
#' @export
get_index_pars <- function(chain, index){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    pars <- as.numeric(chain[index,2:(ncol(chain)-1)])
    names(pars) <- tmpNames
    return(pars)
}
    
#' PDF - Rich's function to print to device without potential for bad errors
#'
#' Prints to pdf, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @export
to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

#' PNG - Rich's function to print to device without potential for bad errors
#'
#' Prints to png, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @export
to.png <- function(expr, filename, ..., verbose=TRUE) {
    if ( verbose )
        cat(sprintf("Creating %s\n", filename))
    png(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}

#' SVG - Rich's function to print to device without potential for bad errors
#'
#' Prints to SVG, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @export
to.svg <- function(expr, filename, ..., verbose=TRUE) {
    if ( verbose )
        cat(sprintf("Creating %s\n", filename))
    svg(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}


#' Protect function
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
#' @param f the function to be protected
#' @return the protected function
#' @export
#' @useDynLib serosolver
protect <- function(f){
    function(...){
        tryCatch(f(...),error=function(e){
            message("caught error: ", e$message)
            -10000000
        })
    }
}

#' Convert to unit scale
toUnitScale <- function(x, min, max){
    return((x-min)/(max-min))
}

#' Convert from unit scale to original scale
fromUnitScale <- function(x,min,max){
    return(min + (max-min)*x)
}

#' @export
describe_proposals <- function(){
    print("The following text describes the proposal step for updating infection histories.")
    print("Version 1: performs N `flip` proposals at random locations in an individual's infectoin history, switching 1->0 or 0->1. Otherwise, swaps the contents of two random locations")
    print("Version 2: samples from a beta binomial with alpha and beta specified by the parTab input. Only proposes one move at a time")
    print("Version 3: samples from a beta binomial with alpha and beta specified by the parTab input. Proposes nInfs moves at a time for add/remove, or when swapping, swaps locations up to moveSize time steps away")
    print("Version 4: samples from Adam's original proposal")
    print("Version 5: samples directly from the current value of lambda")
    print("Version 6: gibbs sampling of infection histories as in Indian Buffet Process papers")
}

#' @export
logistic_transform <- function(x, maxX){
  return(maxX/(1 + exp(-x)))
}
#' @export
logit_transform <- function(p, maxX){
  return(log(p/(maxX-p)))
}


#' @export
pad_alphas_and_betas <- function(parTab, n_times){
    alpha_row <- parTab[parTab$names == "alpha",]
    beta_row <- parTab[parTab$names == "beta",]

    for(i in 1:(n_times-1)){
        parTab <- rbind(parTab, alpha_row, beta_row)
    }
    parTab    
}



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
setup_titredat_for_posterior_func <- function(titreDat, antigenicMap, ageMask, n_alive){
    strain_isolation_times <- antigenicMap$inf_years
    number_strains <- length(strain_isolation_times)
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))
    
    measured_strain_indices <- match(titreDat$virus, antigenicMap$inf_years) - 1 ## For each virus tested, what is its index in the antigenic map?
    infection_strain_indices <- match(strain_isolation_times, strain_isolation_times) -1 ## For each virus that circulated, what is its index in the antigenic map?

    ## Get unique measurement sets for each individual at
    ## each sampling time for each repeat
    ## ie. each row of this is a unique blood sample taken
    samples <- unique(titreDat[,c("individual","samples","run")])
    sample_times <- samples$samples ## What were the times that these samples were taken?
    individuals <- samples$individual ## Who are the individuals that these samples correspond to?
    n_indiv <- length(unique(individuals))

    ## Firstly, how many rows in the titre data correspond to each unique individual, sample and titre repeat?
    ## ie. each element of this vector corresponds to one set of titres that need to be predicted
    nrows_per_blood_sample <- NULL
    for(i in 1:nrow(samples)){
        nrows_per_blood_sample <- c(nrows_per_blood_sample, nrow(samples[titreDat$individual == samples[i,"individual"] &
                                                                         titreDat$samples == samples[i,"samples"] &
                                                                         titreDat$run == samples[i,"run"],]))
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
        nrows_per_individual_in_data <- c(nrows_per_individual_in_data, nrow(titreDat[titreDat$individual == individual,]))
    }
    cum_nrows_per_individual_in_data <- cumsum(c(0,nrows_per_individual_in_data))

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
    return(list("individuals"=individuals,
                "antigenicMapMelted"=antigenicMapMelted,
                "strain_isolation_times"=strain_isolation_times,
                "infection_strain_indices"=infection_strain_indices,
                "sample_times"=sample_times,
                "rows_per_indiv_in_samples"=rows_per_indiv_in_samples,
                "nrows_per_individual_in_data"=nrows_per_individual_in_data,
                "cum_nrows_per_individual_in_data"=cum_nrows_per_individual_in_data,
                "nrows_per_blood_sample"=nrows_per_blood_sample,
                "measured_strain_indices"=measured_strain_indices,
                "n_indiv"=n_indiv,
                "ageMask"=ageMask,
                "strainMask"=strainMask,
                "n_alive"=n_alive,
                "DOBs"=DOBs))   
    
}
