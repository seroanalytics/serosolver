#' @export
solve_model_individual <- function(parTab, infectionHistory, sampleTime, antigenicMap, mu_indices=NULL){
    ## Unique strains that an individual could be exposed to
    strains <- unique(antigenicMap$inf_years)
    ## The entry of each strain in the antigenic map table
    strainIndices <- match(strains, strains) - 1
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))

    pars <- parTab$values
    names(pars) <- parTab$names
    ## Work out short and long term boosting cross reactivity - C++ function
    antigenicMapLong <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma1"])
    antigenicMapShort <- create_cross_reactivity_vector(antigenicMapMelted, pars["sigma2"])

    if(is.null(mu_indices)){
        y <- infection_model_indiv(pars, infectionHistory, strains, strainIndices, sampleTime,
                               strainIndices,antigenicMapLong, antigenicMapShort, length(infectionHistory))
    } else {
        mus <- parTab[parTab$identity == 3,"values"]
        print(mus)
        print(mu_indices)
        y <- infection_model_indiv_mus(pars, mus, infectionHistory,
                                       strains, mu_indices, strainIndices, sampleTime,
                                       strainIndices,antigenicMapLong,
                                       antigenicMapShort, length(infectionHistory))
    }
    return(y)

}
