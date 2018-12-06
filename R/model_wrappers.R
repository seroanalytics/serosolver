#' Solves the re-implementation of the model for a single individual
#' @export
solve_model_individual <- function(par_tab, infectionHistory, sampleTime, antigenic_map, mu_indices=NULL){
    ## Unique strains that an individual could be exposed to
    strains <- unique(antigenic_map$inf_years)
    ## The entry of each strain in the antigenic map table
    strainIndices <- match(strains, strains) - 1
    antigenic_map_melted <- c(outputdmatrix.fromcoord(antigenic_map[,c("x_coord","y_coord")]))

    pars <- par_tab$values
    names(pars) <- par_tab$names
    ## Work out short and long term boosting cross reactivity - C++ function
    antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma1"])
    antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, pars["sigma2"])

    if(is.null(mu_indices)){
        y <- infection_model_indiv(pars, infectionHistory, strains, strainIndices, sampleTime,
                               strainIndices,antigenic_map_long, antigenic_map_short, length(infectionHistory))
    } else {
        mus <- par_tab[par_tab$identity == 3,"values"]
        print(mus)
        print(mu_indices)
        y <- infection_model_indiv_mus(pars, mus, infectionHistory,
                                       strains, mu_indices, strainIndices, sampleTime,
                                       strainIndices,antigenic_map_long,
                                       antigenic_map_short, length(infectionHistory))
    }
    return(y)

}

#' Solves the original model as in the original implementation
#' @export
solve_model_individual_original <- function(theta_star, sample, antigenic.map.in, years, infection_history){ 
    dmatrix = 1- theta_star[["sigma1"]]*outputdmatrix.fromcoord.original(theta_star[["sigma1"]],years,antigenic.map.in,linearD=T)
    dmatrix[dmatrix<0]=0
    dmatrix2 = 1- theta_star[["sigma2"]]*outputdmatrix.fromcoord.original(theta_star[["sigma2"]],years,antigenic.map.in,linearD=T)
    dmatrix2[dmatrix2<0]=0

    ind_yearsA <- 1:length(years)
    titredat <- rep(0, length(years))
    ind_years <- years - 1968 + 1
    testyear_index = sample - 1968 + 1
    testyearI <- which(years == sample)
                                        # Check test data available
                                        # Set up test strains
    test.part=ind_years # index of sample strains data available for
    d.ij=dmatrix[test.part,] # Define cross-immunity matrix 1 for sample strain
    d_vector=melt(t(d.ij))$value #melt is by column
    d.ij2=dmatrix2[test.part,] # Define cross-immunity matrix 2 for sample strain
    d_vector2=melt(t(d.ij2))$value #melt is by column
    expect=c_model_original(length(infection_history),length(years),infection_history, theta_star, d_vector,d_vector2,testyear_index);
    #expect=func1(as.numeric(infection_history),titredat,d_vector,d_vector2,theta_star,testyear_index) # Output expectation
    dat <- data.frame(years=years, y=expect)
    return(dat)
}




#' Needed to solve the original implementation of the model
outputdmatrix.fromcoord.original <- function(thetasigma,inf_years,anti.map.in,linearD=F){ #anti.map.in can be vector or matrix - rows give inf_years, columns give location

  # Check if map is 1D or 2D
  #if(length(anti.map.in)==length(inf_years)){
  if(linearD==F){
    # Exponential decay function
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      (dmatrix=sapply(anti.map.in,function(x){exp(-thetasigma*abs(anti.map.in-x))}))
    }else{ # If spline function defined, calculate directly from input
      (dmatrix=apply(anti.map.in,1,function(x){exp(-thetasigma*sqrt(
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
        ))}))
    }
  }else{
    # Linear decay function
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      (dmatrix=sapply(anti.map.in,function(x){y=abs(anti.map.in-x); y   })) 
      
    }else{
                                        # If spline function defined, calculate directly from input
      (dmatrix=apply(anti.map.in,1,function(x){y=sqrt(
        
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
        
      ); y # 1-1*thetasigma* //  y[y<0]=0; HAVE REMOVED BASE
      }))
    }
  }
}
