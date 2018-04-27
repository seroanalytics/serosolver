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
}

#' @export
logistic_transform <- function(x, maxX){
  return(maxX/(1 + exp(-x)))
}
#' @export
logit_transform <- function(p, maxX){
  return(log(p/(maxX-p)))
}
