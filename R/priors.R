#' @export
infHistPrior <- function(pars, infHist, ageMask){
    N <- ncol(infHist) - ageMask  + 1
    a <- pars["alpha"]
    b <- pars["beta"]
    #N <- ncol(infHist)
    #return(rep(-10000,nrow(infHist)))
    priors <- numeric(nrow(infHist))
    for(i in 1:length(priors)){
        priors[i] <- log(dbb_prior(sum(infHist[i,]),N[i],a,b))
    }
    return(priors)
}

#' @export
density_beta_binom<- function(x, N, u, v) {
    (beta(x+u, N-x+v)/beta(u,v))*choose(N,x)
}

#' @export
dbb_prior <- function(x, N, u, v){
    (beta(x+u, N-x+v)/beta(u,v))
}


db <- function(x, a, b){
    x^(a-1)*(1-x^(b-1))/beta(a,b)
}
