#' @export
infHistPrior <- function(pars, infHist){
    a <- pars["alpha"]
    b <- pars["beta"]
    N <- ncol(infHist)
    priors <- log(apply(infHist, 1, function(x) dbb_prior(sum(x), N, a, b))  )
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
