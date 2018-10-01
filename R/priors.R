#' @export
infHistPrior <- function(pars, infHist, ageMask){
    N <- ncol(infHist) - ageMask  + 1
    a <- pars["alpha"]
    b <- pars["beta"]
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


#' @export
inf_mat_prior <- function(infHist, ageMask, alpha, beta1){
    n_alive <- sapply(1:ncol(infHist), function(x) length(ageMask[ageMask <= x]))
    lk <- 0
    for(i in 1:length(n_alive)){
        lk <- lk + log(beta(sum(infHist[,i]) + alpha, n_alive[i]- sum(infHist[,i]) + beta1)/beta(alpha, beta1))
    }
    return(lk)
}
