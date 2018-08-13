#' Convert to unit scale
toUnitScale <- function(x, min, max){
    return((x-min)/(max-min))
}

#' Convert from unit scale to original scale
fromUnitScale <- function(x,min,max){
    return(min + (max-min)*x)
}

likelihood <- function(pars, coin_results, dat){
  #return(-1000000)
  y <- coin_toss_function(pars, coin_results)
  return(sum(dnorm(dat, mean=y, sd=pars[3],log=TRUE)))
}

likelihood_group <- function(pars, coin_results, dat, samps=seq(1,20,by=4)){
  #return(rep(-100000,nrow(coin_results)))
  y <- coin_toss_group(pars, coin_results)
  liks <- numeric(nrow(dat))
  for(i in 1:nrow(y)){
    liks[i] <- sum(dnorm(dat[i,samps],mean=y[i,samps],sd=pars[3],log=TRUE))
  }
  return(liks)
}

inf_mat_prior <- function(infHist, alpha, beta){
  N <- nrow(infHist)
  apply(infHist, 2, function(x) log(beta(sum(x) + alpha, N - sum(x) + beta)/beta(alpha, beta)))
}


prior_gibbs <- function(pars){
  a <- dunif(pars[1],0,10,log=TRUE)
  b <- dunif(pars[2],0,1,log=TRUE)
  c <- dunif(pars[3],0,5,log=TRUE)
  return(a+b+c)
}

prior_probs <- function(pars, probs){
  a <- dunif(pars[1],0,10,log=TRUE)
  b <- dunif(pars[2],0,1,log=TRUE)
  c <- dunif(pars[3],0,5,log=TRUE)
  #a <- dnorm(pars[1], 0, 1000,log=TRUE)
  #b <- dnorm(pars[2],0,1000,log=TRUE)
  #c <- dnorm(pars[3],0,1000,log=TRUE)
  d <- sum(dbeta(probs, 1,1,log=TRUE))
  return(a+b+c+d)
}


prior_group <- function(coin_probs, coin_results){
  prob <- numeric(nrow(coin_probs))
  for(i in 1:nrow(coin_probs)){
    prob[i] <- sum(log(coin_probs[i,]^coin_results[i,] * (1-coin_probs[i,])^(1-coin_results[i,])))
  }
  return(prob)
}


hyper_prior_group <- function(probs, coin_results){
  prob <- apply(coin_results, 1, function(x) sum(log(probs^x * (1-probs)^(1-x))))
  return(prob)
}


posterior_group <- function(pars, probs, coin_results, dat, samps){
  return(sum(likelihood_group(pars, coin_results, dat, samps)) + sum(hyper_prior_group(probs, coin_results)) + prior(pars, probs))
}
