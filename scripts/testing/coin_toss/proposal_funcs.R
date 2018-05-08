proposal_theta <- function(pars, fixed, covMat, covMat0, beta=0.05){
  proposed <- pars
  proposed[fixed] <- (1-beta)*MASS::mvrnorm(n=1,mu=proposed[fixed], Sigma=5.6644*covMat/length(fixed)) + 
    beta*MASS::mvrnorm(n=1,mu=proposed[fixed], Sigma=0.001*covMat0/length(fixed))
  return(proposed)
}

coin_proposal <- function(coin_results, k=1){
  #proposed <- sample(c(0,1),length(coin_results),replace=TRUE)
  locs <- sample(length(coin_results),k)
  proposed <- coin_results
  proposed[locs] <- !coin_results[locs]
  return(proposed)
}

coin_proposal_symmetric_group <- function(coin_results, k=1, indivs=1){
  proposed <- coin_results
  for(j in indivs){
    if(runif(1) < 0.5){
      locs <- sample(length(coin_results[j,]), k)
      proposed[j,locs] <- !coin_results[j,locs]
    } else {
      locs <- sample(length(coin_results[j,]),2)
      loc1 <- locs[1]
      loc2 <- locs[2]
      tmp <- proposed[j,loc1]
      proposed[j,loc1] <- proposed[j,loc2]
      proposed[j,loc2] <- tmp
    }
  }
  return(proposed)
}

foi_proposal <- function(fois, infMat, is, alpha, beta, n){
  proposed <- fois
  if(length(is) > 1){
    infs <- colSums(infMat[,is])
    proposed[is] <- rbeta(length(is), alpha+infs, beta + n[is]-infs)
  } else {
    infs <- sum(infMat[,years])
    proposed[is] <- rbeta(1, alpha+infs, beta + n[is]-infs)
  }
  forward <- sum(dbeta(proposed[is], alpha,beta,log=TRUE))
  back <- sum(dbeta(fois[is],alpha,beta,log=TRUE))
  ratio <- back - forward
  return(proposed)
  return(list(proposed,ratio))
}

coin_prob_proposal <- function(coin_probs, k=1, indivs=1){
  proposed_coin_prob <- coin_probs
  for(j in indivs){
    if(runif(1) < 1){
      locs <- sample(length(coin_probs[j,]), k)
      proposed_coin_prob[j,locs] <- runif(k)
    } else {
      locs <- sample(length(coin_results[j,]),2)
      loc1 <- locs[1]
      loc2 <- locs[2]
      tmp <- proposed_coin_prob[j,loc1]
      proposed_coin_prob[j,loc1] <- proposed_coin_prob[j,loc2]
      proposed_coin_prob[j,loc2] <- tmp
    }
  }
  return(proposed_coin_prob)
}
