univ_proposal <- function(values, lower_bounds, upper_bounds,steps, index){
  rtn <- values
  #rtn[index] <- rnorm(1,values[index],steps[index])
  #return(rtn)
  mn <- lower_bounds[index]
  mx <- upper_bounds[index]
  
  rtn <- values
  
  x <- toUnitScale(values[index],mn,mx)
  
  stp <- steps[index]
  
  rv <- runif(1)
  rv <- (rv-0.5)*stp
  x <- x + rv
  
  ## Bouncing boundary condition
  #if (x < 0) x <- -x
  #if (x > 1) x <- 2-x
  
  ## Cyclical boundary conditions
  if (x < 0) x <- 1 + x	
  if (x > 1) x <- x - 1
  
  if(x < 0 | x > 1){
    print("Stepped outside of unit scale. Something went wrong...")
  }
  
  rtn[index] <- fromUnitScale(x,mn,mx)
  rtn
}

univ_proposal_theta <- function(pars, steps, index){
  proposed <- pars
  proposed[index] <- proposed[index] + rnorm(1, mean=0, sd=exp(steps[index])^2)
  return(proposed)
}

scaletuning1 <- function(step, popt,pcur){
  if(pcur ==1) pcur <- 0.99
  if(pcur == 0) pcur <- 0.01
  step = (step*qnorm(popt/2))/qnorm(pcur/2)
  if(step > 1) step <- 1
  step <- max(0.00001, step)
  return(step)
}

proposal_theta <- function(pars, fixed, covMat, covMat0, beta=0.05){
  proposed <- pars
  proposed[fixed] <- (1-beta)*MASS::mvrnorm(n=1,mu=proposed[fixed], Sigma=5.6644*covMat/length(fixed)) + 
    beta*MASS::mvrnorm(n=1,mu=proposed[fixed], Sigma=0.001*covMat0/length(fixed))
  return(proposed)
}

# Given the number of infections in each year, generates a proposal for
# the probs of infection
foi_proposal <- function(fois, infMat, is, alpha=1, beta=1, n){
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

## If using individual coin probs per individual per year
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
