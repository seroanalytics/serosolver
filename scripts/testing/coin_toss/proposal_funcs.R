# Just flips the entries of k locations
coin_proposal_simple <- function(coin_results, k=1){
  locs <- sample(length(coin_results),k)
  proposed <- coin_results
  proposed[locs] <- !coin_results[locs]
  return(proposed)
}

# Most complicated sampler - uses gibbs sampling for each feature to propose move with correct probability
coin_proposal_gibbs <- function(pars, coin_results, dat, 
                                alpha=1, beta = 1, 
                                indiv_propn=1,year_propn=1,
                                swap_distance=1,
                                swapPropn=0.5){
  proposed <- coin_results
  years <- ncol(coin_results)
  
  ## For each individual
  indivs <- sample(1:nrow(coin_results), floor(indiv_propn*nrow(coin_results)))
  iter <- accepted <- 0
  
  for(indiv in indivs){
    ## Get this individual's infection history
    x <- proposed[indiv,]
    if(runif(1) > swapPropn){
      
      ## Swap with adjacent year
      loc1 <- sample(1:length(x),1)
      loc2 <- loc1 + sample(c(swap_distance,-swap_distance),1)
      if(loc2 < 1) loc2 <- loc2 + length(x)
      if(loc2 > length(x)) loc2 <- loc2 %% length(x)
      
      if(x[loc1] != x[loc2]){
        tmp <- x[loc1]
        x[loc1] <- x[loc2]
        x[loc2] <- tmp
        old_prob <- likelihood(pars, proposed[indiv,],dat[indiv,])
        new_prob <- likelihood(pars, x,dat[indiv,])
        log_prob <- min(new_prob - old_prob,0)
        if(is.finite(log_prob) & log(runif(1)) < log_prob){
          proposed[indiv,] <- x 
        }
      }
    } else {
      ## For each year
      years <- sample(1:ncol(coin_results),floor(year_propn*ncol(coin_results)))
      for(year in years){
        iter <- iter + 1
        ## Get number of other infections in this year
        m <- sum(proposed[-indiv,year])
        number <- runif(1)
        #prob <- (m + alpha/ncol(coin_results))/(length(x[-indiv]) + ncol(coin_results))
        ## Make the infection state for this individual 1 or 0 depending on prior and other individuals, gibbs
        prob <- (m + alpha)/(length(proposed[-indiv,year]) + alpha+ beta)
        if(number < prob){
          x[year] <- 1
        } else {
          x[year] <- 0
        }
        ## If proposing a change, need to check likelihood ratio
        if(x[year] != proposed[indiv,year]){
          old_prob <- likelihood(pars, proposed[indiv,],dat[indiv,])
          new_prob <- likelihood(pars, x,dat[indiv,])
          log_prob <- min(new_prob - old_prob,0)
          if(is.finite(log_prob) & log(runif(1)) < log_prob){
            proposed[indiv,year] <- x[year] 
            accepted <- accepted + 1
          }
        }
      }
    }
  }
  #message(paste0("Acceptance rate: ",accepted/iter))
  return(proposed)
}

coin_proposal_by_year <- function(coin_results, n_years=1, n_indiv=1, swapPropn=0.5){
  ver <- runif(1)
  proposed <- coin_results
  if(ver < swapPropn){
    #print("swap years")
    locs <- sample(ncol(proposed),2)
    tmp <- proposed[,locs[1]]
    proposed[,locs[1]] <- proposed[,locs[2]]
    proposed[,locs[2]] <- tmp
  } else {
    #print("swap individuals")
    indivs <- sample(1:nrow(proposed),n_indiv)
    locs <- sample(ncol(proposed),n_years)
    proposed[indivs,locs] <- !proposed[indivs,locs]
  }
  return(proposed)
}

coin_proposal_add_or_subtract <- function(coin_results, n_years=1, n_indivs=1, swapPropn){
  ver <- runif(1)
  proposed <- coin_results
  ## Add or remove infections
  if(ver < 1){
    ver1 <- runif(1)
    if(ver1 < swapPropn){
      years <- sample(1:ncol(proposed),n_years)
      proposed[,years] <- sapply(years, function(x){
        locs <- sample(1:nrow(proposed), n_indivs)
        proposed[locs,x] <- 1
        proposed[, x]
      })
    } else {
      years <- sample(1:ncol(proposed),n_years)
      proposed[,years] <- sapply(years, function(x){
        locs <- sample(1:nrow(proposed), n_indivs)
        proposed[locs,x] <- 0
        proposed[, x]
      })
    }
  }
  proposed
}


# Choose n indivs and change contents for k locations, or swap two years
coin_proposal_symmetric_group <- function(coin_results, k=1, indivs=1, swapPropn=0.5){
  proposed <- coin_results
  for(j in indivs){
    if(runif(1) < swapPropn){
      locs <- sample(length(coin_results[j,]), k)
      proposed[j,locs] <- !coin_results[j,locs]
      #probs <- rbeta(k, 1,1)
      #proposed[j,locs] <- rbinom(k,1,probs)
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

## For each individual and year (ny/ni), resample an infection/clear based on infection prob
coin_proposal_probs <- function(coin_results, coin_probs, ny=1, ni=1){
  proposed <- coin_results
  indivs <- sample(1:nrow(coin_results), ni)
  for(indiv in indivs){
    years <- sample(1:ncol(coin_results),ny)
    for(year in years){
      if(runif(1) < coin_probs[year]){
        proposed[indiv, year] <- 1
      } else {
        proposed[indiv, year] <- 0
      }
    }
  }
  return(proposed)
}
