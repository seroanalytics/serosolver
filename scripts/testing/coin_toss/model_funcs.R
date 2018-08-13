coin_toss_group <- function(pars, coin_results){
  indivs <- nrow(coin_results)
  y <- t(apply(coin_results, 1, function(x) coin_toss_function(pars, x)))
  return(y)
}

coin_toss_function <- function(pars, coin_results){
  mu <- pars[1]
  sigma <- pars[2]
  res <- numeric(length(coin_results))
  for(toss in seq_along(coin_results)){
    ## Boost to actual coin
    res[toss] <- res[toss] + mu*coin_results[toss]
    ## Boost to adjacent coins
    adjacent_toss <- c(toss + 1, toss - 1)
    adjacent_toss <- adjacent_toss[adjacent_toss > 0 & adjacent_toss <= length(res)]
    res[adjacent_toss] <- res[adjacent_toss] + sigma*mu*coin_results[toss]
  }
  return(res)
}

measurement_error_group <- function(pars, dat){
  error <- pars[3]
  dat <- t(apply(dat, 1, function(x) rnorm(length(x), x, error)))
  return(dat)
}

extract_number_infections_from_chain <- function(x, n, by.year=FALSE){
  colnames(x) <- c("sampno",seq(1,n,by=1),"indiv")
  if(!by.year) sums <- ddply(as.data.frame(x),~sampno, function(x) sum(x[,2:(ncol(x)-1)]))
  else{
    x <- reshape2::melt(as.data.frame(x),id.vars=c("sampno","indiv"))
    x[,"value"] <- as.numeric(x[,"value"])
    sums <- ddply(as.data.frame(x),.(sampno,variable), function(x) sum(x$value))
    sums <- dcast(sums, sampno~variable)  
  }
  return(sums)
}
