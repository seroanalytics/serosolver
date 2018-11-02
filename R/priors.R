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

fit_beta_prior <- function(chain_samples,parName="",error_tol=999999999,try_attempts=10){
  data <- density(chain_samples)
  data <- data.frame(x=data$x, y=data$y)
  x <- data$x
  
  f <- function(pars){
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dbeta(x, shape1, shape2)
    return(sum((out - data$y)^2))
  }
  res <- list(value=error_tol+1)
  try_no <- 1
  while(res$value > error_tol & try_no < try_attempts){
    message(cat("Try no: ", try_no))
    res <- optim(runif(2,1,100),f,method="Nelder-Mead",control=list(abstol=1e-8,reltol=1e-8))    
    try_no <- try_no + 1
  }
  if(try_no == try_attempts){
    message("Fitting did not work - does this distribution make sense?")
    return(list(par=c(NA,NA)))
  }
  plot(dbeta(x,res$par[1],res$par[2])~x,type='l',col="blue",
       main=paste0("Beta dist fit to posterior of ",parName),
       xlab="Value",
       ylab="Density")
  legend(x=x[1],y=max(data$y),box.lty=0, lty=1,col=c("blue","red"),legend=c("Model fit","Posterior density"))
  lines(data,col="red")
  return(res)
}

fit_normal_prior <- function(chain_samples, parName="",error_tol=999999999,try_attempts=10){
  data <- density(chain_samples)
  data <- data.frame(x=data$x, y=data$y)
  x <- data$x
  
  f <- function(pars){
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dnorm(x, mean=shape1, sd=shape2)
    return(sum((out - data$y)^2))
  }
  start_mean <- mean(chain_samples)
  start_sd <- sd(chain_samples)
  res <- list(value=error_tol+1)
  try_no <- 1
  while(res$value > error_tol & try_no < try_attempts){
    message(cat("Try no: ", try_no))
    res <- optim(c(start_mean, start_sd),f,method="Nelder-Mead",control=list(abstol=1e-8,reltol=1e-8))
    try_no <- try_no + 1
  }
  if(try_no == try_attempts){
    message("Fitting did not work - does this distribution make sense?")
    return(list(par=c(NA,NA)))
  }
  plot(dnorm(x,res$par[1],res$par[2])~x,type='l',col="blue",
       main=paste0("Normal dist fit to posterior of ",parName),
       xlab="Value",
       ylab="Density")
  legend(x=x[1],y=max(data$y),box.lty=0, lty=1,col=c("blue","red"),legend=c("Model fit","Posterior density"))
  lines(data,col="red")
  return(res)
}

#' @export
find_beta_parameters <- function(mean, var){
    #if(mean < 0 | mean > 1) stop("Mean is outside of bounds (0, 1)")
    #if(var < 0 | var > 0.25^2) stop("Var is outside of bounds (0, 0.25^2)")
    alpha <- ((1 - mean)/(var) - (1/mean))*mean^2
    beta <- alpha*(1/mean - 1)
    x <- seq(0,1,by=0.001)
    y <- dbeta(x, alpha,beta)
    plot(y~x)
    return(list(alpha=alpha,beta=beta))
}
