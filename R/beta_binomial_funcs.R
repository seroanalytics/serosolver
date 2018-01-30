#' @export
dbb <- function(x, N, u, v) {
  beta(x+u, N-x+v)/beta(u,v)*choose(N,x)
}
#' @export
pbb <- function(q, N, u, v) {
  sapply(q, function(xx) sum(dbb(0:xx, N, u, v)))
}

#' @export
qbb <- function(p, N, u, v) {
  pp <- cumsum(dbb(0:N, N, u, v))
  sapply(p, function(x) sum(pp < x))
}

#' @export
rbb <- function(n, N, u, v) {
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}

#' @export
bb_mean <- function(n, alpha, beta){
  return(n*alpha/(alpha+beta))
}

#' @export
bb_var <- function(n, alpha, beta){
  top <- n*alpha*beta*(alpha+beta+n)
  bot <- ((alpha+beta)^2) *(alpha+beta+1)
  return(top/bot)
}

#' @export
hist_rbb <- function(n, mean, var){
  pars <- find_a_b(n,mean,var)
  a <- pars["a"]
  b <- pars["b"]
  hist(rbb(10000,n,a,b),breaks=seq(-1,n,by=1))
}

#' @export
find_a_b <- function(n, mean, var){
  y <- mean
  z <- var
  x <- (n-y)/y
  top <- z*(1+x)^2 - (n^2)*x
  bot <- n*x*(1+x) - z*((1+x)^3)
  a2 <- top/bot
  b2 <- x*a2
  return(c("a"=a2,"b"=b2))
}

#' @export
find_bb_2 <- function(n1, a1, b1, n2){
  y <- bb_mean(n1,a1,b1)
  z <- bb_var(n1,a1,b1)
  x <- (n2-y)/y
  top <- z*(1+x)^2 - (n2^2)*x
  bot <- n2*x*(1+x) - z*((1+x)^3)
  a2 <- top/bot
  b2 <- x*a2
  return(c("a2"=a2,"b2"=b2))
  
}
