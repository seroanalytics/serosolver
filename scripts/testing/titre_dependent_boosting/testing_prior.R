prior <- function(mu,gamma, x){
  res <- sapply(x,function(y) 
    if(mu < 1 + y*gamma | mu > 1 - y*gamma){
      FALSE
      } else {
      TRUE
      })
  res
}

test <- function(gamma){
  return(1/(1-6*gamma))
}

prior1 <- function(mu,gamma,switch_lim=6){
  if(gamma > 0 & mu*gamma >(mu-1)/switch_lim) return(1)
  if(gamma < 0 & mu*(-gamma) > (mu-1)/switch_lim) return(1)
  return(0)
}

prior <- function(mu, gamma, switch_lim){
  if(gamma > (mu-1)/switch_lim) return(1)
  if(gamma < (1-mu)/switch_lim) return(1)
  return(0)
}

titre_dependent1 <- function(mu, gamma, x,switch_lim=6){
  x[x >= switch_lim] <- switch_lim
  return(mu*(1-gamma*x))
}

titre_dependent <- function(mu, gamma, x, switch_lim=6){
  x[x >= switch_lim] <- switch_lim
  return(mu - gamma*x)  
}


switch_lim <- 6
mu_max <- 16
plot(titre_dependent1(2.5,0.1,seq(0,10,by=0.1),switch_lim=switch_lim)~seq(0,10,by=0.1),
     type='l',ylim=c(0,200),col="blue")

values <- expand.grid("mu"=seq(1,mu_max,by=0.1),"gamma"=seq(-1,1,by=0.01))
x <- seq(0,10,by=0.1)
accept <- numeric(nrow(values))
for(i in 1:nrow(values)){
  y <- titre_dependent1(values[i,1],values[i,2],seq(0,10,by=0.1),switch_lim=switch_lim)
  test <- accept[i] <- prior1(values[i,1],values[i,2],switch_lim=switch_lim)
  if(test == 0){
    #print(values[i,])
    lines(y~x,col="blue")
  } else if(test==1){
    lines(y~x,col="red")
  } else {
    lines(y~x,col="purple")
  } 
}
values <- as.data.frame(cbind(values,accept))
abline(h=1)
abline(h=mu_max)
abline(v=switch_lim)

ggplot(values) + geom_tile(aes(x=mu,y=gamma,fill=as.factor(accept)))
