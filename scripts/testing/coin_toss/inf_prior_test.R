
n <- 2000
Zs <- sample(c(0,1),n,prob=c(0.9,0.1),replace=TRUE)
indiv_prob <- sum(log(dbb(Zs, 1, alpha, beta)))
pars <- expand.grid(seq(1,100,by=1),seq(1,100,by=1))
#pars <- matrix(runif(2*1000,0,100),nrow=1000)
colnames(pars) <- c("alpha","beta")
liks <- NULL
liks1 <- NULL
probs <- runif(100*100,0,1)
for(i in 1:nrow(pars)){
  liks[i] <- sum(log(dbb(Zs, 1, pars[i,1],pars[i,2])))  
  liks1[i] <- sum(log((probs[i]^Zs)*(1-probs[i])^(1-Zs)))
}
pars1 <- tibble(alpha=pars[,1],beta=pars[,2],liks=liks)
pars1 <- pars1[order(pars1$beta, pars1$alpha),]
ggplot(pars1,aes(x=alpha,y=beta)) + geom_raster(aes(fill=exp(liks)),interpolate=TRUE) + 
  scale_fill_gradient(low="blue",high="red")+
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0))

alpha=pars[which.max(liks),1]
beta=pars[which.max(liks),2]

par(mfrow=c(1,2))
plot(density(rbb(10000,n,alpha,beta)/n),xlim=c(0,1))
plot(exp(liks1)~probs)



dbernoulli <- function(psi, Z){
  return(psi^Z * (1-psi)^(1-Z))
}

liks <- NULL
probs <- NULL
Zs <- NULL
for(i in 1:10000){
  prob <- rbeta(1,1,1)
  #prob <- 0.5
  probs[i] <- prob
  Z <- rbinom(1,1,prob)
  Zs[i] <- Z
  #Z <- 0
  liks[i] <- dbernoulli(prob, Z)*dbeta(prob,1,1)
}
plot(liks~probs,ylim=c(0,1))

