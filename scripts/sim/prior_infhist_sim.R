devtools::load_all("~/Documents/Fluscape/serosolver")
samps <- 1000
n <- 25
indiv <- 50
storage <- matrix(ncol=n+1,nrow=(indiv*samps))
index <- 1
store_indiv <- matrix(ncol=indiv+1,nrow=samps)
store_year <- matrix(ncol=n+1,nrow=samps)
tmp <- matrix(ncol=n,nrow=indiv)
for(i in 1:samps){
  print(index)
  tmp <- sapply(1:n, function(x) rbinom(indiv,1,prob=rbeta(1,1,1)))
  store_indiv[i,] <- c(rowSums(tmp),i)
  store_year[i,] <- c(colSums(tmp),i)
  tmp <- cbind(tmp, i)
  storage[index:(index+indiv-1),] <- tmp
  index <- index + indiv
}

library(reshape2)
library(ggplot2)
store_indiv <- data.frame(store_indiv)
colnames(store_indiv)[ncol(store_indiv)] <- "sampno"
melted_store_indiv <- melt(store_indiv,id.var="sampno")
ggplot(melted_store_indiv) + geom_histogram(aes(x=value),binwidth=1) + facet_wrap(~variable)



store_year <- data.frame(store_year)
colnames(store_year)[ncol(store_year)] <- "sampno"
melted_store_year <- melt(store_year,id.var="sampno")
ggplot(melted_store_year) + geom_histogram(aes(x=value),binwidth=1) + facet_wrap(~variable)
