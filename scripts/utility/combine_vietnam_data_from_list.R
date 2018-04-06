library(ggplot2)
library(plyr)
setwd("~/Documents/Fluscape/serosolver/OLD")
load("HaNam_data.RData")
overall <- NULL
for(indiv in 1:n_part){
  tmp <- NULL
  for(i in seq_along(test_years)){
   tmp[[i]] <- t(test.list[[indiv]][[i]][1:3,])
  } 
  overall[[indiv]] <- as.data.frame(do.call("rbind",tmp))
  row.names(overall[[indiv]]) <- NULL
  overall[[indiv]]$individual <- indiv
  colnames(overall[[indiv]]) <- c("samples","titre","virus","individual")
}
overall <- do.call("rbind",overall)

finalDat <- plyr::ddply(overall,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))
finalDat <- finalDat[order(finalDat$individual, finalDat$run, finalDat$samples, finalDat$virus),]
finalDat <- finalDat[complete.cases(finalDat),]
write.table(finalDat,"~/Documents/Fluscape/serosolver/data/vietnam.csv",sep=",",row.names=FALSE)
ggplot(finalDat[finalDat$run == 1 & finalDat$individual %in% 1:5,]) + 
  geom_point(aes(x=virus,y=titre)) + 
  facet_grid(individual~samples)

