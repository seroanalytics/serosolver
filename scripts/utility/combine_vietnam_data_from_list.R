library(ggplot2)
library(plyr)
setwd("~/Documents/Fluscape/flu-model/sero_model/R_datasets/")
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
write.table(finalDat,"~/Documents/Fluscape/serosolver/data/real/vietnam_data.csv",sep=",",row.names=FALSE)
ggplot(finalDat[finalDat$run == 1 & finalDat$individual %in% 1:5,]) + 
  geom_point(aes(x=virus,y=titre)) + 
  facet_grid(individual~samples)

dat <- read.csv("~/Documents/Fluscape/flu-model/sero_model/datasets/HaNamCohort.csv",stringsAsFactors=FALSE)
dat[dat == "*"] <- NA
dat[,3:ncol(dat)] <- apply(dat[,3:ncol(dat)], 2,function(x) as.character(as.numeric(x)))
dat <- reshape2::melt(dat,id.vars=c("Subject.number","Sample.year"))
dat <- dat[complete.cases(dat),]
dat$value <- as.numeric(dat$value)
dat$value <- log2(dat$value/5)
dat$variable <- as.character(dat$variable)

all_names <- unique(unlist(lapply(test.list, function(x) lapply(x ,colnames))))

                    
