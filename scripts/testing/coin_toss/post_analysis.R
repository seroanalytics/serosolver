setwd("~/net/home/coinflip/outputs1")
library(coda)
files <- list.files()

for(dir in files){
  print(dir)
  setwd("~/net/home/coinflip/outputs1")
  setwd(dir)
  all_files <- list.files()
  print(length(all_files))
  chains <- NULL
  if(length(all_files) > 2){
  for(file in all_files){
    chain <- as.mcmc(data.table::fread(file,data.table=FALSE))
    chains[[file]] <- chain[,2:ncol(chain)]
  }
  chains <- mcmc.list(chains)
  ##tryCatch({
  #  print(gelman.diag(chains)$mpsrf)
  #}, error=function(e){})
  tryCatch(print(mean(effectiveSize(chains))))
  pdf(paste0(dir,"_plot.pdf"))
  plot(chains)
  dev.off()
  setwd("~/net/home/coinflip/outputs1")
 
  }
}