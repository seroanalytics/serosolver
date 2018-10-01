############################
## 1. Generate k subsets
############################
## Read in Fluscape data
fluscapeDat <- read.csv("data/real/fluscape_data.csv",stringsAsFactors=FALSE)
fluscapeAges <- read.csv("data/real/fluscape_ages.csv")

## Remove individuals with NA for DOB
na_indiv <- fluscapeAges[which(is.na(fluscapeAges$DOB)),"individual"]
fluscapeDat <- fluscapeDat[-na_indiv,]
fluscapeAges <- fluscapeAges[-na_indiv,]
## Take random subset of individuals
indivs <- unique(fluscapeDat$individual)

titreDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
ages <- fluscapeAges[fluscapeAges$individual %in% indivs,]
titreDat$individual <- match(titreDat$individual, indivs)
ages$individual <- match(ages$individual, indivs)
titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group")]

## Generate sample groups
flds <- createFolds(1:nrow(titreDat), k = 10, list = TRUE, returnTrain = FALSE)

## Isolate titres for the training and test set
titreDat$index <- 1:nrow(titreDat)

for(i in 1:length(flds)){
  ## Get training data - titres less those in the sample group
  testTitres <- titreDat[!(titreDat$index %in% flds[[i]]),]

  ## Make sure that, if we've randomly removed some individuals, we are
  ## only considering individuals in the training set
  test_individuals <- unique(testTitres$individual)
  testAges <- ages[ages$individual %in% test_individuals,]
  
  ## As we cannot infer any titres for individuals with no inferred infection
  ## history, the "full titre set" is just those unique individuals in the
  ## training data
  fullTestTitres <- titreDat[titreDat$individual %in% test_individuals, ]
  
  subsetTitres <- titreDat[titreDat$index %in% flds[[i]],]
  subset_individuals <- unique(subsetTitres$individual)
  subsetAges <- ages[ages$individual %in% subset_individuals,]
  

  write.table(flds[[i]],paste0("~/net/home/serosolver/data_CV/indices_",i,".csv"),sep=",",row.names=FALSE)
  write.table(test_individuals,paste0("~/net/home/serosolver/data_CV/individuals_",i,".csv"),sep=",",row.names=FALSE)
  write.table(testTitres,paste0("~/net/home/serosolver/data_CV/testTitres_",i,".csv"),sep=",",row.names=FALSE)
  write.table(fullTestTitres,paste0("~/net/home/serosolver/data_CV/fullTestTitres_",i,".csv"),sep=",",row.names=FALSE)
  write.table(testAges,paste0("~/net/home/serosolver/data_CV/testAges_",i,".csv"),sep=",",row.names=FALSE)
  write.table(subsetTitres,paste0("~/net/home/serosolver/data_CV/subsetTitres",i,".csv"),sep=",",row.names=FALSE)
  write.table(subsetAges,paste0("~/net/home/serosolver/data_CV/subsetAges",i,".csv"),sep=",",row.names=FALSE)
}

############################
## 2. Generate subsets by individual, 
## so that it is ensured that each individual is missing some data
############################
k <- 10
folds <- vector("list",k)
for(indiv in indivs){
  tmpDat <- titreDat[titreDat$individual ==indiv,]
  tmpFolds <- createFolds(tmpDat$index, k = k, list = TRUE, returnTrain = FALSE)
  for(j in 1:k){
    folds[[j]] <- c(folds[[j]],tmpDat[tmpFolds[[j]],"index"])
  }
}

for(i in 1:length(folds)){
  fold <- folds[[i]]
  ## Get training data - titres less those in the sample group
  testTitres <- titreDat[!(titreDat$index %in% fold),]
  
  ## Make sure that, if we've randomly removed some individuals, we are
  ## only considering individuals in the training set
  test_individuals <- unique(testTitres$individual)
  testAges <- ages[ages$individual %in% test_individuals,]
  
  ## As we cannot infer any titres for individuals with no inferred infection
  ## history, the "full titre set" is just those unique individuals in the
  ## training data
  fullTestTitres <- titreDat[titreDat$individual %in% test_individuals, ]
  
  write.table(fold,paste0("~/net/home/serosolver/data_CV/by_indiv/byindiv_indices_",i,".csv"),sep=",",row.names=FALSE)
  write.table(test_individuals,paste0("~/net/home/serosolver/data_CV/by_indiv/byindiv_individuals_",i,".csv"),sep=",",row.names=FALSE)
  write.table(testTitres,paste0("~/net/home/serosolver/data_CV/by_indiv/byindiv_testTitres_",i,".csv"),sep=",",row.names=FALSE)
  write.table(fullTestTitres,paste0("~/net/home/serosolver/data_CV/by_indiv/byindiv_fullTestTitres_",i,".csv"),sep=",",row.names=FALSE)
  write.table(testAges,paste0("~/net/home/serosolver/data_CV/by_indiv/byindiv_testAges_",i,".csv"),sep=",",row.names=FALSE)
}

############################
## 3. Generate subset without each virus
############################
for(virus in unique(titreDat$virus)){
  ## Get training data - titres less those in the sample group
  testTitres <- titreDat[titreDat$virus != virus,]
  
  ## Make sure that, if we've randomly removed some individuals, we are
  ## only considering individuals in the training set
  test_individuals <- unique(testTitres$individual)
  testAges <- ages[ages$individual %in% test_individuals,]
  
  ## As we cannot infer any titres for individuals with no inferred infection
  ## history, the "full titre set" is just those unique individuals in the
  ## training data
  fullTestTitres <- titreDat[titreDat$individual %in% test_individuals, ]
  
  write.table(fold,paste0("~/net/home/serosolver/data_CV/by_strain/bystrain_indices_",virus,".csv"),sep=",",row.names=FALSE)
  write.table(test_individuals,paste0("~/net/home/serosolver/data_CV/by_strain/bystrain_individuals_",virus,".csv"),sep=",",row.names=FALSE)
  write.table(testTitres,paste0("~/net/home/serosolver/data_CV/by_strain/bystrain_testTitres_",virus,".csv"),sep=",",row.names=FALSE)
  write.table(fullTestTitres,paste0("~/net/home/serosolver/data_CV/by_strain/bystrain_fullTestTitres_",virus,".csv"),sep=",",row.names=FALSE)
  write.table(testAges,paste0("~/net/home/serosolver/data_CV/by_strain/bystrain_testAges_",virus,".csv"),sep=",",row.names=FALSE)
}
