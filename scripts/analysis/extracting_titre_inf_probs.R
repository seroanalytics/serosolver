setwd("/media/james/F40474B404747B84/Storage/serosolver/outputs_26032018/real_lambda_1")
devtools::load_all("~/Documents/Fluscape/serosolver")
library(data.table)
library(ggplot2)
chain <- read.csv("real_lambda_1_1_chain.csv")
infChain <- data.table::fread("real_lambda_1_1_infectionHistories.csv")
fluscapeDat <- read.csv("~/Documents/Fluscape/serosolver/data/real/fluscape_data.csv",stringsAsFactors=FALSE)
fluscapeAges <- read.csv("~/Documents/Fluscape/serosolver/data/real/fluscape_ages.csv")
## Remove individuals with NA for DOB
na_indiv <- fluscapeAges[which(is.na(fluscapeAges$DOB)),"individual"]
fluscapeDat <- fluscapeDat[-na_indiv,]
fluscapeAges <- fluscapeAges[-na_indiv,]
## Take random subset of individuals
indivs <- unique(fluscapeAges$individual)
titreDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
ages <- fluscapeAges[fluscapeAges$individual %in% indivs,]
titreDat$individual <- match(titreDat$individual, indivs)
ages$individual <- match(ages$individual, indivs)
titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group")]

indivs <- 1:1000

titreDat1 <- expand.grid(individual=indivs, samples=1967:2014,virus=1968:2015,titre=0,run=1,group=1)
titreDat1 <- titreDat1[order(titreDat1$individual, titreDat1$samples,titreDat1$virus),]
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",stringsAsFactors=FALSE)
fit_dat <- read.csv("~/Documents/Fluscape/serosolver/data/antigenic_maps/antigenicMap_vietnam.csv")

#Rprof(tmp<-tempfile())
titres <- get_titre_predictions(chain, infChain, titreDat1, indivs, fit_dat, ages, parTab, 100, TRUE, NULL, TRUE)
#Rprof(NULL)
#library(proftools)
#plotProfileCallGraph(readProfileData(tmp),score = "total")

titres1 <- cbind(titreDat1,titres[[4]])
titres1 <- titres1[,c("individual","samples","virus",colnames(titres[[4]]))]
titres1 <- melt(titres1, id.vars=c("individual","samples","virus"))
colnames(titres1) <- c("individual","samples","virus","sampno","titre")
titres1$sampno <- as.numeric(as.character(titres1$sampno))
infChain1 <- infChain
colnames(infChain1) <- c("individual","virus","infected","sampno")
infChain1 <- infChain1[infChain1$sampno %in% unique(titres1$sampno),]
infChain1 <- infChain1[infChain1$individual %in% indivs,]
viruses <- 1968:2015
infChain1$virus <- viruses[infChain1$virus]
infChain1$samples <- infChain1$virus - 1
titres2 <- titres1[,c("individual","virus","sampno","titre","samples")]

## Which year's sample are we interested in for which virus if there's an infection?
pre_inf_years <- data.frame("samples"=1967:2014,"virus"=1968:2015)
#not_pre_inf_years <- expand.grid("samples"=1967:2014,"virus"=1968:2015)
#not_pre_inf_years <- dplyr::setdiff(not_pre_inf_years,pre_inf_years)
titres3 <- merge(data.table(titres2), data.table(pre_inf_years))
## Have only titres for years preceding infections
titres_pre_infection_years <- merge(titres3, infChain1, by=c("individual","virus","sampno","samples"))

## Now we need titres for years that didn't precede infections
titres_pre_infection_years_less <- titres_pre_infection_years[,colnames(titres3), with=FALSE]
titres_not_pre_infection_years <- dplyr::setdiff(titres3,titres_pre_infection_years_less)
titres_not_pre_infection_years$infected <- 0
all_titres <- rbind(titres_pre_infection_years, titres_not_pre_infection_years[,colnames(titres_pre_infection_years),with=FALSE])
all_titres <- data.table(all_titres)
rows <- sample(1:nrow(all_titres),1000)
thinned <- all_titres[rows,]
library(popbio)
logi.hist.plot(all_titres$titre,all_titres$infected,boxp=FALSE,type="hist",col="gray")

model <- glm(infected~poly(titre,2),family=binomial(link='logit'),data=all_titres)
y <- predict(model, x=data.frame(titre=seq(0,10,by=0.1)),type="resp")
plot(y)

summary(model)

colnames(infChain1)[2] <- "samples"
#all_titres$samples <- all_titres$samples + 1
indiv <- 3
ggplot(titres2[titres2$individual == indiv &
                    titres2$sampno == unique(titres2$sampno)[1],]) + 
  geom_line(aes(x=virus,y=titre))+ 
  geom_point(data=titreDat[titreDat$individual == indiv,],aes(x=virus,y=titre)) +
  geom_vline(data=infChain1[infChain1$individual == indiv & 
                              infChain1$sampno == unique(titres2$sampno)[1],],aes(xintercept=samples)) +
  facet_wrap(~samples)



wow <- NULL
for(indiv in unique(infChain1$individual)){
  print(indiv)
  tmp <- infChain1[infChain1$individual == indiv,]
  for(i in unique(tmp$sampno)){
    tmp1 <- subset(tmp, sampno == i)
    distances <-diff(tmp1$samples)
    wow <- c(wow, distances)
  }
}
hist(wow)
