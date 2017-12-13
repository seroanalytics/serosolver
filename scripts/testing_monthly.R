library(ggplot2)
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map.csv")

strains <- unique(antigenicMap$inf_years)
numberStrains <- length(strains)
strainIsolationTimes <- strains
#strainIsolationTimes <- sample(strains,15)
#strainIsolationTimes <- strainIsolationTimes[order(strainIsolationTimes)]
virusIndices <- match(strainIsolationTimes,strains)-1

antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))

pars <- parTab$values
mynames <- parTab$names
names(pars) <- parTab$names

antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
antigenicMapShort[antigenicMapShort < 0] <- 0

infectionHistory <- sample(c(0,0,0,1),numberStrains,replace=T)
infectionHistory <- rep(0,numberStrains)
infectionHistory[which(strains==1968*12)] <- 1
infectionHistory[which(strains==1990*12)] <- 1
infectionHistory[which(strains==2009*12)] <- 1


sampleTimes <- c(1980,1990,2000,2010)*12

nsamps <- 30
strainIsolationTimesA <- sample(strainIsolationTimes,nsamps)
strainIsolationTimesA <- strainIsolationTimesA[order(strainIsolationTimesA)]

strainIsolationTimes1 <- rep(strainIsolationTimesA, length(sampleTimes))
dataIndices <- rep(length(strainIsolationTimesA),length(sampleTimes))

theta <- pars
strainIsolationTimes <- strainIsolationTimes1

library(microbenchmark)
a <- infectionHistory[which(infectionHistory == 1)]
b <- strains[which(infectionHistory == 1)]
c <- match(strains[which(infectionHistory==1)],strains)-1
d <- match(strainIsolationTimesA,strains)-1
microbenchmark(infection_model_indiv(pars, 
                          a, 
                          b,
                           c,
                           24109, 
                           d, 
                           antigenicMapLong, 
                           antigenicMapShort, 
                           length(strains)))


microbenchmark(infection_model_indiv_OLD(pars, 
                                     infectionHistory, 
                                     24109,
                                     strainIsolationTimesA,
                                     d,
                                     antigenicMapLong, 
                                     antigenicMapShort, 
                                     length(strains)))


strainIsolationTimes1 <- rep(strains, length(sampleTimes))
dataIndices <- rep(length(strains),length(sampleTimes))
microbenchmark(simulate_individual(pars, infectionHistory,sampleTimes,dataIndices,strainIsolationTimes1,match(strainIsolationTimes1, strains)-1,antigenicMapLong,antigenicMapShort,strains))

y <- as.data.frame(y)
colnames(y) <- c("sample","virus","titre")
infDat <- data.frame(x=strains[which(infectionHistory==1)])
sampDat <- data.frame(x=sampleTimes)
ggplot(y) + geom_point(aes(x=virus,y=titre)) + 
            geom_vline(data=infDat,aes(xintercept=x)) + 
  geom_vline(data=sampDat,aes(xintercept=x),col="red")+
  facet_wrap(~sample,ncol=1)+
  theme_bw()

circulationTimeIndices = match(strains, strains)-1
measurementMapIndices = match(strainIsolationTimes1,strains)-1
y <- titre_data_individual(pars, infectionHistory, strains,circulationTimeIndices , sampleTimes, dataIndices, measurementMapIndices,strainIsolationTimes1,
                           antigenicMapLong,antigenicMapShort,length(strains))
samples <- rep(sampleTimes, each=length(strains))

dat <- data.frame(sample=samples,virus=strainIsolationTimes1,titre=y)

dat1 <- simulate_individual(pars, infectionHistory,sampleTimes,dataIndices,strainIsolationTimes1,match(strainIsolationTimes1, strains)-1,antigenicMapLong,antigenicMapShort,strains)
dat1 <- as.data.frame(dat1)
colnames(dat1) <- c("sample","virus","titre")


infectionHistory <- rep(0,numberStrains)
infectionHistory[which(strains==1968*12)] <- 1
infectionHistory[which(strains==1990*12)] <- 1
infectionHistory[which(strains==2009*12)] <- 1

sampleTimes <- seq(1970,2010,by=10)*12
strains <- unique(antigenicMap$inf_years)
strainIsolationTimesA <- sample(strains,length(strains))
strainIsolationTimesA <- strainIsolationTimesA[order(strainIsolationTimesA)]

dataIndices <- rep(length(strainIsolationTimesA),length(sampleTimes))
strainIsolationTimes1 <- as.numeric(rep(strainIsolationTimesA, length(sampleTimes)))
infectionHistory <- as.integer(infectionHistory)
strains <- as.numeric(strains)
circulationTimeIndices <- seq_along(strains)-1
circulationTimeIndices <- as.integer(circulationTimeIndices)
sampleTimes <- as.numeric(sampleTimes)
dataIndices <- as.integer(dataIndices)
measurementMapIndices <- as.numeric(match(strainIsolationTimes1,strains)) - 1
nStrains <- length(strains)
pars["mu"] <- 3.75
pars["mu_short"]  <- 5.3
dat1 <- simulate_individual(pars, infectionHistory,sampleTimes,
                            dataIndices,strainIsolationTimes1,
                            measurementMapIndices,
                            antigenicMapLong,antigenicMapShort,strains)
dat1 <- as.data.frame(dat1)
colnames(dat1) <- c("sample","virus","titre")
#dat1$titre <- floor(dat1$titre)

r_likelihood <- function(expected, data, theta){
  largeI <- data > theta["MAX_TITRE"]
  smallI <- data <= 0
  restI <- data > 0 & data <= theta["MAX_TITRE"]
  
  large <- pnorm(theta["MAX_TITRE"], expected[largeI],theta["error"],lower.tail=FALSE,log.p=TRUE)
  small <- pnorm(1, expected[smallI],theta["error"],lower.tail=TRUE,log.p=TRUE)
  rest <- log(pnorm(data[restI]+1,expected[restI],theta["error"],lower.tail=TRUE,log.p=FALSE) - 
                pnorm(data[restI],expected[restI], theta["error"],lower.tail=TRUE,log.p=FALSE))
  return(sum(large, small, rest))
}
conciseInfHist <- infectionHistory[infectionHistory > 0]
conciseInfTimes <- strains[which(infectionHistory > 0)]
infectionMapIndices <- circulationTimeIndices[which(infectionHistory > 0)]
#samplingTime <- 2011*12

#pars["wane"] <- 0.03
#dat <- infection_model_indiv(pars, conciseInfHist,conciseInfTimes,infectionMapIndices,samplingTime,match(strains,strains)-1,antigenicMapLong,antigenicMapShort,nStrains)
liks <- NULL
data <- floor(dat1$titre)
system.time(
for(i in 1:1000){
  pars["mu"] <- (i-1)/100
  
  #y <- infection_model_indiv(pars, conciseInfHist,conciseInfTimes,infectionMapIndices,samplingTime,match(strains,strains)-1,antigenicMapLong,antigenicMapShort,nStrains)
  titres <- titre_data_individual(pars, infectionHistory, strains,circulationTimeIndices , sampleTimes, 
                                  dataIndices, measurementMapIndices,strainIsolationTimes1,
                           antigenicMapLong,antigenicMapShort,nStrains)
  liks[i] <- r_likelihood(titres,data,pars)
  #liks[i] <- likelihood_titre(titres, dat1$titre, pars)
  #liks[i] <- sum(dnorm(x=dat1$titre, mean=titres, sd=1,log=TRUE))
}
)
plot(liks)

bestPars <- pars
bestPars["mu"] <- (which.max(liks)-1)/100
titres <- titre_data_individual(bestPars, infectionHistory, strains,circulationTimeIndices , sampleTimes, dataIndices, 
                                measurementMapIndices,strainIsolationTimes1,
                                antigenicMapLong,antigenicMapShort,nStrains)
fitDat <- dat1
fitDat$titre <- titres

pars["mu"] <- 3.75
actualTitres <- titre_data_individual(pars, infectionHistory, strains,circulationTimeIndices , 
                                      sampleTimes, dataIndices, 
                                      measurementMapIndices,strainIsolationTimes1,
                                antigenicMapLong,antigenicMapShort,nStrains)
actualDat <- dat1
actualDat$titre <- actualTitres


infDat <- data.frame(x=strains[which(infectionHistory==1)])
sampDat <- data.frame(x=sampleTimes)
ggplot()+ 
  geom_line(data=fitDat, aes(x=virus,y=titre),col="green")+
  geom_line(data=actualDat,aes(x=virus,y=titre),col="purple") +
  geom_point(data=dat1,aes(x=virus,y=floor(titre)),col="blue",size=0.1) +
  geom_vline(data=infDat,aes(xintercept=x)) + 
  geom_vline(data=sampDat,aes(xintercept=x),col="red")+
  facet_wrap(~sample,ncol=1)


