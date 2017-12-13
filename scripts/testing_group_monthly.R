setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map.csv")
data <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)
ages <- read.csv("data/fluscape_ages.csv")
if(!is.null(n_indiv)){
  indivs <- sample(unique(data$individual),n_indiv)
  data <- data[data$individual %in% indivs,]
  ages <- ages[ages$individual %in% indivs,]
}
data$virus <- data$virus*12
data$samples <- data$samples*12

f <- create_post_func(parTab, data, antigenicMap, NULL)

liks <- NULL
for(i in 1:100){
  pars["mu"] <- (i-1)/10
  liks[i] <- f(pars, infectionHistories)
}
plot(liks)


sampleTimes <- c(2007,2014,2014)*12
strainIsolationTimesA <- sample(strains, 15,replace=FALSE)
strainIsolationTimes <- rep(strainIsolationTimesA, length(sampleTimes))
virusIndices <- match(strainIsolationTimes,strains)-1
dataIndices <- rep(length(strainIsolationTimesA), length(sampleTimes))
y <- simulate_individual(pars, 
                    infectionHistories[1,], 
                    sampleTimes, 
                    dataIndices, 
                    strainIsolationTimes, 
                    virusIndices, 
                    antigenicMapLong,
                    antigenicMapShort, 
                    strains)



infTimes <- NULL
indiv <- NULL
for(i in 1:nrow(infectionHistories[1:5,])){
  tmp <- strains[which(infectionHistories[i,]==1)]
  infTimes <- c(infTimes, tmp)
  indiv <- c(indiv,rep(i, length(tmp)))
}
infDat <- data.frame(individual=indiv, x=infTimes)

ggplot(y) + geom_line(aes(x=virus,y=titre)) + geom_vline(data=infDat,aes(xintercept=x))+
  facet_grid(individual~samples) +
  theme_bw()

