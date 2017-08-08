# - - - - - - - - - - - - - - - - - - - - -
# Model of serological dynamics
# github.com/adamkucharski/serology-model
#
# Main execution code
# - - - - - - - - - - - - - - - - - - - - -

# Load libraries
library(reshape2)
library(mvtnorm)
library(MASS)
library(coda)
library(RColorBrewer)
library(magrittr)
library(plot3D)
library(colorspace)

library(foreach)
library(doMC)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()


setwd("~/Documents/serology-model/R/") # set R directory
rm(list=ls(all=TRUE)) # Clear environment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

source("sero_functions.R")
source("posterior_analysis.R")

compile.c() # Compile c code

antigenic_map = read.csv("../data/antigenic_map.csv") # load antigenic locations


# >>> Run code up to here to set everything up


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE MODEL
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation study and inference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Generate simulated data and infer parameters -- simulation parameters are defined in sero_functions.R
# flu.type defines which dataset format (i.e. test strains, test years) the simulated data will produce

kk=1
simulation.infer(seed_i=kk,mcmc.iterations=5e4, flu.type="H3HN", strain.fix=T,
                   fit.map=antigenic_map,vp1=0.4,linearFn=T) 



# Plot convergence for MCMC chains for H3 Vietnam simulated data
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",simDat=T,year_test=c(2007:2012),
                            linearFn=T,loadpick = c(1:4))


# Plot simulation study posteriors and attack rate comparisons for simulation plots
for(kk in 1:3){
  
  plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),flu.type="H3HN",year_test=c(2007:2012),plotmap=F,fr.lim=T,linearFn=T) #H3 Vietnam
  #plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),year_test=c(2009),plotmap=F,fr.lim=T) #H3 FluScape
  
}

# Plot simulation study titres against inferred model
kk=1
plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3",simDat=T,year_test=c(2007:2012),btstrap=10,linearFn=T) #H3 Vietnam
#plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3FS",simDat=T,year_test=c(2009),btstrap=10) #H3 FluScape

# Plot estimated vs true attack rates from simulated data (Fig 3B)

plot.multi.true.vs.estimated(burnCut=0.25,flu.type="H3HN",simDat=T,loadpick=c(1:10))

