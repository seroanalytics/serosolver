# - - - - - - - - - - - - - - - - - - - - -
# Main execution code
#
# Model of serological dynamics
# github.com/adamkucharski/serosolver
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
if(Sys.info()["user"]=="adamkuchars") { setwd("~/Documents/serosolver/R") }  # set R directory


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
# Parameters currently specified within the function (see "sero_functions.R")

kk=1
simulation.infer(seed_i=kk, # Set seed
                 mcmc.iterations=20,  # MCMC iterations
                 flu.type="H3HN", # Flu name
                 fit.map=antigenic_map, # Antigenic map input
                 fix.param=NULL, # Add noise to input parameters?
                 vp1=0.4 # Proportion of infection histories to resample each step
                 ) 


# Plot convergence for MCMC chains for simulated data
plot.multi.chain.posteriors(burnCut=0.25,
                            flu.type="H3HN",
                            simDat=T,
                            year_test=c(2007:2012),
                            loadpick = c(1) # Which simulation to plot
                            )


# Plot simulation study posteriors and attack rate comparisons for simulation plots
for(kk in 1){
  
  plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),flu.type="H3HN",year_test=c(2007:2012),plotmap=F,fr.lim=T)

}

# Plot simulation study titres against inferred model
kk=1
plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3",simDat=T,year_test=c(2007:2012),btstrap=10)


