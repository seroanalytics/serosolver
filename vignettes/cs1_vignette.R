## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo=TRUE,
  comment = "#>",
  eval=TRUE
)

## ----message=FALSE, warning=FALSE, r,eval=TRUE---------------------------
# Required to run serosolver
#devtools::install_github("seroanalytics/serosolver")
library(serosolver)
library(plyr)
library(data.table)

## Required for this analysis
library(reshape2)
library(foreach)
library(doParallel)
library(bayesplot)
library(coda)
library(ggplot2)
library(viridis)

# set up cluster
set.seed(1234)
cl <- makeCluster(5)
registerDoParallel(cl)

## Note that this vignette was generated on a Windows machine,
## and the setup for parallelisation is different on a Linux machine
## for Linux machine:
#  library(doMC)
#  library(doRNG)
#  registerDoMC(cores=5)

## ---- eval=TRUE----------------------------------------------------------
filename <- "case_study_1"
resolution <- 4 ## set to 4 for quarterly resolution
sample_years <- 2009:2012

serosolver::describe_proposals()
prior_version <- 2

## ------------------------------------------------------------------------
## Read in titre data
# unvaccinated
  input_dat_path <- system.file("extdata", "HKdata_h1n1_unvac.csv", package = "serosolver")
  input_dat <- read.csv(file = input_dat_path, header = TRUE)
# vaccinated
# input_dat_path2 <- system.file("extdata", "HKdata_h1n1_vac.csv", package = "serosolver")
# input_dat_vac <- read.csv(file = input_dat_path2, header = TRUE)

indivs <- unique(input_dat$individual) #all individuals

# Subset data for indivs
titre_dat <- input_dat[input_dat$individual %in% indivs,c("individual","virus","titre","samples","DOB")]
titre_dat$individual <- match(titre_dat$individual, indivs)

titre_dat <- unique(titre_dat)
titre_dat <- plyr::ddply(titre_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x),"group"=1))
print(head(titre_dat))

## ------------------------------------------------------------------------
## Read in raw coordinates
antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)

## Convert to form expected by serosolver
antigenic_map <- generate_antigenic_map(antigenic_coords, resolution)

## Restrict entries to years of interest. Entries in antigenic_map determine
## the times that individual can be infected ie. the dimensions of the infection
## history matrix.
antigenic_map <- antigenic_map[antigenic_map$inf_years>=(sample_years[1]*resolution+1) & antigenic_map$inf_years<=sample_years[4]*resolution,]

## Change all coordinates to the same as the first pair of coordinates (because we assume the virus is the same at every time point)
antigenic_map[2:nrow(antigenic_map),c('x_coord','y_coord')] <- antigenic_map[1,c('x_coord','y_coord')]
print(head(antigenic_map))

strain_isolation_times <- unique(antigenic_map$inf_years)

## ------------------------------------------------------------------------
par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(par_tab_path, stringsAsFactors=FALSE)

## Set parameters for beta and alpha to 1
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1,1)
## Maximum recordable log titre in these data is 9
par_tab[par_tab$names == "MAX_TITRE","values"] <- 9

## Remove phi parameters, as these are integrated out under prior version 2
par_tab <- par_tab[par_tab$names != "phi",]

## Fix all long term parameters to 0
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"fixed"] <- 1 #mu, tau, sigma1, and sigma2 are fixed
par_tab[par_tab$names %in% c("mu","tau","sigma1","sigma2"),"values"] <- 0 # set these values to 0

