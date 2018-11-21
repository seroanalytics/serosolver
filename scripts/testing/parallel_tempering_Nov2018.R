######################
## JAMES HAY 13.08.2018 - jameshay218@gmail.com
## This script fits the serosolver antibody kinetics model to the fluscape HI titre data
## This particular script uses a gibbs proposal step to resample infection histories
## which integrates out the annual force of infection parameters.
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

source("~/net/home/serosolver/scripts_Oct2018/submit_run_PT.R")

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
output_wd <- "testing_cluster_script_Oct2018"
runName <- "parallel_tempering_test"

## How many individuals to fit to?
n_indiv <-50
indivs <- 1:25
indivs <- NULL
buckets <- 1

alpha <- beta <- 1

realData <- FALSE
infHist_file <- "~/net/home/serosolver/data_Oct2018/fluscape_sim_annual/fluscape_sim_annual_infHist_1.csv"
AR_file <- "~/net/home/serosolver/data_Oct2018/fluscape_sim_annual/fluscape_sim_annual_AR_1.csv"
parTab_file <- "~/Documents/Fluscape/serosolver/inputs/parTab_base.csv"
n_alive_file <- "NULL"
titreDat_file <- "~/net/home/serosolver/data_Oct2018/fluscape_sim_annual/fluscape_sim_annual_dat_1.csv"
measurement_indices <- mu_indices <-  NULL
n_alive <- NULL
realPars <- NULL
random_effects_measurement <- FALSE
create_prior_f <- NULL
mcmcPars <- NULL
sampd <- NULL

## Read in and generate the antigenic map to read strain relationships from
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
strainIsolationTimes <- unique(fit_dat$inf_years)

temperatures <- c(seq(1,1.2,by=0.05),seq(1.3,1.9,by=1),seq(2,10,by=1))
#Rprof(tmp<-tempfile())
res <- submit_serosolver_run_PT(output_wd=output_wd, runName=runName, chainNo=1,buckets=buckets,
                             parTab_file=parTab_file, titreDat_file=titreDat_file,
                             infHist_file=infHist_file, AR_file=AR_file,alpha=alpha,beta=beta,
                             version=2,solve_likelihood=TRUE,
                             create_prior_f=NULL,
                             fit_dat=fit_dat, n_alive_file=n_alive_file,use_n_alive=FALSE,
                             realData=realData,
                             mu_indices=mu_indices,
                             measurement_indices = measurement_indices,measurement_random_effects=FALSE,
                             indivs=indivs,
                             mcmcPars=list("save_block"=100,"thin"=1,"thin2"=10,"iterations"=50000,"adaptive_period"=20000,
                                        "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0,"year_swap_propn"=0.2,swapPropn=0.5,
                                        "temperature"=temperatures, "parallel_tempering_iter"=5,
                                        "inf_propn"=0.05,"histSampleProb"=0.5,moveSize=2,histOpt=0,popt=0.44,opt_freq=1000),
                             inf_propn=0.5,
                             histSampleProb = 0.5,
                             hist_switch_prob=0.5,
                             year_swap_propn=0.2,
                             inf_hist_func=3,
                             message_slack=FALSE,
                             message_slack_pars=NULL,
                             save_to_google=FALSE,
                             google_pars=list(key="15yNgMEXq9yH8OCekF7nnmQC0JwTocMhftnz72-nSvQM",
                                              token=NULL)
  )
  #Rprof(NULL)
  #summaryRprof(tmp)

