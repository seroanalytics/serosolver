###################################################
## Simulate data ##
###################################################
devtools::load_all("Q:/serosolver-parallel_tempering")

path <- "Q:/serosolver_testing/"
n_indiv <- 500   # number of individuals
n_samp <- 4      # number of samples per individual
buckets <- 1     # time buckets (1 - annumal, 2 - semi-annual, 4 - quarterly, 12 - monthly)
repeats <- 1     # number of sample repeats
yearMin <- 1968  # mminimum year for circulating viruses
yearMax <- 2015  # maximum year for circulating viruses
sampMin <- 2010  # minimum time for samples to be collected
sampMax <- 2015  # maximum time for samples to be collected 
ageMin <- 10     # minimum age
ageMax <- 75     # maximum age
dataType <- 1    # 1 - heterologous, 2 - homologous

# Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
parTab <- read.csv(paste0(path,"parTab_base.csv"),stringsAsFactors=FALSE)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
parTab[parTab$names == "wane","values"] <- 1
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets

parTab <- parTab[parTab$names != "lambda",] # remove lambda row of parameter table

yearMin <- yearMin*buckets
yearMax <- yearMax*buckets
samplingTimes <- seq(sampMin, sampMax, by=1)
ageMin <- ageMin*buckets
ageMax <- ageMax*buckets

# Read in and generate the antigenic map to read strain relationships from
antigenicMap <- read.csv(paste0(path,"Fonville2014AxMapPositionsApprox.csv"),stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- fit_dat[fit_dat$inf_years >= yearMin & fit_dat$inf_years <= yearMax,]
if (dataType==2){ # for H1N1 data make all coordinates the same 
  fit_dat[2:nrow(fit_dat),c('x_coord','y_coord')] <- fit_dat[1,c('x_coord','y_coord')]
}

strainIsolationTimes <- unique(fit_dat$inf_years)

simInfPars=c("mean"=0.15/buckets,"sd"=0.5,"bigMean"=0.5/buckets/2,"logSD"=1)
# use SIR
# attackRates  <- simulate_ars_buckets(strainIsolationTimes, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
# attackRates  <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
# attackRates  <- simulate_ars_spline(strainIsolationTimes, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"], knots,theta)
attackRates <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])

# Simulate data
dat <- simulate_data(par_tab = parTab,
                     group = 1, 
                     n_indiv = n_indiv, 
                     buckets = buckets,
                     strain_isolation_times = strainIsolationTimes, 
                     sampling_times = samplingTimes, 
                     nsamps = n_samp,
                     antigenic_map = fit_dat,
                     titre_sensoring = 0,
                     age_min = ageMin*buckets, 
                     age_max = ageMax*buckets,
                     attack_rates = attackRates,
                     repeats = repeats,
                     mu_indices = NULL,
                     measurement_indices = NULL,
                     add_noise = TRUE)
  
titreDat <- merge(dat$data,dat$ages,by='individual') # add DOB column

if (dataType==2){
  mydat <- dat$data[dat$data$virus == yearMin,]
  titreDat <- merge(mydat,dat$ages,by='individual') # add DOB column
}

# write to .csv file
write.csv(titreDat,file = paste0(path,"titreDat.csv"),row.names = FALSE)
