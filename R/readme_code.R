## Load in example parameter values and antigenic map
data(example_par_tab) 
data(example_antigenic_map)

## Get all possible infection times
strain_isolation_times <- unique(example_antigenic_map$inf_years)

## Vector of strains that have titres (note only one representative strain per time)
sampled_viruses <- seq(min(strain_isolation_times), max(strain_isolation_times), by=2)

## Times at which serum samples can be taken
sampling_times <- 2010:2015

## Number of serum samples taken
n_samps <- 2

## Simulate some random attack rates
attack_rates <- runif(length(strain_isolation_times), 0.05, 0.15)

## Simulate a full serosurvey with these parameters
all_simulated_data <- simulate_data(par_tab=example_par_tab, group=1, n_indiv=50,    
                                  strain_isolation_times=strain_isolation_times,
                                  measured_strains=sampled_viruses,
                                  sampling_times=2010:2015, nsamps=n_samps, 
                                  antigenic_map=example_antigenic_map, 
                                  age_min=10,age_max=75,
                                  attack_rates=attack_rates, repeats=2)

## Pull out the simulated titre data and infection histories
titre_dat <- all_simulated_data$data
ages <- all_simulated_data$ages
example_inf_hist <- all_simulated_data$infection_histories
example_titre_dat <- merge(titre_dat, ages)

## Run the MCMC. We have to remove phi, as running prior version 2.
par_tab <- example_par_tab[example_par_tab$names != "phi",]
res <- run_MCMC(par_tab, example_titre_dat, example_antigenic_map, filename="test", version=2,
              mcmc_pars=c(save_block=1000,thin=10,thin_hist=100,move_sizes=3,swap_propn=0.5,inf_propn=0.5))

## Read in the MCMC chains and plot!
chain <- read.csv(res$chain_file)
inf_chain <- data.table::fread(res$history_file)
plot(coda::as.mcmc(chain[chain$sampno > 10000,c("mu","mu_short","wane")]))      


