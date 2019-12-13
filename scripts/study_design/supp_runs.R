###Supplementary runs 
code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

buckets <- 1

## Antigenic map for cross reactivity parameters
antigenic_map <- read.csv(file.path(code.dir,"data/fonville_map_approx.csv"),
                          stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, buckets=1)

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968,]
strain_isolation_times <- unique(fit_dat$inf_years)

## Read in parameter table
par_tab <- read.csv(file.path(code.dir,"inputs/par_tab_study_design.csv"))

## Specify paramters controlling the MCMC procedure
mcmc_pars <- c("iterations"=500000,"adaptive_period"=500000, "burnin"=50000,
               "hist_sample_prob"=0.5,"thin"=50,"thin_hist"=200,
               "swap_propn"=0.5,"hist_switch_prob"=0.2,"year_swap_propn"=0.5)

for(filename in c("titre_boost_dat", "titre_AR_dat")){
  titre_dat <- read.csv(paste("~/Documents/Github/serosolver/data/study_design/",filename,".csv",sep=""))
  res <- run_MCMC(par_tab = par_tab, titre_dat = titre_dat, 
                  antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                  mvr_pars = NULL, start_inf_hist = NULL, filename=paste("chains/", filename,sep=""),
                  CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = NULL,
                  version = 2,  
                  fast_version = TRUE)
}

##Plot of boosting parameter estimates for boosting
chain <- read.csv("~/Documents/GitHub/serosolver/chains/titre_boost_dat_chain.csv")
chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]

min_mu_short <- par_tab$values[which(par_tab$names == "mu_short")]
max_mu_short <- par_tab$values[which(par_tab$names == "mu_short")] * 2

png("mu_short_supp.png",width=3000,height=2000,res=300,units='px')
hist(chain$mu_short, xlab =  expression(paste("short boost, ",mu[2])), freq = F,
     main = "", xlim = c(min(min_mu_short, min(chain$mu_short)), max(max_mu_short, max(chain$mu_short))))
abline(v = min_mu_short, lwd = 2, col = "gray", lty = 2)
abline(v = max_mu_short, lwd = 2, col = "gray", lty = 2)
abline(v =(min_mu_short + max_mu_short)/2, lwd = 2, col = "gray")
dev.off()

##Attack rate plot for AR data
filename <- "titre_AR_dat"
titre_dat <- read.csv(paste("~/Documents/Github/serosolver/data/study_design/",filename,".csv",sep=""))

inf_chain <- data.table::fread(paste("chains/",filename,"_infection_histories.csv",sep=""))
inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]

#the first 1 : 125 ids are the youngest age
p_1 <- plot_attack_rates(infection_histories = inf_chain[which(inf_chain$i %in% 1:125), ], titre_dat_AR[which(titre_dat_AR$individual %in% 1:125),], strain_isolation_times, by_val = 10)+
  geom_hline(yintercept = 0.3, linetype="dashed", color = "gray") +
  ggtitle("Individuals aged 1 - 10 years")

p_2 <- plot_attack_rates(infection_histories = inf_chain[which(inf_chain$i %in% 126:250), ], titre_dat_AR[which(titre_dat_AR$individual %in% 126:250),], strain_isolation_times, by_val = 10)+
  geom_hline(yintercept = 0.15, linetype="dashed", color = "gray") +
  ggtitle("Individuals aged 11 - 20 years")

p <- plot_attack_rates(infection_histories = inf_chain, titre_dat_AR, strain_isolation_times, by_val = 10)+
  geom_hline(yintercept = 0.225, linetype="dashed", color = "gray") +
  ggtitle("All individuals")


p_all <- do.call(grid.arrange,c(list(p_1, p_2, p), nrow = 1))
ggsave("AR_supp.png",p_all, width = 20, height = 6)
