## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("seroanalytics/serosolver", build_vignettes=TRUE)
#  library(serosolver)

## ----eval=FALSE----------------------------------------------------------
#    titre_dat <- data.frame(individual=c(rep(1,4),rep(2,4)),
#                            samples=c(8039,8040,8044,8047,8039,8041,8045,8048),
#                            virus=c(rep(8036,8)),
#                            titre=c(0,0,7,7,0,5,6,5),
#                            run=c(rep(1,8)),
#                            DOB=c(rep(8036,8)),
#                            group=c(rep(1,8))
#                           )
#    knitr::kable(head(titre_dat))

## ----echo=FALSE,results='asis'-------------------------------------------
  options(scipen=999)
  par_tab_path <- system.file("extdata", "parTab_base.csv", package = "serosolver")
  par_tab <- read.csv(par_tab_path, stringsAsFactors=FALSE)
  par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1,1)
  knitr::kable(par_tab)

## ----eval=FALSE----------------------------------------------------------
#  # generate starting parameter values
#    start_tab <- par_tab
#      for(i in 1:nrow(start_tab)){
#        if(start_tab[i,"fixed"] == 0){
#          start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"],
#                                    start_tab[i,"upper_start"])
#        }
#      }
#  

## ----antigenic_map, echo=FALSE, fig.cap="Assumed antigenic locations of historical strains in model between 1968 and 2012 (from Kucharski *et al.* @kucharski2018).", out.width = '75%'----
knitr::include_graphics("antigenic_map.tiff")

## ---- eval=FALSE---------------------------------------------------------
#  ## Read in raw coordinates
#  antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
#  antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)
#  
#  ## Convert to form expected by serosolver
#  antigenic_map <- generate_antigenic_map(antigenic_coords, buckets = 4)
#  
#  # unique strain circulation times
#    strains_isolation_times <- unique(antigenic_map$inf_years)

## ---- eval=FALSE---------------------------------------------------------
#    unique_indiv <- titre_dat[!duplicated(titre_dat$individual),]
#    ageMask <- create_age_mask(unique_indiv$DOB, strains_isolation_times)

## ---- eval=FALSE---------------------------------------------------------
#    start_inf <- setup_infection_histories_new_2(titre_dat, strains_isolation_times)

## ---- eval=FALSE---------------------------------------------------------
#  # run MCMC
#    res <- run_MCMC(par_tab = start_tab, titre_dat = titre_dat, antigenic_map = antigenic_map,
#                    start_inf_hist = start_inf, CREATE_POSTERIOR_FUNC = create_posterior_func,
#                    version = 2)

## ---- eval=FALSE---------------------------------------------------------
#     mcmcPars <- c("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=1000,"thin"=1,
#                   "adaptive_period"=10000, "save_block"=1000, "thin2"=10, "histSampleProb"=0.5,
#                   "switch_sample"=2, "burnin"=0, "inf_propn"=0.5, "moveSize"=5, "histOpt"=1,
#                   "swapPropn"=0.5,"hist_switch_prob"=0.5,"year_swap_propn"=0.5)
#  
#    res <- run_MCMC(par_tab = start_tab, titre_dat = titre_dat, antigenic_map = antigenic_map,
#                    mcmcPars = mcmcPars, start_inf_hist = start_inf,
#                    CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2)

## ----eval=FALSE----------------------------------------------------------
#  # Density/trace plots
#    chain1 <- read.csv(res$chain_file)
#    chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
#    plot(coda::as.mcmc(chain1))

## ----eval=FALSE----------------------------------------------------------
#  # Read in infection history file from MCMC output
#    inf_chain <- data.table::fread(res$history_file,data.table=FALSE)
#  # Remove adaptive period and burn-in
#    inf_chain <- inf_chain[inf_chain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
#    inf_chain1 <- setDT(inf_chain)
#  # Define year range
#    xs <- min(strains_isolation_times):max(strains_isolation_times)
#  # Plot inferred attack rates
#    arP <- plot_attack_rates(infectionHistories = inf_chain1, dat = titre_dat, ages = titre_dat[,c('individual', 'DOB')], yearRange = xs)

## ---- eval=FALSE---------------------------------------------------------
#  # Plot infection histories
#    IH_plot <- plot_infection_histories(chain = chain1, infection_histories = inf_chain,
#                                        titre_dat = titre_dat, individuals = c(1:5),
#                                        antigenic_map = antigenic_map, par_tab = start_tab1)

## ----eval=FALSE----------------------------------------------------------
#  # Plot inferred antibody titres
#    titre_preds <- get_titre_predictions(chain = chain1, infection_histories = inf_chain,
#                                         titre_dat = titre_dat, individuals = c(1:311),
#                                         antigenic_map = antigenic_map, par_tab = start_tab1)
#    to_use <- titre_preds$predictions
#  
#    titre_pred_p <- ggplot(to_use[to_use$individual %in% 151:200,])+
#                    geom_line(aes(x=samples, y=median))+
#                    geom_point(aes(x=samples, y=titre))+
#                    geom_ribbon(aes(x=samples,ymin=lower, ymax=upper),alpha=0.2,col='red')+
#                    facet_wrap(~individual)

## ----eval=FALSE----------------------------------------------------------
#    y <- generate_cumulative_inf_plots(inf_chain,burnin = 0,51:75,nsamp=100,
#                                       strain_isolation_times = strains_isolation_times,
#                                       pad_chain=FALSE,number_col = 2,subset_years = NULL)

## ----eval=FALSE----------------------------------------------------------
#  # Table of antibody kinetics parameters
#    myresults <- matrix(c(rep(0,3*7)),nrow=3)
#    rownames(myresults) <- c("mu_short","wane","error")
#    colnames(myresults) <- c("mean","sd","2.5%","25%","50%","75%","97.5%")
#  
#    myresults[,"mean"] <- round(apply(chain1[,c("mu_short","wane","error")],2,mean),3)
#    myresults[,"sd"] <- round(apply(chain1[,c("mu_short","wane","error")],2,sd),3)
#    myresults[,3:7] <- t(round(apply(chain1[,c("mu_short","wane","error")],2,
#                                     quantile,probs=c(0.025,0.25,0.5,0.75,0.975)),3))

