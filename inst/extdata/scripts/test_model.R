## Test antibody kinetics models
setwd("~/Documents/GitHub/serosolver")
devtools::load_all("~/Documents/GitHub/serosolver")
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)


par_tab <- read.csv("~/Documents/GitHub/serosolver/inst/extdata/par_tab_base.csv")

n_indiv <- 9
possible_exposure_times <- seq(1,25,by=1)

example_inf_hist <- simulate_infection_histories(runif(length(possible_exposure_times),0,0.15),possible_exposure_times,possible_exposure_times,rep(1,5))[[1]]

example_inf_hist <- matrix(0, nrow=n_indiv,ncol=25)
example_inf_hist[,1] <- 1
example_inf_hist[,10] <- 1

antibody_data <- as.data.frame(expand_grid(individual=1:n_indiv, biomarker_group=1,sample_time=possible_exposure_times,biomarker_id=25,repeat_number=1,birth=1,measurement=0))
#antibody_data <- antibody_data %>% mutate(demographic_group = if_else(individual %in% 1:3, 1, 2))
antibody_data <- antibody_data %>% mutate(age_group = if_else(individual %in% 1:6, if_else(individual %in% 1:3, 1,2), 3))
antibody_data <- antibody_data %>% mutate(sex = if_else(individual %in% 1:5, 1,2))

#antibody_data$demographic_group <- 1
#antibody_data <- antibody_data %>%mutate(demographic_group = if_else(individual %in% 1:6, if_else(individual %in% 1:3, 1,2), 3))

f <- create_posterior_func(par_tab,antibody_data,NULL,function_type=4,
                           possible_exposure_times = possible_exposure_times,
                           start_level="none",use_demographic_groups=c("sex"))

antibody_data$measurement <- f(par_tab$values, example_inf_hist)
ggplot(antibody_data) + 
  geom_line(aes(x=sample_time,y=measurement)) +
  facet_wrap(~individual)

## Create a set of scaling parameters for each demographic variable
demographics <- antibody_data %>% select(age_group, sex) %>% distinct()
n_demographic_groups <- nrow(demographics)
n_pars <- nrow(par_tab[par_tab$par_type == 1,])
## Add stratification column to par_tab? Tell model which variables to stratify this parameter by e.g., sex, age group etc
par_tab$stratification <- NA
par_tab[par_tab1$names %in% c("boost_long"),"stratification"] <- c("age_group, sex")
par_tab[par_tab1$names %in% c("wane_short"),"stratification"] <- c("age_group")


setup_stratification_table <- function(par_tab, demographics){
  n_pars <- nrow(par_tab[par_tab$par_type == 1,])
  
  ## Creates an estimated parameter entry for each 
  stratifications <- unique(unlist(sapply(par_tab$stratification,function(x) strsplit(x,", "))))
  stratifications <- stratifications[!is.na(stratifications)]
  
  scale_table <- vector(mode="list",length=length(stratifications))
  for(i in seq_along(stratifications)){
    stratification <- stratifications[i]
    scale_table[[i]] <- matrix(1, nrow=length(unique(demographics[,stratification])),ncol=n_pars)
  }
  names(scale_table) <- stratifications
  
  ## First row is base case
  ## Subsequent rows, check if they are estimated. If so, flag as estimated. Otherwise, is fixed
  index <- 2
  strat_par_names <- NULL
  for(j in 1:nrow(par_tab[par_tab$par_type == 1,])){
    stratification_par <- par_tab$stratification[j]
    if(!is.na(stratification_par)){
      strats <- strsplit(stratification_par,", ")[[1]]
      for(strat in strats){
        n_groups <- length(unique(demographics[,strat]))
        for(x in 2:n_groups){
          scale_table[[strat]][x,j] <- index
          strat_par_names[[index]] <- paste0(par_tab$names[j],"_coef_",strat,"_",x)
          index <- index + 1
        }
      }
    }
  }
  scale_pars <- c(rnorm(index-2,0,0.1))
  names(scale_pars) <- unlist(strat_par_names)
  return(list(scale_table, scale_pars))
}

stratification_pars <- setup_stratification_table(par_tab, demographics)
scale_table <- stratification_pars[[1]]
scale_pars <- stratification_pars[[2]]
## Say I have one parameters, mu
pars <- par_tab[par_tab$par_type == 1,"values"]
names(pars) <- par_tab[par_tab$par_type == 1,"names"]
pars_all <- c(pars, scale_pars)
scale_par_indices <- (length(pars)+1):(length(pars) + length(scale_pars))
theta_indices <- 1:length(pars)


transform_parameters <- function(pars, scale_table, theta_indices,scale_par_indices,demographics){
  scale_pars <- c(0,pars[scale_par_indices])
  theta_pars <- pars[theta_indices]
  
  n_demographic_groups <- nrow(demographics)
  n_strats <- ncol(demographics)
  
  theta <- matrix(0, nrow=n_demographic_groups,ncol=length(theta_pars))
  ## For each parameter
  for(i in seq_along(theta_pars)){
    ## Need to calculate the value for each demographic group
    for(j in 1:n_demographic_groups){
      ## Each demographic group has its own values for each stratification -- sum contribution of all of these
      scales <- 0
      ## For each stratification
      for(x in 1:n_strats){
        ## Get value of variable for this stratification
        tmp_strat <- demographics[j,x]
        ## Get index in scale_table for this stratification level for this parameter
        par_index <- scale_table[[x]][tmp_strat,i]
        scales <- scales + scale_pars[par_index]
      }
      theta[j,i] <- exp(log(theta_pars[i]) + scales)
    }
  }
  theta
}
theta <- transform_parameters(pars_all, scale_table, theta_indices, scale_par_indices,demographics)
microbenchmark::microbenchmark(transform_parameters(pars_all, scale_table, theta_indices, scale_par_indices,demographics))
