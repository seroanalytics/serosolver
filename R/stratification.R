
#' Add scaling parameters to par_tab
#'
#' Adds entries to par_tab to stratify parameters by requested stratification levels. antibody_data and timevarying_demographics are used to find how many levels of each stratification have been requested
#' @param par_tab the parameter table, including a column called stratification which is NA if no stratification is requested, for with a character value matching an entry in antibody_data or timevarying_demographics
#' @param antibody_data the antibody data data frame, see \code{\link{example_antibody_data}}. If NULL, then uses timevarying_demographics
#' @param timevarying_demographics a data frame of timevarying demographics, with columns individual, time and any stratification variables. If NULL, then uses antibody_data
#' @param scale_par_lower the lower bound of any used scale parameters
#' @param scale_par_upper the upper bound of any used scale parameters
#' @return the updated par_tab
#' @export
add_scale_pars <- function(par_tab, antibody_data, timevarying_demographics=NULL, scale_par_lower=-25,scale_par_upper=25){
  ## Check if timevarying demographics are used. If so, then use these to create demographic table and add scale parameters
  if(!is.null(timevarying_demographics)){
    timevarying_demographics <- timevarying_demographics %>% dplyr::arrange(individual, time)
    timevarying_demographics <- as.data.frame(timevarying_demographics)
    demographics <- create_demographic_table(timevarying_demographics,par_tab)
    ## Otherwise, base stratification parameters on the antibody data
  } else {
    demographics <- create_demographic_table(antibody_data,par_tab)
  }
  stratification_pars <- setup_stratification_table(par_tab, demographics)
  scale_table <- stratification_pars[[1]]
  scale_pars <- stratification_pars[[2]]
  
  ## If stratifying by anything
  if(length(scale_pars) > 0){
    ## Creates an estimated parameter entry for each 
    strsplit1 <- function(x){
      if(!is.na(x)){
        return(strsplit(x,"_"))
      } else {
        return(NA)
      }
    }
    
    names_split <- lapply(names(scale_pars), function(x) strsplit1(x))
    biomarker_groups <- unlist(lapply(names_split, function(x) as.numeric(x[[1]][which(x[[1]] == "biomarker") + 1])))
    
    tmp_par_tab <- data.frame(names=names(scale_pars), values=rnorm(length(scale_pars),0,0.1),fixed=0,
                              lower_bound=scale_par_lower,upper_bound=scale_par_upper,
                              lower_start=-0.25,upper_start=0.25,
                              par_type=4,biomarker_group=biomarker_groups,stratification=NA)
    par_tab <- bind_rows(par_tab, tmp_par_tab)
  }
  return(par_tab)
}


#' Create a demographics data frame of unique stratification levels
#'
#' Takes all stratifications requested in par_tab and checks antibody_data for all unique levels of that stratification. If no stratifications are requested, then returns a data frame with a single row and a column called "all" with value 0.
#' @param antibody_data the antibody data data frame, see \code{\link{example_antibody_data}}. 
#' @param par_tab the parameter table, including a column called stratification which is NA if no stratification is requested, for with a character value matching an entry in antibody_data
#' @return a data frame of unique stratification level combinations. Each column is a stratification variable, and each row is a unique combination of levels
#' @export
create_demographic_table <- function(antibody_data, par_tab){
  strsplit1 <- function(x){
    if(!is.na(x)){
      return(strsplit(x,", "))
    } else {
      return(NA)
    }
  }
  ## Skip any infection history prior parameters
  skip_pars <- c("infection_model_prior_shape1","infection_model_prior_shape2")
  
  ## Creates an estimated parameter entry for each 
  stratifications <- unique(unlist(sapply(par_tab[!(par_tab$names %in% skip_pars),"stratification"],function(x) strsplit1(x))))
  stratifications <- stratifications[!is.na(stratifications)]
  
  if(any(!stratifications %in% colnames(antibody_data))){
    stop(paste0("Trying to stratify by variable in par_tab, but this variable not in antibody_data: ", 
                paste(stratifications[!stratifications %in% colnames(antibody_data)], collapse=", ")
       ))
  }
  if(length(stratifications) == 0){
    demographics <- data.frame(all=0)
  } else {
    demographics <- antibody_data %>% select(all_of(stratifications)) %>% distinct()
    if(any(apply(demographics, 2, function(x) length(unique(x))) < 2)){
      message("Error - trying to stratify by variable in par_tab, but <2 levels for this variable in antibody_data")
    }
  }  
  
  demographics <- demographics %>% dplyr::arrange(across(everything()))
  return(demographics)
}



#' Create indexing table for parameter stratifications
#'
#' Creates an indexing table to tell serosolver which scale parameter to add to which model parameter.
#' @param par_tab the parameter table, including a column called stratification which is NA if no stratification is requested, for with a character value matching an entry in antibody_data or timevarying_demographics
#' @param unique_demographic_combinations table of unique stratification combinations. Columns give stratification variable, rows give each unique combinations
#' @return a list with two entries: 1) a list of matrices with entries for each stratification level, with number of columns equal to the number of parameters, and number of rows equal to the number of stratification levels. Each entry is the index of the scale parameter to use for that stratification level and parameter; 2) a vector of scale parameters
#' @export
setup_stratification_table <- function(par_tab, unique_demographic_combinations){
  unique_demographic_combinations <- as.data.frame(unique_demographic_combinations)
  use_par_tab <- par_tab[par_tab$par_type %in% c(1,3),]
  n_pars <- nrow(use_par_tab)
  ## Creates an estimated parameter entry for each 
  strsplit1 <- function(x){
    if(!is.na(x)){
      return(strsplit(x,", "))
    } else {
      return(NA)
    }
  }
  
  ## Skip any infection history prior parameters
  skip_pars <- c("infection_model_prior_shape1","infection_model_prior_shape2")
  ## Creates an estimated parameter entry for each 
  stratifications <- unique(unlist(sapply(par_tab[!(par_tab$names %in% skip_pars),"stratification"],function(x) strsplit1(x))))
  stratifications <- stratifications[!is.na(stratifications)]
  
  scale_table <- vector(mode="list",length=length(stratifications))
  
  ## If no stratifications, create dummy table
  if(nrow(unique_demographic_combinations) == 1){
    scale_table[[1]] <- matrix(0, nrow=1,ncol=n_pars)
    scale_pars <- NULL
  } else {
    for(i in seq_along(stratifications)){
      stratification <- stratifications[i]
      scale_table[[i]] <- matrix(0, nrow=length(unique(unique_demographic_combinations[,stratification])),ncol=n_pars)
    }
    names(scale_table) <- stratifications
    
    ## First row is base case
    ## Subsequent rows, check if they are estimated. If so, flag as estimated. Otherwise, is fixed
    index <- 2
    strat_par_names <- NULL
    
    ## Skip any infection history prior parameters
    skip_pars <- c("infection_model_prior_shape1","infection_model_prior_shape2")
    
    for(j in 1:nrow(use_par_tab)){
      stratification_par <- use_par_tab$stratification[j]
      if(!is.na(stratification_par) & !(use_par_tab$names[j] %in% skip_pars)){
        strats <- strsplit(stratification_par,", ")[[1]]
        for(strat in strats){
          unique_demo_strats <- unique(unique_demographic_combinations[,strat])
          unique_demo_strats_names <- unique_demo_strats[!is.na(unique_demo_strats)]
          n_groups <- length(unique_demo_strats_names)
          for(x in 2:n_groups){
            scale_table[[strat]][x,j] <- index-1
            strat_par_names[[index]] <- paste0(use_par_tab$names[j],"_biomarker_",use_par_tab$biomarker_group[j],"_coef_",strat,"_",unique_demo_strats_names[x])
            index <- index + 1
          }
        }
      }
    }
    scale_pars <- c(rnorm(index-2,0,0.1))
    names(scale_pars) <- unlist(strat_par_names)
    
  }
  return(list(scale_table, scale_pars))
}


#' Get stratification groups used to group infection history model
#' 
#' Takes antibody_data and demographics and returns a list with the stratification groups used to group the antibody kinetics model. If timevarying_demographics are provided, then these are used to create the demographic groups. If not, then antibody_data is used.
#' @param par_tab the parameter table. See \code{\link{example_par_tab}} for an example.
#' @param antibody_data the antibody data data frame. See \code{\link{example_antibody_data}} for an example.
#' @param timevarying_demographics optional. a data frame of timevarying demographics, with columns individual, time and any stratification variables. If NULL, then uses antibody_data
#' @param demographic_groups optional. a data frame of demographic groups, with columns for each demographic group and rows for each unique combination of demographic groups. If NULL, then create this from antibody_data or timevarying_demographics
#' @return a list with three entries: 1) the names of the demographic groups used, 2) a data frame of demographic groups, with each row a unique combination of demographic groups, and 3) a boolean indicating whether timevarying demographics were used
#' @export
get_demographic_groups <- function(par_tab, antibody_data, timevarying_demographics=NULL,demographic_groups=NULL){
  ## Setup data vectors and extract
  if(!is.null(timevarying_demographics)){
    timevarying_demographics <- timevarying_demographics %>% arrange(individual, time)
    timevarying_demographics <- as.data.frame(timevarying_demographics)
    use_timevarying_demographics <- TRUE
    if(is.null(demographic_groups)){
      demographic_groups <- create_demographic_table(timevarying_demographics,par_tab)
    }
  } else {
    use_timevarying_demographics <- FALSE
    if(is.null(demographic_groups)){
      demographic_groups <- create_demographic_table(antibody_data,par_tab)
    }
  }
  use_demographic_groups <- colnames(demographic_groups)
  if(length(use_demographic_groups) == 1 && use_demographic_groups == "all") use_demographic_groups <- NULL
  demographic_groups <- demographic_groups %>% dplyr::arrange(across(everything()))
  return(list(use_demographic_groups=use_demographic_groups, demographic_groups=demographic_groups,timevarying_demographics=use_timevarying_demographics))
}

#' Internal function -- merges antibody_data and demographics such that antibody_data has correct variables
align_antibody_demographic_dat <- function(antibody_data, demographics=NULL){
  
  if(!is.null(demographics)){
    overlapping_colnames <- intersect(colnames(antibody_data),colnames(demographics))
    overlapping_colnames <- overlapping_colnames[!(overlapping_colnames %in% c("individual","birth"))]
    if(length(overlapping_colnames) > 0){
      message(paste0("Warning: antibody_data and demographics have overlapping column names: ", paste(overlapping_colnames, collapse=", "), ". Using the values given in demographics."))
    }
    if("time" %in% colnames(demographics)){
      antibody_data <- suppressMessages(antibody_data %>% select(-overlapping_colnames) %>% left_join(demographics %>% dplyr::rename(sample_time = time)))
    } else {
      antibody_data <- suppressMessages(antibody_data %>% select(-overlapping_colnames) %>% left_join(demographics))
    }
  }
  antibody_data
}

#' Align stratification levels for serosolver infection history and antibody kinetics models
#' 
#' 
add_stratifying_variables <- function(antibody_data, timevarying_demographics=NULL, par_tab, use_demographic_groups=NULL){
    # Any stratification of population attack rates?
  ## Pull out any parameters related to attack rates
  population_group_strats <- par_tab %>% filter(names %like% "infection_model_prior" | names == "phi") %>% 
    pull(stratification) %>% unique()
  
  if(length(population_group_strats) > 1){
    stop("Error - trying to stratify infection model parameters by different variables, All entries for infection_model_prior or phi should have the same stratifications")
  }
  
  strsplit1 <- function(x){
    if(!is.na(x)){
      return(strsplit(x,", "))
    } else {
      return(NA)
    }
  }
  
 ## Get unique population-level stratifications for infection history model
  skip_pars <- c("infection_model_prior_shape1","infection_model_prior_shape2","phi")
  population_group_strats <- unique(unlist(sapply(par_tab[(par_tab$names %in% skip_pars),"stratification"],function(x) strsplit1(x))))
  population_group_strats <- population_group_strats[!is.na(population_group_strats)]
  if(length(population_group_strats) == 0) population_group_strats <- NA
  
  if(!is.null(timevarying_demographics)){
    antibody_data <- align_antibody_demographic_dat(antibody_data, timevarying_demographics)
  }
  ## Get unique demographic groups -- these correspond to different groupings for the antibody kinetics model and can be different to the population group
  ## If no demographic groups requested, and no labeling is included in antibody data, assume all individuals in the same grouping
  ## If no demographic groups supplied, set all to 1
  if(is.null(use_demographic_groups) & !("demographic_group" %in% colnames(antibody_data))){
    antibody_data$demographic_group <- 1
    if(!is.null(timevarying_demographics)){
      timevarying_demographics$demographic_group <- 1
    }
    demographics <- NULL
  } else {
    ## Otherwise, check if timevarying demographics used and update
    ## If timevarying demography is used, then set the demographic groups for each time
    if(!is.null(timevarying_demographics)){
      if("demographic_group" %in% colnames(timevarying_demographics)) use_demographic_groups <- "demographic_group"
      
      ## Get unique demographic group combinations
      demographics <- timevarying_demographics %>% 
        dplyr::select(all_of(use_demographic_groups)) %>% 
        distinct() %>% 
        dplyr::arrange(across(everything())) %>%
        dplyr::mutate(demographic_group = 1:n())
      
      ## Assign to timevarying demographics
      timevarying_demographics <- timevarying_demographics  %>% left_join(demographics,by=use_demographic_groups)
    } else {
      if("demographic_group" %in% colnames(antibody_data)) use_demographic_groups <- "demographic_group"
      ## Otherwise, set demographic_group just based on antibody_data
      demographics <- antibody_data %>% dplyr::select(all_of(use_demographic_groups)) %>% distinct() %>% dplyr::mutate(demographic_group = 1:n())
    }
    ## Merge into antibody data to get correct demographic groups at sample times
    antibody_data <- antibody_data %>% left_join(demographics,by=use_demographic_groups)
  }
  ## Now check for population group (attack rate stratifying variable)
  ## If nothing specified, set all to 1
  if(any(is.na(population_group_strats)) & !("population_group" %in% colnames(antibody_data))){
    antibody_data$population_group <- 1
    population_groups <- NULL
    if(!is.null(timevarying_demographics)){
      timevarying_demographics$population_group <- 1
    }
  } else {
    ## Otherwise, check it it's timevarying
    if(!is.null(timevarying_demographics)){
      if("population_group" %in% colnames(timevarying_demographics)) population_group_strats <- "population_group"
      ## It it is, assign unique combinations a unique population_group ID
      population_groups <- timevarying_demographics %>% 
        dplyr::select(all_of(population_group_strats))%>% 
        distinct() %>% 
        dplyr::arrange(across(everything())) %>%
        dplyr::mutate(population_group = 1:n()) %>%
        drop_na()
      ## Merge into timevarying_demographics
      timevarying_demographics <- timevarying_demographics  %>% left_join(population_groups,by=population_group_strats)
    } else {
      if("population_group" %in% colnames(antibody_data)) population_group_strats <- "population_group"
      
      ## If not timevarying, then unique combinations are just based on antibody_data
      population_groups <- antibody_data %>% 
        dplyr::select(all_of(population_group_strats))%>% 
        distinct() %>% 
        dplyr::mutate(population_group = 1:n()) %>%
        drop_na()
    }
    antibody_data <- antibody_data %>% left_join(population_groups,by=population_group_strats)
    
  }
  if(!is.null(timevarying_demographics)){
    indiv_group_indices <- timevarying_demographics %>% select(individual, time, demographic_group) %>% distinct() %>% pull(demographic_group)
    indiv_pop_group_indices <- timevarying_demographics %>% select(individual, time, population_group) %>% distinct() %>% pull(population_group)
    ## Get demographic group at birth. If isn't there, get the demographic group at the earliest time
    ## Demographic group at earliest time
    demographics_start <- timevarying_demographics %>% group_by(individual) %>% filter(time == min(time)) %>% pull(demographic_group)
    
    ## See if birth demographic groups available
    birth_demographics <- timevarying_demographics %>% group_by(individual) %>% 
      select(individual, time, birth, demographic_group) %>% 
      filter(time == birth) %>% 
      select(individual,demographic_group) %>%
      ungroup()
    demographics_start[birth_demographics$individual] <- birth_demographics$demographic_group
    indiv_group_indices <- c(rbind(demographics_start, matrix(indiv_group_indices, ncol = length(unique(antibody_data$individual)))))
  } else {
    indiv_group_indices <- antibody_data %>% select(individual, demographic_group) %>% distinct() %>% pull(demographic_group)
    indiv_pop_group_indices <- antibody_data %>% select(individual, population_group) %>% distinct() %>% pull(population_group)
  }
  
  indiv_group_indices <- indiv_group_indices - 1
  indiv_pop_group_indices <- indiv_pop_group_indices - 1
  return(list(antibody_data=antibody_data,
              timevarying_demographics=timevarying_demographics,
              demographics=demographics,
              population_groups=population_groups,
              population_group_strats=population_group_strats,
              indiv_group_indices=indiv_group_indices,
              indiv_pop_group_indices=indiv_pop_group_indices
  ))
}
