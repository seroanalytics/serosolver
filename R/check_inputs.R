#' Check parTab for simulate_data
#'
#' Checks the entries of parTab used in simulate_data
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @return nothing at the moment
#' @export
check_parTab <- function(parTab,mcmc=FALSE){
  
  ## Checks that should happen in simulate_data and run_MCMC
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Checks for waneType
  # Extract waneType
  waneType <- pars["wane_type"]
  if(is.na(waneType)) stop('wane_type is missing from parTab (specify 0 for linear decrease or 1 for piecewise linear) ') ## If user has not entered wane_type in parTab
  # Check that additional parameters are present
  if(waneType==1&(!("kappa"%in%names(pars))|!("t_change"%in%names(pars)))) stop('Parameters needed for wane_type=1 (piecewise linear) are missing')
  
  ## Additional checks that should happen in run_MCMC
  if(mcmc==TRUE){
    ## Check that all optional parameters are fixed, if not, fix them
    op_pars<-parTab[which(parTab$type==0),]
    if(all(op_pars$fixed==1)==FALSE) stop('All optional parameters must be fixed')
  }
  
}

#' Checks the entries of data used in run_MCMC
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @return nothing at the moment
#' @export
check_data <- function(data){
  ## Check that all columns are present
  col.names<-c('individual','samples','virus','titre')
  ## If there are any missing columns (NOTE: not checking if group or run are present)
  if(all(col.names%in%colnames(data))!=TRUE){
    missing.cols<-col.names[which(col.names%in%colnames(data)==FALSE)] ## Find the missing column names
    stop(paste(c('The following column(s) are missing from data: ',missing.cols),collapse = " "))
  }

}

