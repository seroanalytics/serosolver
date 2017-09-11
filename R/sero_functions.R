# - - - - - - - - - - - - - - - - - - - - -
# Main model functions
#
# Model of serological dynamics
# github.com/adamkucharski/serology-model
# - - - - - - - - - - - - - - - - - - - - -


# General functions ---------------------------------------------------------

# - - - - - - - - - - - - - - - -
# Compile c code
compile.c<-function(){
  require("Rcpp")
  setwd("./c_code")
  system("R CMD SHLIB c_model.c")
  dyn.load("c_model.so") # Note edit to remove ./ for cluster runs
  setwd("..")
}

# - - - - - - - - - - - - - - - -
# Define function to call C for expected titre 

func1 <- function(x,titredat,dd,dd2,theta,testyear_index) {
  if (!is.numeric(x)){stop("argument x must be numeric")}
  out <- .C("c_model",
            n=as.integer(length(x)), # Number potential infections
            itot=as.integer(sum(x)), # Total infections
            nsample=as.integer(length(titredat)), # Number of serological samples
            x=as.double(x), # Proposed infection history
            x1=as.double(rep(0,length(x))), # Vector for cumulative infection calc
            titre=as.double(titredat), # Titre data (not currently used)
            titrepred=as.double(rep(0,length(titredat))), # Set vector for titre prediction
            dd=as.double(dd), # Vector for long-term cross-reaction
            dd2=as.double(dd2), # Vector for short-term cross-reaction
            ntheta=as.integer(length(theta)), # number of parameters
            theta=as.double(theta), # Parameter vector
            inputtestyr=as.integer(testyear_index) # Index of serological test years
  )
  return(out$titrepred)
}

# Simulation model ---------------------------------------------------------
# This generates a simulated model dataset

simulate_data<-function(test_years,
                        historytabPost=NULL, # This imposes a particular history
                        inf_years,
                        strain_years,
                        n_part=20,
                        thetastar=theta0,
                        p.inf=0.2,
                        seedi=1,
                        antigenic.map.in=NULL,
                        pmask=NULL){
  
  # DEBUG pickyr=1; test_years=test.yr[pickyr]; historytabPost=hist.sample; thetastar=theta.max; p.inf=0.1; antigenic.map.in=NULL ; pmask=NULL
  
  # - - - -
  # DEFINE PARAMETERS
  # Fix parameters depending on what is fitted and not
  if(sum(pmask=="muShort")>0){thetastar[["muShort"]]=1e-10} # Set short term boosting ~ 0 if waning not fitted
  if(sum(pmask=="map.fit")>0){ thetastar[["sigma"]]=1} # Set cross-reactivity = 1 and don't fit if antigenic map also fitted (to avoid overparameterisation)
  if(sum(pmask=="sigma2")>0){ thetastar[["sigma2"]]=thetastar[["sigma"]] } # Fix equal if sigma same for both 
  
  part.n=c(1:n_part) # Store vector of participants
  age.yr=sample(1:80,n_part,replace = TRUE)  # Set ages
  test.n=length(test_years)  # Number of years of serological tests
  inf.n=length(inf_years) # Nuber of potential years of infection
  nstrains=length(strain_years) # Number of test strains
  sample.index=strain_years-min(inf_years)+1 # Test strain index relative to infections
  theta.sim.out=thetastar # Store input theta so can save later
  historytabSim2=historytabPost # Store input history so can save later
  
  # Check inputs are valid
  if(sum(max(test_years)==inf_years)==0){
    stop("need infection years >= test years")
    return
  }
  
  # Define antigenic map as strains on a line if no input
  if(is.null(antigenic.map.in)){antigenic.map.in=inf_years} # If no specified antigenic map, use linear function by year
  
  # Define cross-reaction matrices
  dmatrix = 1- thetastar[["sigma"]]*outputdmatrix.fromcoord(inf_years,antigenic.map.in)
  dmatrix[dmatrix<0]=0 # Fix negative values at 0
  dmatrix2 = 1- thetastar[["sigma2"]]*outputdmatrix.fromcoord(inf_years,antigenic.map.in)
  dmatrix2[dmatrix2<0]=0 # Fix negative values at 0
  
  # - - - -
  # SIMULATE ATTACK RATES
  # Set per year incidence, to create correlation between participant infection histories
  log.sd=1
  if(length(p.inf)==1){
    attack.yr=rlnorm(inf.n,meanlog=log(p.inf)-log.sd^2/2,sdlog=log.sd)
  }else{
    attack.yr=p.inf
  }
  
  # Simulate random infection history for each infection year
  if(is.null(historytabPost)){
    historytabSim=matrix(0,ncol=inf.n,nrow=n_part)
    for(jj in 1:inf.n){
      #hist0=(runif(inf.n)<attack.yr)+0
      alive=((max(test_years)-age.yr)<=inf_years[jj]) # Work out who was alive
      historytabSim[sample(part.n[alive],round(length(part.n[alive])*attack.yr[jj])),jj]=1
    }
  }else{
    historytabSim=historytabPost
  }
  
  # - - - -
  # Simulate titres for each participant
  # ii=participant | jj=test year
  
  test.list=list() # empty list to store values
  
  for(ii in 1:n_part){
    
    i.list=list() # empty list to store values
    historyii=historytabSim[ii,] # Pick simulated infection history
    
    for(jj in 1:test.n){
      
      d.ij=dmatrix[sample.index,] # Define cross-immunity matrix for sample strain
      d_vector=melt(t(d.ij))$value # Melt for use in func1()
      
      d.ij2=dmatrix2[sample.index,] # Define cross-immunity matrix for sample strain
      d_vector2=melt(t(d.ij2))$value # Melt for use in func1()
      
      testyr=test_years[jj] # Pick test year to simulate
      testyearI=c(1:inf.n)[inf_years==testyr] # Identify which strains being tested
      
      expect=func1(historyii,sample.index,d_vector,d_vector2, thetastar,testyearI) # Output expected titre
      
      # Observation model - convert to observed titres
      titredat=sapply(expect,function(x){ floor( rnorm(1,mean=x,sd=thetastar[["error"]]) ) })
      titredat[titredat<0]=0
      titredat=sapply(titredat,function(x){min(x,8)}) # Censor titres at 0
      
      # TO DO - ADD AGE  # FIXTHIS
      
      # Store outputs
      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat,
                         strain_years,
                         sample.index,
                         age.yr[ii]
                         
      )
    }
    test.list[[ii]]=i.list
  }
  test.listSim = test.list # Store simulated strain data
  
  # Export data
  if(is.null(historytabPost)){
    save(test_years,inf_years,strain_years,n_part,test.listSim,theta.sim.out, age.yr,antigenic.map.in,historytabSim,file=paste("output_simulation/Simulated_data_",seedi,".RData",sep=""))
  }else{
    save(test_years,inf_years,strain_years,n_part,test.listSim,theta.sim.out, age.yr,antigenic.map.in,historytabSim2,file=paste("output_simulation/Simulated_dataPost_",seedi,".RData",sep=""))
  }
}


# Inference model ---------------------------------------------------------

# - - - - - - - - - - - - - - - -
# Calculate antigenic distance matrix from antigenic map data

outputdmatrix.fromcoord <- function(inf_years,anti.map.in){ #anti.map.in can be vector or matrix - rows give inf_years, columns give location

    # Calculate antigenic distances
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      
      # If 1D antigenic 'line' defined, calculate distances directory from input
      (dmatrix=sapply(anti.map.in,function(x){y=abs(anti.map.in-x); y   })) 
      
    }else{ # If 2D antigenic map defined, calculate distances directory from input
      (dmatrix=apply(anti.map.in,1,function(x){y=sqrt(
        
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
        
      ); y 
      }))
    }

}




# - - - - - - - - - - - - - - - -
# Define likelihood function given expected titre and data - include uniform error term

likelihood.titre<-function(expect,titredat,theta){
  
  largett=(titredat > 8)  # Identify censored titres in data (>=8)
  smalltt=(titredat <= 0)  # Identify censored titres in data (>=8)
  
  # Define likelihood using censored normal distribution
  # First sum up titres 0 < . < 8
  p_jkMID =  ( sum( log(  pnorm(as.numeric(titredat[!largett & !smalltt])+1, mean = expect[!largett & !smalltt], sd=theta[["error"]], log = FALSE) 
                     - pnorm(as.numeric(titredat[!largett & !smalltt]), mean = expect[!largett & !smalltt], sd=theta[["error"]], log = FALSE)  ) ) )
  
  # Calculate titres >=8
  p_jkSML = sum( (  pnorm(1, mean = expect[smalltt], sd=theta[["error"]], log = T) ) )
  
  # Calculate titres = 0
  p_jkLRG = sum( (  pnorm(8, mean = expect[largett], sd=theta[["error"]], log = T, lower.tail = F) ) ) 

  # Add everything together
  p_jk = p_jkSML + p_jkMID + p_jkLRG

  p_jk

}

# - - - - - - - - - - - - - - - -
# Calculate likelihood of titres for given participant and test year

estimatelik<-function(ii,jj,historyii,dmatrix,dmatrix2,theta_star,test.list,testyearI){ 
  # ii=participant | jj=test year
  
  # jj=jj_year[kk];historyii=as.numeric(history_star[ii,]);testyearI=testyear_index[kk]
  
  # Extract relevant strain and test ID
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  

  if(length(test.jj[,1])==1){ # Check test data is available
    0
    }else{
    
    # Set up test strains
    test.part=as.numeric(test.jj[4,]) # index of sample strains data available for
    titredat=test.jj[2,] # Define titre data
    
    d.ij=dmatrix[test.part,] # Define cross-reaction matrix 1 for sample strain
    d_vector=melt(t(d.ij))$value #melt is by column
    
    d.ij2=dmatrix2[test.part,] # Define cross-reaction matrix 2 for sample strain
    d_vector2=melt(t(d.ij2))$value #melt is by column

    expect=func1(historyii,titredat,d_vector,d_vector2,theta_star,testyearI) # Output expectation

    lik = likelihood.titre(expect,titredat,theta_star) # Calculate likelihood

    lik
    
  }
  
}

# - - - - - - - - - - - - - - - -
# Set arbitrary initial condition for infection history - aim is to ensure likelihood is valid

setuphistIC<-function(ii,jj,inf.n,test.list,testyear_index, test_years, inf_years){ # ii=participant | jj=test year
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  spyear=unique(as.numeric(test.jj[3,])) # year of samples taken
  
  hist0=rep(0,inf.n)   
  #hist0[sample(c(1:inf.n),round(0.1*inf.n))]=1
  
  # Check test data available - may be issue if age column added too
  if(length(test.jj[,1])>1){
    
    # Set up test strains
    titredat=as.numeric(test.jj[2,]) # Define titre data
    
    # Use simple cutoff for titres -- set high titres = 1 in history if >=4
    for(i in 1:length(spyear)){
      if(max(titredat[(as.numeric(test.jj[3,])==spyear[i])])>=4 & runif(1)>0.2 ){
        hist0[(inf_years==spyear[i])]=1
      }
    }

  }
  
  min.range = max(1,testyear_index[1]-10) # Add an infection within past 10 years (useful for initial sampling)
  inf_index = inf_years-min(inf_years)+1
  inf_pick = sample(c(1:inf.n)[inf_index>=min.range & inf_index<= testyear_index[1]],1)  # pick strain within plausible region to add
  if(sum(hist0[inf_years < min(test_years)])==0){hist0[inf_pick]=1} # Make sure at least one infection occurs before earliest test year
  hist0
  
}




# - - - - - - - - - - - - - - - -
# Resample infection history in MH algorithm - included ageA input in case needed later

SampleHistory<-function(historyA,pick,inf.n,inf_years,age.mask){

  for(ii in pick){ # Resample subset of individuals
    
    rand1=runif(1)
    x=historyA[ii,age.mask[ii]:inf.n] # Only resample years individual was alive
    
    infvector=c(1:length(x))
    infvector2=rev(infvector)
    
    # Remove infection
    if(rand1<1/3){
      infectID=infvector[(as.numeric(x)>0)]
      if(length(infectID)>0){
        x[sample(c(infectID),1)]=0 # Why double? DEBUG
      }
    }
    
    # Add infection
    if(rand1>1/3 & rand1<2/3){
      ninfecID=infvector[(as.numeric(x)==0)]
      if(length(ninfecID)>0){
        x[sample(c(ninfecID),1)]=1
      }
    }
    
    # Move infection position
    if(rand1>2/3){
      infectID=infvector[(as.numeric(x)>0)]
      ninfecID=infvector[(as.numeric(x)==0)]
      
      if(length(infectID)>0 & length(ninfecID)>0){
        x[sample(c(infectID),1)]=0
        x[sample(c(ninfecID),1)]=1
      }
    }

    historyA[ii,age.mask[ii]:inf.n]=x # Only =1 if individual was alive
    
  } # end loop over individuals
  
  historyA
  
}

# - - - - - - - - - - - - - - - -
# Resample theta in MH algorithm

SampleTheta<-function(theta_initial,m,covartheta,covarbasic,nparam){
  
  # sample from multivariate normal distribution - no adaptive sampling
  theta_star = as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=covarbasic)))
  
  # sample from multivariate normal distribution - include adaptive samples (Roberts & Rosenthal, 2009)
  #theta_star = 0.05*as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=(2.38^2/nparam)*covarbasic))) +
  #              0.95*as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=(2.38^2/nparam)*covartheta)))
  
  names(theta_star)=names(theta_initial)
  
  # reflective boundary condition for max boost=10
  mu1=min(20-theta_star[["mu"]],theta_star[["mu"]])
  theta_star[["mu"]]=ifelse(mu1<0,theta_initial[["mu"]],mu1)
  
  #mu2=min(20-theta_star[["muShort"]],theta_star[["muShort"]])
  #theta_star[["muShort"]]=ifelse(mu2<0,theta_initial[["muShort"]],mu2)
  
  # reflective boundary condition for wane function = max is 1 for now # DEBUG
  wane2=min(2-theta_star[["wane"]],theta_star[["wane"]])
  theta_star[["wane"]]=ifelse(wane2<0,theta_initial[["wane"]],wane2)
  
  #print(rbind(theta_initial,theta_star1,theta_star2))
  return(thetaS=theta_star)
}

# - - - - - - - - - - - - - - - - 
# Acceptance probability in MH algorithm

ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  # Flat priors on theta => symmetric update probability
  calc.lik = exp(marg_likelihood_star-marg_likelihood)
  calc.lik[calc.lik>1]=1 
  calc.lik
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Metropolis-Hastings algorithm

run_mcmc<-function(
  test.yr,
  inf_years,
  strain_years,
  n_part,
  test.list,
  theta,
  runs,
  varpart_prob,
  hist.true=NULL,
  switch1=2,
  seedi=1,
  pmask=NULL,
  antigenic.map.in=NULL, # define specific map structure (or initial structure if fitting)
  age.mask = F, # prevent infection before birth
  flu_type = NULL
  ){
  
  # DEBUG SIMULATION params <<<
   # test.yr=define.year; test.list=test.listSim; theta=theta0;runs=1e2; varpart_prob=vp1;hist.true=historytabSim;  switch1=10; pmask=pmask0;
   
   # seedi=loadseed; antigenic.map.in = NULL; flu_type = flu.type;
  
  # DEBUG set params <<<
  # hist.true=NULL; test.yr=c(2009); runs=200; switch1=10; varpart_prob=0.05 ;   seedi=1; pmask=NULL ; antigenic.map.in=NULL; flu_type="H3HN"
  
  time.1 = Sys.time() # DEBUG TIME 1
  
  if(is.null(antigenic.map.in)){antigenic.map.in=inf_years} # if no input map, assume 1D
  test.n = length(test.yr); inf.n=length(inf_years); nstrains=length(strain_years) # Set up summary parameters
  sample.index = strain_years-min(inf_years)+1
  test.listPost=test.list
  #historyii=rbinom(inf.n, 1, 0.1) # DEBUG dummy infection history
  
  # Predefine index parameters to speed up code
  jj_year=c(1:test.n); testyear_index = test.yr - min(inf_years) + 1
  sample.n=length(jj_year)
  
  
  # Extract ages and create mask 
  if(age.mask==T){  # FIXTHIS
    age.list = array(unlist(test.list),dim=c(5,length(strain_years),n_part))[5,1,]
    age.mask = sapply(age.list,function(x){if(is.na(x)){1}else{match(max(min(inf_years),test_years[1]-x),inf_years)  }  })
  }else{
    age.mask = rep(1,n_part)
  }

  
  # Make adjustments depending on what is fitted and not
  if(sum(pmask=="muShort")>0){theta[["muShort"]]=1e-10} # Set short term boosting ~ 0 if waning not fitted
  #if(sum(pmask=="map.fit")>0){ theta[["sigma"]]=1; pmask=c(pmask,"sigma","sigma2")} # Set cross-reactivity = 1 and don't fit if antigenic map also fitted (to avoid overparameterisation)
  if(sum(pmask=="sigma2")>0){ theta[["sigma2"]]=theta[["sigma"]] } # Fix equal if sigma same for both 
  if(sum(pmask=="error")>0){ theta[["error"]]=1e-10 } # Fix equal if sigma same for both 
  if(sum(pmask=="tau1")>0){ theta[["tau1"]]=1e-10 } # Fix equal if sigma same for both 
  if(sum(pmask=="tau2")>0){ theta[["tau2"]]=1e-10 } # Fix equal if sigma same for both 
  
  # Preallocate memory
  nparam=length(theta); npcov=rep(1,nparam); npcov[match(pmask,names(theta))]=0 # mask specified parameters
  cov_matrix_theta0 = diag(npcov)
  cov_matrix_thetaA=cov_matrix_theta0
  
  thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta)); colnames(thetatab)=names(theta)
  thetatab[1,]=theta
  
  historytab=matrix(NA,nrow=n_part,ncol=inf.n)
  historytabCollect=historytab
  age.tab=matrix(NA,nrow=n_part,ncol=1)
  map.tab=antigenic.map.in
  
  dmatrix0 = outputdmatrix.fromcoord(inf_years,anti.map.in=map.tab) # Arrange antigenic map into cross-reaction matrix
  dmatrix20 = outputdmatrix.fromcoord(inf_years,anti.map.in=map.tab) # Arrange antigenic map into cross-reaction matrix
  
  # Pick plausible initial conditions -- using all test years
  if(is.null(hist.true)){
    for(ii in 1:n_part){
      histIC=NULL
      for(kk in 1:length(jj_year)){
        histIC=rbind(histIC,setuphistIC(ii,jj_year[kk],inf.n,test.list,testyear_index,test_years, inf_years))
      }
      histA=as.numeric(colSums(histIC)>0) # combine all histories
      histA0=histA*0
      histA0[c(age.mask[ii]:inf.n)]=histA[c(age.mask[ii]:inf.n)] # Remove history before individuals born
      
      historytab[ii,]=histA0
    }
  } else { historytab=hist.true }

  colnames(historytab)=as.character(inf_years)

  # Preallocate matrices
  likelihoodtab=matrix(-Inf,nrow=(runs+1),ncol=n_part)
  accepttabT=NULL
  accepttabH=NULL
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Run MCMC
  
  for (m in 1:runs){
    
    #time.3 = Sys.time() # DEBUG Set Time 3
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.01
      #cov_matrix_theta=epsilon0*cov_matrix_thetaA
      cov_matrix_basic=epsilon0*cov_matrix_theta0
      varpart_prob0=varpart_prob
    }else{
      epsilon0=max(0.00001,min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m)))
      #cov_matrix_theta=epsilon0*cov_matrix_thetaA
      cov_matrix_basic=epsilon0*cov_matrix_theta0
      #varpart_prob0=max(0.02,min(0.25,exp(log(varpart_prob0)+(accept_rateH-0.234)*0.999^m))) # resample max of 25%, min of 2%
    }
    
    # - - - - - - - - - - - - - - - -
    # Resample parameters

    if(m %% switch1==0 | m==1){ # m==1 condition as have to calculate all liks on first step
      theta_star = SampleTheta(thetatab[m,], m,cov_matrix_theta,cov_matrix_basic,nparam=sum(cov_matrix_theta0)) #resample theta
      
      if(sum(pmask=="sigma2")>0){ theta_star[["sigma2"]]=theta_star[["sigma"]] } # Fix equal if sigma same for both 
      
      map_star=map.tab
      history_star = historytab
      pickA=c(1:n_part)
      
    }else{
      pickA=NULL
      pickA=sample(n_part, ceiling(varpart_prob0*n_part)) # check that not length zero (i.e. at least one person sampled)
      history_star = SampleHistory(historytab,pickA,inf.n,inf_years,age.mask) #resample history
      theta_star = thetatab[m,]
    }

    # Define cross-reaction matrices with new parameters
    dmatrix =  1-theta_star[["sigma"]] *dmatrix0;  dmatrix[dmatrix<0]=0 # Arrange antigenic map into cross-reaction matrix
    dmatrix2 = 1-theta_star[["sigma2"]]*dmatrix20; dmatrix2[dmatrix2<0]=0 # Arrange antigenic map into cross-reaction matrix
    
    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Only calculate for updated history
    
    lik_val=likelihoodtab[m,]
    for(ii in pickA){
      # Set history to zero after test date
      lik.ii=rep(NA,sample.n)
      for(kk in 1:sample.n){
        #For DEBUG: set params <<<  ii=1;kk=2;historyii=as.numeric(history_star[ii,])
        lik.ii[kk]=estimatelik(ii,jj_year[kk],as.numeric(history_star[ii,]),dmatrix,dmatrix2,theta_star,test.list,testyear_index[kk])
      }
      lik_val[ii]=sum(lik.ii)
    }
    
    # - - - - - - - - - - - - - - - -
    # Metropolis Hastings step

    # History sample step
    if( (m %% switch1 != 0) & m>1){
      
      # Calculate piecewise likelihood
      output_prob = ComputeProbability(likelihoodtab[m,pickA],lik_val[pickA]) # Only calculate for selected
      pickCP = pickA[runif( length(pickA) ) < output_prob]
      
      historytab[pickCP,] = history_star[pickCP,]
      likelihoodtab[m+1,] = likelihoodtab[m,]
      likelihoodtab[m+1,pickCP] = lik_val[pickCP]
      thetatab[m+1,] = thetatab[m,]
      
    } # end history step

    # Theta sample step
    if( (m %% switch1==0) | m==1){
      
      # Estimate probability of update
      output_prob = ComputeProbability(sum(likelihoodtab[m,]),sum(lik_val)) 
      if(is.na(output_prob) & m==1){stop(paste('check initial parameter values',theta_star[["error"]]))} # check initial likelihood is valid
        
      if(runif(1) < output_prob){
        thetatab[m+1,] = theta_star
        accepttabT=c(accepttabT,1)
        likelihoodtab[m+1,] = lik_val
        
      }else{
        thetatab[m+1,] = thetatab[m,]
        likelihoodtab[m+1,] = likelihoodtab[m,]
        if(m %% switch1==0){accepttabT=c(accepttabT,0)}
      }
    } # End theta sample step

    # Calculate acceptance rate for theta
    if(m<max(100)){
      accept_rateT=0.234 # Target acceptance
    }else{
      accept_rateT = sum(accepttabT)/length(accepttabT)
    }

    # Store infection history every 20 runs
    if(m %% min(runs,20) ==0){
      historytabCollect=rbind(historytabCollect,historytab)
    }

    # Store outputs every 500 runs
    if(m %% min(runs,500) ==0){
      print(c(m,accept_rateT,varpart_prob0,round(sum(likelihoodtab[m,])))) # DEBUG HERE
      save(likelihoodtab,thetatab,inf_years,n_part,test.listPost,historytab,historytabCollect,age.tab,test.yr,switch1,file=paste("output_posterior/outputR_f",paste(test.yr,"_",collapse="",sep=""),"s",seedi,".RData",sep=""))
    }

  } #End MCMC runs loop
  
}




# Example simulation and inference ----------------------------------------
# Generates simulated data and estimates original parameters

simulation.infer <- function(seed_i,
                             mcmc.iterations=1e3, # MCMC iterations
                             flu.type="H3HN", # Flu name
                             vp1=0.2, # Proportion of infection histories to resample each step
                             fit.map = NULL, # Antigenic map input
                             fix.param="vary.init" # Add noise to input parameters?
                             ) {

  #DEBUG seed_i=1; mcmc.iterations=40; flu.type="H3HN"; fix.param ="vary.init"; vp1 =0.2; fit.map = antigenic_map
  
  loadseed=paste("SIM_",seed_i,sep="")
  
  # Define simulation parameters
  
  npartM=70 # Number of participants
  define.year=c(2007:2012) # Years serology taken
  strain_years.in = c(1968:2012) # Years of circulation of test strains included in serological analysis
  inf_years.in = c(1968:2012) # Years in which potential infections could occur
  antigenic.map0 = fit.map[,c("x_coord","y_coord")] # Define antigenic map
  pmask0=c("tau1") # Omit antigenic seniority boosting from model
  
  # Define theta
  #    short-term cross-reactivity
  thetaSim = c(mu=2, # long-term boost
               tau1=0.02, # AGS back-boost (deprecated)
               tau2=0.05, # AGS suppresion
               wane=1, # waning rate
               sigma=0.3, # long-term cross-reactivity
               muShort=2, # short-term boost
               error=1, # Measurement error
               sigma2=0.1 # short-term cross-reactivity
               ) 
  # Note: parameter ordering fixed in C code input
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # SIMULATION MODEL

  # Simulate annual attack rates based on lognormal
  sd.param = 0.5
  attack.yr = rlnorm(inf_years.in,meanlog=log(0.15)-sd.param^2/2,sdlog=sd.param) # lognormal random attack rate
  attack.yr[1] = rlnorm(1,meanlog=log(0.5)-(sd.param/2)^2/2,sdlog=(sd.param/2)) # Make first year larger (i.e. 1968 pandemic)
  write.csv(attack.yr,paste("output_simulation/sim_attackS",seed_i,".csv",sep="")) # Store outputs
  
  # Simulate titres
  
  simulate_data(test_years=define.year, # Serology test years
                inf_years=inf_years.in, # Infection years
                strain_years=strain_years.in, # Test strain years
                n_part=npartM, # N participants
                thetastar=thetaSim, # Specific parameter set
                antigenic.map.in = antigenic.map0, # Specify antigenic map to use
                #pmask=c("wane","sigma2"), # Specify what parameters to mask
                p.inf=attack.yr, # Annual attack rate
                seedi=loadseed # Seed
                )
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # INFERENCE MODEL
  # Run MCMC for above simulated data set
  
  load(paste("output_simulation/Simulated_data_",loadseed,".RData",sep="")) # Load simulation data for inference step that follows
  
  # Set initial theta -- add noise for multiple MCMC chains
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,sigma2=0.1)
  theta0[["mu"]]=2+ if(sum(fix.param=="vary.init")>0){runif(1,c(-1,1))}else{0} # basic boosting
  theta0[["tau1"]]=0.1+ if(sum(fix.param=="vary.init")>0){0.03*runif(1,c(-1,1))}else{0} # back-boost
  theta0[["tau2"]]=0.05+ if(sum(fix.param=="vary.init")>0){0.03*runif(1,c(-1,1))}else{0} # suppression via AGS
  theta0[["wane"]]= 1 + if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # -log(0.5)/1 # short term waning - half life of /X years
  theta0[["sigma"]]=0.3+ if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # long-term cross-reaction
  theta0[["sigma2"]]=0.1+ if(sum(fix.param=="vary.init")>0){0.05*runif(1,c(-1,1))}else{0} # short-term cross-reaction
  theta0[["muShort"]]=2 + if(sum(fix.param=="vary.init")>0){runif(1,c(-1,1))}else{0} # short term boosting
  theta0[["error"]]= 1 + if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # measurement error
  theta=theta0
  
  # RUN MCMC
  # Note: NEED TO RE-INITIALISE DATAFRAME IF REPEAT RUN (i.e. reload dataset above)
  run_mcmc(
    test.yr=define.year,
    inf_years,
    strain_years,
    n_part,
    theta=theta0,
    test.list=test.listSim, # use simulated data as input
    runs=mcmc.iterations, # number of MCMC runs
    varpart_prob=vp1,
    hist.true= NULL, # Set =historytabSim to start with correct infection histories
    switch1=2, # ratio of infection history resamples to theta resamples. This is fixed
    pmask=pmask0, # specify parameters to fix
    seedi=loadseed, # seed ID
    antigenic.map.in = antigenic.map0, # Define random initial map to fit
    flu_type = flu.type # note which subtype
    )
  
}



