# - - - - - - - - - - - - - - - - - - - - -
# Analysis of model outputs
#
# Model of serological dynamics
# github.com/adamkucharski/serology-model
# - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot posteriors and compare MCMC output to titre data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Convert vector to median and 95% CrI string
c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

# Convert vector to median and 95% CrI numeric
c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot MCMC posterior distributions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

plot.posteriors<-function(simDat=F,loadseed=1,flu.type="",year_test=c(2007:2012),plotmap=F,fr.lim=F,plot.corr=F,linearFn=F){

  # simDat=F;loadseed=2;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; plot.corr=F; linearFn=T
  
  # simDat=T;loadseed="SIM_1";year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; plot.corr=F; linearFn=T
  
  if(simDat==F){loadseed=paste(loadseed,"_",flu.type,sep="")}
  if(simDat==T){load(paste("output_simulation/Simulated_data_",loadseed,".RData",sep=""))}

  load(paste("output_posterior/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep=""))

  # Define lengths and sizes of inputs
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  maxlik=max(lik.tot[1:runsPOST])
  #plot(as.data.frame(thetatab)$sigma[1:runsPOST],type="l",ylab="parameter")

  # - - - - - - - 
  # Calculate ESS post burn-in
  runs1=ceiling(0.25*runsPOST)
  calculate.ESS<-function(runs1){
    thetaT=as.data.frame(thetatab)[runs1:runsPOST,]; ltheta=length(thetaT[["mu"]]); thin.theta=thetaT[seq(1,ltheta,switch1),]
    ESS.calc=effectiveSize(thin.theta); ESS.calc
  }

  ESS.calc=calculate.ESS(runs1)
  thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
  ltheta=length(thetaT[["mu"]])
  thin.theta=thetaT[seq(1,ltheta,switch1),] # Needs to be = switch1 as this is how often theta is reasmpled
  ESS.label2<-function(x){signif(as.numeric(x),3)}

  
  # Plot posterior distributions
  
  par(mfrow=c(5,2),mar = c(4,4,1,1),mgp=c(1.8,0.6,0))
  colA=rgb(0.8,0.8,0.8)
    
    param.names = names(theta)
    param.labels= names(theta)
    par(mfcol=c(length(param.names),length(param.names)),mar = c(3,3,1,1),mgp=c(1.8,0.5,0))
    
    thinner.theta=thin.theta[sample(length(thin.theta$mu),1000,replace=T),]
    thinner.theta[["muShort"]] = thinner.theta[["muShort"]] # Adjust scaling
    
    for(ii in 1:length(param.names)){
      for(jj in 1:length(param.names)){
        if(ii<jj){
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
        }
        if(ii>jj){
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.3, xlab=param.labels[ii], ylab=param.labels[jj])
        }
        if(ii==jj){
          hist(thinner.theta[[param.names[ii]]],xlab=param.labels[ii],main="")
        }
        
      }
    }

    dev.copy(png,paste("plot_outputs/CorrelationPlot_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".png",sep=""),units="cm",width=25,height=25,res=200)
    dev.off()
  
  
  # Output table of parameters
  table.param <- rbind(
    c("mu",c.text(thin.theta[["mu"]]),ESS.label2(ESS.calc[["mu"]]) ),
    c("muShort",c.text(thin.theta[["muShort"]]  ),ESS.label2(ESS.calc[["muShort"]]) ), # **Corrected for additional +1 in fitting code**
    c("sigma",c.text(thin.theta[["sigma"]]) ,ESS.label2(ESS.calc[["sigma"]]) ),
    c("sigma2",c.text(thin.theta[["sigma2"]]),ESS.label2(ESS.calc[["sigma2"]]) ),
    c("error",c.text(thin.theta[["error"]]) ,ESS.label2(ESS.calc[["error"]])),
    c("tau2",c.text(thin.theta[["tau2"]]),ESS.label2(ESS.calc[["tau2"]]) ),
    c("wane",c.text(thin.theta[["wane"]]) ,ESS.label2(ESS.calc[["wane"]])),
    c("sigmadrop",c.text(thin.theta[["mu"]]*(thin.theta[["sigma"]])),ESS.label2(ESS.calc[["sigma"]]) ),
    c("sigma2drop",c.text(thin.theta[["muShort"]]*(thin.theta[["sigma2"]])) ,ESS.label2(ESS.calc[["sigma2"]])),
    c("wane_yr",c.text(thin.theta[["muShort"]]*(thin.theta[["wane"]])),ESS.label2(ESS.calc[["wane"]]) ),
    #c("errorCorrect",c.text( pnorm(ceiling(thin.theta[["mu"]]) + 1,mean=thin.theta[["mu"]],sd=thin.theta[["error"]])-pnorm(floor(thin.theta[["mu"]]),mean=thin.theta[["mu"]],sd=thin.theta[["error"]]) ) ,ESS.label2(ESS.calc[["error"]]) )
    c("errorCorrect1",c.text( pnorm( 3,mean=1.5,sd=thin.theta[["error"]])-pnorm(2,mean=2.5,sd=thin.theta[["error"]]) ) ,ESS.label2(ESS.calc[["error"]]) ),
    c("errorCorrect2",c.text( pnorm( 4,mean=1.5,sd=thin.theta[["error"]])-pnorm(1,mean=2.5,sd=thin.theta[["error"]]) ) ,ESS.label2(ESS.calc[["error"]]) )
    
    )
  table.param = table.param %>% data.frame()
  names(table.param)=c("parameter","estimate","ESS")
  
  write.csv(table.param,paste("plot_outputs/param_table",paste(year_test,collapse="_"),"_",loadseed,".csv",sep="") )

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run multi-chain diagnostics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

plot.multi.chain.posteriors<-function(simDat=F,flu.type="H3HN",loadpick=c(1:4),burnCut=0.25,year_test=c(2007:2012),plotmap=F,fr.lim=F,linearFn=F,runsPOST=NULL){
  
  # DEBUG simDat=T;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; loadpick=c(1); burnCut=0.25; loadseed=1; linearFn=T; runsPOST=NULL
  
  storeMu = NULL
  storeMu2 = NULL
  storeLik = NULL
  storeSigma = NULL
  storeSigma2 = NULL
  storeError = NULL
  storeTau2 = NULL
  storeWane = NULL
  
  # Load simulation data
  load(paste("output_simulation/Simulated_data_SIM_",1,".RData",sep=""))
  theta.true=theta.sim.out
  
  # Load MCMC posteriors
  for(loadseed in loadpick){
  
    load(paste("output_posterior/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseedA,".RData",sep=""))
    
    # Define lengths and sizes of inputs
    lik.tot=rowSums(likelihoodtab)
    if(is.null(runsPOST)){
      runsPOST = min(8e5, length(lik.tot[lik.tot!=-Inf]) ) # Check length of posterior
    }
    maxlik=max(lik.tot[1:runsPOST])
    
    # - - - - - - - 
    # Impose burn-in
    runs1=ceiling(0*runsPOST)
    thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
    ltheta=length(thetaT[["mu"]])
    thin.theta=thetaT[seq(1,ltheta,switch1),]
    lik.totA=lik.tot[runs1:runsPOST]
    thin.lik = lik.totA[seq(1,ltheta,switch1)]
    
    storeLik = rbind(storeLik,thin.lik)
    
    storeMu = rbind(storeMu,thin.theta[["mu"]])
    storeMu2 = rbind(storeMu2,thin.theta[["muShort"]])
    storeSigma = rbind(storeSigma,thin.theta[["sigma"]])
    storeSigma2 = rbind(storeSigma2,thin.theta[["sigma2"]])
    storeWane = rbind(storeWane,thin.theta[["wane"]])
    storeError = rbind(storeError,thin.theta[["error"]])
    storeTau2 = rbind(storeTau2,thin.theta[["tau2"]])
  }
  
  # Plot MCMC chain data
  par(mfrow=c(4,2),xpd=F,mar = c(3.5,3.5,1,1),mgp=c(1.8,0.6,0))
  
  # Orange, blue, green, pink for plots
  col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  colA=rgb(0.8,0.8,0.8)
  Ctrue="black"; Wtrue=2 # True parameter colour and line
  
  maxlik=max(storeLik)
  plot(1:length(storeLik[1,]),storeLik[1,],type="l",col="white",xlab="iteration",ylab="likelihood",ylim=c(maxlik-ifelse(simDat==F,1000,500),maxlik+100) ); 
  for(ii in 1:length(loadpick)){ 
    lines(storeLik[ii,],type="l",col=col.list[[ii]] , ylim=c(maxlik-ifelse(simDat==F,1000,500),maxlik+100)) 
    }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(maxlik-ifelse(simDat==F,1000,500),maxlik+100),col="gray",lty=2)

  plot(storeWane[1,],type="l",col="white",xlab="iteration",ylab=expression(omega),ylim=c(0,1.2)); 
  for(ii in 1:length(loadpick)){ lines(storeWane[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,1.2),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["wane"]],theta.true[["wane"]]),col=Ctrue,lwd=Wtrue) }

  plot(storeMu[1,],type="l",col="white",xlab="iteration",ylab=expression(mu[1]),ylim=c(0,4)); 
  for(ii in 1:length(loadpick)){ lines(storeMu[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,4),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST),c(theta.true[["mu"]],theta.true[["mu"]]),col=Ctrue,lwd=Wtrue) }
  
  if(!(flu.type=="H3FS")){ 
    plot(storeMu2[1,],type="l",col="white",xlab="iteration",ylab=expression(mu[2]),ylim=c(0,4)); 
    for(ii in 1:length(loadpick)){ lines(storeMu2[ii,],type="l",col=col.list[[ii]]) }
    lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,4),col="gray",lty=2)
    if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["muShort"]],theta.true[["muShort"]]),col=Ctrue,lwd=Wtrue) }
  }
  
  plot(storeSigma[1,],type="l",col="white",xlab="iteration",ylab=expression(sigma[1]),ylim=c(0,0.4)); 
  for(ii in 1:length(loadpick)){ lines(storeSigma[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,0.4),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["sigma"]],theta.true[["sigma"]]),col=Ctrue,lwd=Wtrue) }
  
  if(!(flu.type=="H3FS" )){ 
    plot(storeSigma2[1,],type="l",col="white",xlab="iteration",ylab=expression(sigma[2]),ylim=c(0,0.4)); 
    for(ii in 1:length(loadpick)){ lines(storeSigma2[ii,],type="l",col=col.list[[ii]]) }
    lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,0.4),col="gray",lty=2)
    if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["sigma2"]],theta.true[["sigma2"]]),col=Ctrue,lwd=Wtrue) }
  }
  
  plot(storeError[1,],type="l",col="white",xlab="iteration",ylab=expression(epsilon),ylim=c(0,2.5)); 
  for(ii in 1:length(loadpick)){ lines(storeError[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,2.5),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["error"]],theta.true[["error"]]),col=Ctrue,lwd=Wtrue) }
  
  plot(storeTau2[1,],type="l",col="white",xlab="iteration",ylab=expression(tau),ylim=c(0,0.15)); 
  for(ii in 1:length(loadpick)){ lines(storeTau2[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,0.15),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["tau2"]],theta.true[["tau2"]]),col=Ctrue,lwd=Wtrue) }
  
  dev.copy(png,paste("plot_outputs/MCMC_chains",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".png",sep=""),units="cm",width=20,height=10,res=300)
  dev.off()
  
  
}
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot expected titres using sampled posterior estimates against true titres 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# plot.posterior.titres(loadseed="SIM",simDat=T,define.year=c(2007:2012))

plot.posterior.titres<-function(loadseed=1,year_test=c(2007:2012),flu.type,simDat=F,btstrap=5,plotRes=F){
  
  # simDat=F; flu.type="H1" ; year_test= c(2009:2011); btstrap=5; 
  # simDat=T;loadseed=1;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; plot.corr=F; linearFn=T;btstrap=5
  
  # Load simulated data
  load(paste("output_simulation/Simulated_data_",loadseed,".RData",sep=""))
  test.list=test.listSim
  hist.true=historytabSim
  theta.true=theta.sim.out
  
  load(paste("output_posterior/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep="")) # Note that this includes test.listPost
  
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  hist.switch= runsPOST/ (length(historytabCollect[,1])/n_part - 1)
  runs1=ceiling(0.2*runsPOST)

  # Set up matrices to store -- need btstrap >1
  n.strains=length(strain_years) # this loads from main_model.R
  n.test=length(test.yr)
  n.inf=length(inf_years)
  
  store.mcmc.test.data=array(NA, dim=c(btstrap,n_part,n.strains,n.test,2)) # Store expected titres for each test year
  store.mcmc.hist.data=array(NA, dim=c(btstrap,n_part,n.inf,n.test)) # Store history for each test year
  
  # - - - - - - - - - - - -
  # Sample from MCMC runs to get stored matrices of expected titre and estimated infection years

  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1) # Pick random history
    
    pickAhist=ceiling(pickA/hist.switch)+1 # Extract from history table 
    
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
    
  for(pickyr in 1:n.test){
    
    # Note here that inf_years and strain_years are loads from main_model.R
    # Output expected titre 
    simulate_data(test.yr[pickyr],historytabPost=hist.sample,
                  inf_years,
                  strain_years,
                  n_part,thetastar=theta.max,p.inf=0.1,
                  antigenic.map.in=antigenic_map)
    
    load("output_simulation/Simulated_dataPost_1.RData") # Load simulated data
    
    # Compile data for this participant
    for(ii0 in 1:n_part){
      
      sim.titre=test.listSim[[ii0]][[1]] # sort sample years - drawn from simulated data above
      hist.sampleB=hist.sample;  hist.sampleB[,as.numeric(colnames(hist.sample))>test.yr[pickyr]]=0 # don't show infections after test year
      store.mcmc.test.data[sampk,ii0,,pickyr,1]= min(inf_years)-1+sort(sim.titre["sample.index",]) # Sampled strain years
      s.titre=sim.titre["titredat",order(sim.titre["sample.index",])]
      store.mcmc.test.data[sampk,ii0,,pickyr,2]=s.titre # Sampled expected titre
      #store.mcmc.test.data[sampk,ii0,,pickyr,2]=rpois(length(s.titre),lambda=s.titre) # Sampled expected titre - include Poisson noise
      store.mcmc.hist.data[sampk,ii0,,pickyr]=hist.sampleB[ii0,] # Sampled history
      
    }  # end loop over participants
    
  } # end loop over test years
  } # end bootstrap loop
  
  # - - - - - - - - - - - - - - - - - - - 
  # Plot histogram of estimated vs true titre
  
  if(plotRes==T){ # Also plot residuals if specified
  
    compTab = NULL
    likTab = NULL
    compTabNULL = NULL
    for(pickyr in 1:n.test){
      for(ii0 in 1:n_part){ # Check more than one strain per year
        if(!is.null(dim(store.mcmc.test.data[,ii0,,pickyr,2]))){
          estT = apply(store.mcmc.test.data[,ii0,,pickyr,2],2,median)
        }else{
          estT = sapply(store.mcmc.test.data[,ii0,,pickyr,2],median)
        }
        truT = test.list[[ii0]][[pickyr]][2,]
        
        # Calculate likelihood -- if entries = NA it's because strains are missing
        
        if(!is.na(truT[1])){likTab=c(likTab,likelihood.titre(estT,truT,theta.max))}
        
        # store actual error and MC poisson error
        compTab = c(compTab, floor(estT) - truT)
        compTabNULL = c(compTabNULL, estT - sapply(estT,function(x){rpois(1,lambda=x)})  )
      }
    }
    
    compTab = as.numeric(compTab); compTab = compTab[!is.na(compTab)]; compTabNULL = compTabNULL[!is.na(compTab)]
    
    par(mfrow=c(1,1),mar = c(5,5,1,1))
    
    breaks0=seq(round(min(c(compTab,compTabNULL)))-0.5,round(max(c(compTab,compTabNULL)))+1.5,1)
    hist(compTab , breaks = breaks0 , freq=F, col =rgb(0.8,0.8,0.8),main=NULL,xlab="expected titre - observed titre"); #title(LETTERS[1],adj=0)
    #hist(compTabNULL , breaks = breaks0 , freq=F, col =rgb(0.8,0.8,0.8),main=NULL,xlab="expected titre - Pois(expected titre)"); #title(LETTERS[2],adj=0)
    
    stats.store = rbind( c(sum( abs(compTab)<=1 )/length(compTab), sum( abs(compTab)<=2 )/length(compTab), sum( abs(compTab)<=3 )/length(compTab) ), 
                         c(sum( abs(compTabNULL)<=1 )/length(compTabNULL), sum( abs(compTabNULL)<=2 )/length(compTabNULL), sum( abs(compTabNULL)<=3 )/length(compTabNULL)  )
                  )
    colnames(stats.store)=c("within 1","within 2","within 3")
    rownames(stats.store)=c("model","null")
    write.csv(stats.store,paste("plot_outputs/titre_compare/Store_residuals_",loadseed,".csv",sep=""))
    
    dev.copy(pdf,paste("plot_outputs/titre_compare/Residuals_",loadseed,"_lin",linearFn,".pdf",sep=""),width=5,height=4)
    dev.off()

  
  }
  
  # - - - - - - - - - - - - - - - - - - - 
  # Plot figures from MCMC posteriors

    colN=if(n.test>1){n.test}else{5}
    loopN=if(n.test>1){5}else{25}
    par(mfrow=c(5,colN),mar = c(2,2,1.5,0.2),mgp=c(2,0.5,0))

    yrange1=c(-0.2,9)
    
  for(ii0 in 1:n_part){
      
    for(pickyr in 1:n.test){
        
      simtitreX=store.mcmc.test.data[,ii0,,pickyr,1] # pick out years
      simtitreY=store.mcmc.test.data[,ii0,,pickyr,2] # pick out estimates
      hist.sample=store.mcmc.hist.data[,ii0,,pickyr] # for participant ii0 in year pickyr
  
      plot(inf_years,8*hist.sample[1,],type="l",ylim=yrange1,yaxs='i',col='white',xlab="",ylab="",main=paste("Participant ",ii0,", ",year_test[pickyr],sep=""))
      
      # Sample from infection history
      for(jj in 1:n.inf){
          lines(min(inf_years)-1+c(jj,jj),c(yrange1[1]+0.05,yrange1[2]-0.05),col=rgb(0,0,0,sum(hist.sample[,jj])/btstrap),lwd=2) # Plot estimated infections
      }
      
      # Calculate credible interval for expected titres
      if(!is.null(dim(simtitreY))){ # Check more than one strain per year
        medP=apply(simtitreY,2,function(x){median(x)})
        ciP1=apply(simtitreY,2,function(x){quantile(x,0.025)})
        ciP2=apply(simtitreY,2,function(x){quantile(x,0.975)})
        polygon(c(simtitreX[1,],rev(simtitreX[1,])),c(ciP1,rev(ciP2)),lty=0,col=rgb(0,0.3,1,0.2))
        lines(simtitreX[1,],medP,pch=1,col='blue')
        points(simtitreX[1,],medP,pch=19,cex=0.5,col='blue')
      }else{
        medP=median(simtitreY)
        ciP1=quantile(simtitreY,0.025)
        ciP2=quantile(simtitreY,0.975)
        #polygon(c(simtitreX,rev(simtitreX)),c(ciP1,rev(ciP2)),lty=0,col=rgb(0,0.3,1,0.2))
        lines(c(simtitreX[1],simtitreX[1]),c(ciP1,ciP2),pch=1,col='blue')
        points(simtitreX[1],medP,pch=19,cex=0.5,col='blue')
      }
      
      # Plot true titres
      points(min(inf_years)-1+test.list[[ii0]][[pickyr]][4,],test.list[[ii0]][[pickyr]][2,],pch=1,col='red')
      
      # Plot true infections if simulation
      if(simDat==T){
        histSim1=hist.true[ii0,]
        histSim1[inf_years>year_test[pickyr]]=0 # don't show infections after test year
        lenhis=rep(0,length(histSim1))
        for(jj in 1:length(lenhis)){
          lines(min(inf_years)-1+c(jj,jj),c(20*histSim1[jj]-11.5,20*histSim1[jj]-11),col=rgb(0,0.6,0),lwd=3)
        }
      }
      
      if(ii0 %% loopN==0 | ii0==n_part){
        dev.copy(png,paste("plot_outputs/titre_compare/sim",loadseed,"_",ii0,"P.png",sep=""),units="cm",width=25,height=15,res=150)
        dev.off()
      }
    
    }  # end loop over participants
      
  } # end loop over test years

}
