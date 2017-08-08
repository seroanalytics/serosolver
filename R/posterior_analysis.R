# - - - - - - - - - - - - - - - - - - - - -
# Model of serological dynamics
# github.com/adamkucharski/serology-model
#
# Analysis of model outputs
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

# Transform antigenic maps to be uniform
scale.map<-function(map.pick){
  f.m=1 #length(map.pick[,1])
  # translate map to finish at (0,0)
  map.pick[,1]=map.pick[,1]-map.pick[f.m,1]
  map.pick[,2]=map.pick[,2]-map.pick[f.m,2]
  
  # rotate map so final coordinate is (0,0) and penultimate is (0,1) 
  #r.theta=atan(map.pick[(f.m-1),1]/map.pick[(f.m-1),2]) # angle of rotation
  #map.pick[,1]=map.pick[,1]*cos(r.theta)-map.pick[,2]*sin(r.theta)
  #map.pick[,2]=map.pick[,1]*sin(r.theta)+map.pick[,2]*cos(r.theta)
  
  # flip so 3rd from last is positive
  #map.pick[,1]=sign(map.pick[(f.m-2),1])*map.pick[,1]
  
  map.pick
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot MCMC posterior distributions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

plot.posteriors<-function(simDat=F,loadseed=1,flu.type="",year_test=c(2007:2012),plotmap=F,fr.lim=F,plot.corr=F,linearFn=F){

  # simDat=F;loadseed=2;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; plot.corr=F; linearFn=T
  
  # simDat=T;loadseed="SIM_1";year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; plot.corr=F; linearFn=T
  
  if(simDat==F){loadseed=paste(loadseed,"_",flu.type,sep="")}
  if(simDat==T){load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))}

  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,"_lin",linearFn,".RData",sep=""))
  if(flu.type=="H3HN" & simDat==F){
    yob.data=data.frame(read.csv("datasets/HaNam_YOB.csv",header=FALSE)) # Import age distribution 
    n.alive=sapply(inf_years,function(x){sum(yob.data<=x)})}
  
  if(flu.type=="H3FS" & simDat==F){
    yob.data=data.frame(read.csv("datasets/FluScape_YOB.csv",header=FALSE)) # Import age distribution
    n.alive=sapply(inf_years,function(x){sum(yob.data<=x)})
  }
  
  if(simDat==T){
    yob.data=cbind(rep(1,n_part),rep(1,n_part)) # Import age distribution
    n.alive=n_part+0*inf_years
  }
  
  if(simDat==T){ par(mfrow=c(5,2)) }else{ if(flu.type=="H3FS"){par(mfrow=c(3,2))}else{par(mfrow=c(4,2))} }
  par(mar = c(4,4,1,1))
  par(mgp=c(1.8,0.6,0))
  colA=rgb(0.8,0.8,0.8)
  
  # Define lengths and sizes of inputs
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  maxlik=max(lik.tot[1:runsPOST])
  #plot(as.data.frame(thetatab)$sigma[1:runsPOST],type="l",ylab="parameter")

  # - - - - - - - 
  # Calculate ESS by burn-in
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
  ESS.label<-function(x){""}
  # - - - - - - - - - - - - - - -
  # Plot results
  #if(simDat==T){ plot(c(1:runsPOST),lik.tot[1:runsPOST],type="l",ylab="likelihood",xlab="iteration") }
  
  #plot(c(runs1:runsPOST),lik.tot[runs1:runsPOST],type="l",ylab="likelihood",xlab="iteration")
  plot(c(1,2,3))
  
  # Plot histograms of parameters
  #hist(thin.theta[["error"]],main=ESS.label(ESS.calc[["error"]]),col=colA,xlab="error",prob=TRUE,xlim=c(0,ifelse(fr.lim==F,0.1,1.1*max(thin.theta[["error"]])))); if(simDat==T){abline(v=theta.sim.out[["error"]],col="red")}
  hist(thin.theta[["mu"]],main= ESS.label(ESS.calc[["mu"]]),col=colA,xlab="mu",prob=TRUE,xlim=c(0,ifelse(fr.lim==F,4,1.1*max(thin.theta[["mu"]])))); if(simDat==T){abline(v=theta.sim.out[["mu"]],col="red")}
  hist(thin.theta[["sigma"]],main=ESS.label(ESS.calc[["sigma"]]),col=colA,xlab="sigma",xlim=c(0,ifelse(fr.lim==F,0.5,1.1*max(thin.theta[["sigma"]])))); if(simDat==T){abline(v=theta.sim.out[["sigma"]],col="red")}
  if(flu.type=="H3HN"){hist(thin.theta[["muShort"]] ,main=ESS.label(ESS.calc[["muShort"]]),col=colA,xlab="mu_Short",prob=TRUE,xlim=c(0,ifelse(fr.lim==F,4,1.1*max(thin.theta[["muShort"]])))); if(simDat==T){abline(v=theta.sim.out[["muShort"]],col="red")}
  }
  #**Corrected for additional +1 in fitting code**
  if(flu.type=="H3HN"){hist(thin.theta[["sigma2"]],main=ESS.label(ESS.calc[["sigma2"]]),col=colA,xlab="sigma2",xlim=c(0,ifelse(fr.lim==F,0.5,1.1*max(thin.theta[["sigma2"]])))); if(simDat==T){abline(v=theta.sim.out[["sigma2"]],col="red")}
  }
  hist(thin.theta[["error"]],main=ESS.label(ESS.calc[["error"]]),col=colA,xlab="error",prob=TRUE,xlim=c(0,ifelse(fr.lim==F,0.15,1.1*max(thin.theta[["error"]])))); if(simDat==T){abline(v=theta.sim.out[["error"]],col="red")}
  hist(thin.theta[["tau2"]],main=ESS.label(ESS.calc[["tau2"]]),col=colA,xlab="tau2",prob=TRUE,xlim=c(0,ifelse(fr.lim==F,0.15,1.1*max(thin.theta[["tau2"]])))); if(simDat==T){abline(v=theta.sim.out[["tau2"]],col="red")}
  if(flu.type=="H3HN"){hist(thin.theta[["wane"]],main=ESS.label(ESS.calc[["wane"]]),col=colA,xlab="wane",prob=TRUE,xlim=c(0,ifelse(fr.lim==F,3,1.1*max(thin.theta[["wane"]])))); if(simDat==T){abline(v=theta.sim.out[["wane"]],col="red")}
  }
  # Plot distribution of infections
  hist.sample=length(historytabCollect[,1])/n_part # need this sample value because table is stacked
  
  ind.infN=rowSums(historytabCollect[((round(0.2*hist.sample)*n_part)+1):(hist.sample*n_part),])
  
  #hist(ind.infN,breaks=seq(0,max(ind.infN)+1,2),col=colA,xlab="infections",prob=TRUE,main="",xlim=c(0,44)) #paste("mean/med=",signif(mean(ind.infN),2),"/",median(ind.infN),sep="")

  # Adjust for age - TO ADD?
  
  #age.data=sort(max(test_years)-yob.data[,1])
  #years.at.risk=(max(test_years)-min(inf_years))
  #age.data[age.data > years.at.risk] =years.at.risk

  dev.copy(pdf,paste("plot_simulations/posterior",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,"_lin",linearFn,".pdf",sep=""),width=8,height=12)
  dev.off()
  
  
  if(plot.corr==T){
    
    param.names = c("mu","muShort","sigma","sigma2","error","tau2","wane")
    param.labels= c(expression(mu[1]),expression(mu[2]),expression(sigma[1]),expression(sigma[2]),expression(epsilon),expression(tau),expression(omega))
    par(mfcol=c(length(param.names),length(param.names)))
    par(mar = c(3,3,1,1))
    par(mgp=c(1.8,0.5,0))
    
    thinner.theta=thin.theta[sample(length(thin.theta$mu),1000,replace=T),]
    thinner.theta[["muShort"]] = thinner.theta[["muShort"]] # Adjust scaling
    
    for(ii in 1:length(param.names)){
      for(jj in 1:length(param.names)){
        if(ii<=jj){
          
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
        }else{
          
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.3, xlab=param.labels[ii], ylab=param.labels[jj])
        }
        
        
      }
    }

    dev.copy(png,paste("plot_simulations/CorrelationPlot_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".png",sep=""),units="cm",width=25,height=25,res=200)
    dev.off()
  }
  
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
  
  write.csv(table.param,paste("plot_simulations/param_table",paste(year_test,collapse="_"),"_",loadseed,".csv",sep="") )
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot attack rates, scaled by proportion alive for SIMULATED DATA
  
  if(simDat==T){ #Add simulated attack rates and compare distribution
    yob.data=cbind(rep(1,n_part),rep(1,n_part)) # Import age distribution
    n.alive=n_part+0*inf_years
  
  attack=colSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])/(length(ind.infN)*(n.alive/length(yob.data[,1]))) #scale by proportion alive
  attackCI=NULL
  for(jj in 1:length(inf_years)){
    htest <- binom.test(round(n.alive*attack)[jj], n.alive[jj], p = 1,conf.level=0.95)
    meanA=attack[jj]
    conf1=htest$conf.int[1]
    conf2=htest$conf.int[2]
    attackCI=rbind(attackCI,c(meanA,conf1,conf2))
  }
  attackCI=data.frame(attackCI)
  names(attackCI)=c("mean","CI1","CI2")
  colB=rgb(0,0,0.8)
  colG=rgb(0.2,0.2,1)
  plot(inf_years,attackCI$mean,pch=19,col=colG,ylim=c(0,1),xlab="year",ylab="attack rate")
  for(kk in 1:length(inf_years)){ # Iterate across test years
    if(sum(kk==match(year_test,inf_years))>0){colComp="red"}else{colComp=colG}
    points(inf_years[kk],attackCI$mean[kk],pch=19,col=colComp)
    lines(c(inf_years[kk],inf_years[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=colComp)
  }
  
  load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))
  attack.yr=colSums(historytabSim)/n_part
  
  points(inf_years,attack.yr,col="black",lwd=2)
  
  dev.copy(pdf,paste("plot_simulations/postPlots",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=8,height=12)
  dev.off()
  
  #par(mfrow=c(1,1))

  # Plot comparison
    plot(attack.yr,attackCI$mean,pch=19,col=colG,xlim=c(0,0.6),ylim=c(0,0.6),xlab="true attack rate",ylab="estimated attack rate", xaxs="i", yaxs="i")
    lines(c(0,1),c(0,1),col='grey')
    for(kk in 1:length(inf_years)){ # Iterate across test years
      if(sum(kk==match(year_test,inf_years))>0){colComp="red"}else{colComp=colG}
      points(attack.yr[kk],attackCI$mean[kk],pch=19,col=colComp)
      lines(c(attack.yr[kk],attack.yr[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=colComp)
    }
    
    dev.copy(pdf,paste("plot_simulations/postCompare",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=8,height=12)
    dev.off()
    
    # - - - - - - - - - - - - - - - - - - 
    # Calculate and plot four fold rise in data
    sconverttab = NULL
    
    for(kk in 2:(length(test.yr)-1) ){ # Only valid for 2008-2011 (no test strains for 2012)
      pyear2=0
      pyear4=0
      nyear=0
      for(ii in 1:n_part){
        t.part1=test.listSim[[ii]][[kk-1]]
        t.part2=test.listSim[[ii]][[kk]]

          # Check to match test strains
          matchd1d2 = t.part2[3,]==test.yr[kk]
          
          if(length(matchd1d2) > 0){
            
            diffT = t.part2[2,matchd1d2] - t.part1[2,matchd1d2] # Compare titres
            nyear = nyear +1
            if(median(diffT) >= 2){pyear4 = pyear4 + 1}
            if(median(diffT) >= 1){pyear2 = pyear2 + 1}
          }
          
      }
      sconverttab=rbind(sconverttab, c(pyear4/nyear,pyear2/nyear))
    }
    
    pick_r=match(test_years[2:(length(test.yr)-1)],inf_years)
    
    par(mfrow=c(1,1),mar = c(4,4,1,1),mgp=c(1.8,0.6,0))
    
    plot(attack.yr[pick_r],attackCI[pick_r,"mean"],pch=19,col=colG,xlim=c(0,0.6),ylim=c(0,1),xlab="true attack rate",ylab="estimated attack rate", xaxs="i", yaxs="i")
    lines(c(0,1),c(0,1),col='grey')
    for(kk in pick_r){ # Iterate across test years
      if(sum(kk==match(year_test,inf_years))>0){colComp="red"}else{colComp=colG}
      points(attack.yr[kk],attackCI$mean[kk],pch=19,col=colComp,cex=1)
      lines(c(attack.yr[kk],attack.yr[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=colComp)
    }
    
    points(attack.yr[pick_r], sconverttab[,2],pch=1,cex=1.2,col="black")
    points(attack.yr[pick_r], sconverttab[,1],pch=19,cex=1.2,col="black")

    
    dev.copy(pdf,paste("plot_simulations/True_vs_rise",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=8,height=12)
    dev.off()
    
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot attack rates, scaled by proportion alive for TRUE DATA

  par(mfrow=c(1,1))
  par(mar = c(5,5,1,1))
    
  attack = colSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),]) /( length(ind.infN) * ( n.alive/length(yob.data[,1]) ) ) #scale by proportion alive
  attack = sapply(attack,function(x){min(x,1)}) #adjust for short debug runs, as may have large attack rate
  
  attackCI=NULL
  for(jj in 1:length(inf_years)){
    htest <- binom.test(round(n.alive*attack)[jj], n.alive[jj], p = 1,conf.level=0.95)
    meanA=attack[jj]
    conf1=htest$conf.int[1]
    conf2=htest$conf.int[2]
    attackCI=rbind(attackCI,c(meanA,conf1,conf2))
  }
  attackCI=data.frame(attackCI)
  names(attackCI)=c("mean","CI1","CI2")
  colB=rgb(0,0,0.8)
  plot(inf_years,attackCI$mean,pch=19,col=colB,ylim=c(0,1),xlab="year",ylab="attack rate")
  for(kk in 1:length(inf_years)){ # Iterate across test years
     lines(c(inf_years[kk],inf_years[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=colB)
  }
  
  #dev.copy(pdf,paste("plot_simulations/attackSimple",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=6,height=6)
  #dev.off()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot estimated attack rates against virus isolates that year
  
  if(flu.type=="H3HN" & length(year_test)==6 & simDat==F){
    
    # Calculate number of isolates
    collect.list = data.frame(rbind(c(2006,"2006-12-31"),c(2007,"2007-12-31"),c(2008,"2008-12-31"),
    c(2009,"2009-06-30"),c(2010,"2010-04-30"), c(2011,"2011-07-31"),c(2012,"2012-05-31") ),stringsAsFactors = F) # Dat of isolation
    names(collect.list)=c("year","sample")
    collect.list$year=as.numeric(collect.list$year)
    collect.list$sample=as.Date(collect.list$sample)
    
    flu.isolates=data.frame(read.csv("datasets/Vietnam_H3.csv",stringsAsFactors = F)) # Data from http://www.who.int/influenza/gisrs_laboratory/flunet/en/
    flu.isolates$Start_date=as.Date(flu.isolates$Start_date)
    flu.isolates$A_H3[is.na(flu.isolates$A_H3)]=0 # Set blank spaces to zero
    
    # Count samples within region
    isolatetab=NULL
    for(ii in 2:length(test.yr) ){
      isolatetab=c(isolatetab,
                       sum(flu.isolates[flu.isolates$Start_date > collect.list[ii,"sample"] & flu.isolates$Start_date < collect.list[ii+1,"sample"] ,"A_H3"])
                       )
    }
    pick_years = match(test.yr[2:6],inf_years)
    
    # Calculate values based on FOUR-fold rise or more in that year
    sconverttab=NULL
    
    load("R_datasets/HaNam_data.RData")
    
    for(kk in 2:length(test.yr) ){ # Only valid for 2008-2011 (no test strains for 2012)
      pyear2=0
      pyear4=0
      nyear=0
      for(ii in 1:n_part){
        t.part1=test.list[[ii]][[kk-1]]
        t.part2=test.list[[ii]][[kk]]
        
        if(length(t.part1[,1])>1 & length(t.part2[,1])>1){ # Check if year to compare
          # Check to match test strains
          matchd1d2 = intersect(names(t.part1),names(t.part2)[t.part2[3,]==test.yr[kk]])
          
          if(length(matchd1d2) > 0){
            
            diffT = (t.part2[2,match(matchd1d2,names(t.part2))] - t.part1[2,match(matchd1d2,names(t.part1))] ) %>% as.numeric() # Compare titres
            nyear = nyear +1
            #if(max(diffT) >= 2){pyear4 = pyear4 + 1}; if(max(diffT) >= 1){pyear2 = pyear2 + 1}
            if(median(diffT) >= 2){pyear4 = pyear4 + 1}; if(median(diffT) >= 1){pyear2 = pyear2 + 1}
          }
          
        }
      }
      sconverttab=rbind(sconverttab, c(pyear4/nyear,pyear2/nyear))
    }
    
  #}
    
    #attackCIsero=NULL
    #for(jj in 1:5){
    #  if(jj == 5){attackCIsero=rbind(attackCIsero,c(-1,-1,-1))}else{ # As NA in final year
    #    htest <- binom.test(round(sconverttab[jj,1]*sconverttab[jj,2]), sconverttab[jj,2], p = 1,conf.level=0.95)
    #    meanA=sconverttab[jj,1]
    #    conf1=htest$conf.int[1]
    #    conf2=htest$conf.int[2]
    #    attackCIsero=rbind(attackCIsero,c(meanA,conf1,conf2))
    #  }
    #}
    #attackCIsero=data.frame(attackCIsero)
    #names(attackCIsero)=c("mean","CI1","CI2")
    
    # fit distribution to infection data ***
    
    #Poisson.distn = fitdistr(ind.infN,densfun="Poisson")
    #NBin.distn = fitdistr(ind.infN,densfun="negative binomial")
    
    inf.individuals = matrix(ind.infN,ncol=69,byrow=T)
    inf.samples = length(inf.individuals[,1])
    
    inf.individuals.sorted = inf.individuals[,order(colSums(inf.individuals))]
    age.list.censored = 2012 - sapply(yob.data[,1],function(x){max(x,1968)})
    
    per.year=apply(inf.individuals.sorted,2,function(x){c.nume(x)})
    #plot(per.year[1,]/age.list.censored)
    
    # Plot data
    # Replot attack rates
    
    par(mfrow=c(1,3))
    par(mar = c(3,3,1,1))
    par(mgp=c(1.8,0.6,0))
    if(flu.type=="H3HN"){colA0=c(rep(colB,length(inf_years)-5),rep("red",5))}else{colA0=rep(colB,length(inf_years))}
    letN <- 0
    
    if(flu.type=="H3HN"){
      letN <- 0
      #plot(c(1,2))
    }
    
    plot(inf_years,attackCI$mean,pch=19,col=colA0,ylim=c(0,1),xlab="year",ylab="estimated attack rate", yaxs="i")
    for(kk in 1:length(inf_years)){ # Iterate across test years
      lines(c(inf_years[kk],inf_years[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=colA0[kk])
    }
    title(main=LETTERS[1+letN],adj=0)
    
    #if(flu.type=="H3HN"){
    
    
     plot(isolatetab,attackCI$mean[pick_years],pch=19,cex=1.2,col="red",ylim=c(0,1),xlim=c(0,700),xlab="H3 isolates during sample period",ylab="estimated attack rate", xaxs="i", yaxs="i")
      for(kk in 1:length(pick_years)){ # Iterate across test years
        lines(c(isolatetab[kk],isolatetab[kk])+0.5,c(attackCI$CI1[pick_years[kk]],attackCI$CI2[pick_years[kk]]),col="red")
      }
  
      points(isolatetab,sconverttab[,2],pch=1,cex=1.2,col="black")
      points(isolatetab,sconverttab[,1],pch=19,cex=1.2,col="black")
      #for(kk in 1:length(pick_years)){ # Iterate across test years
      #  lines(c(isolatetab[kk],isolatetab[kk])+0.5,c(attackCIsero$CI1[kk],attackCIsero$CI2[kk]),col=rgb(0.5,0.5,0.5))
      #}
      
      print(cor.test(isolatetab,attackCI$mean[pick_years]))
      print(cor.test(isolatetab,sconverttab[,2]))
      print(cor.test(isolatetab,sconverttab[,1]))
      
      title(main=LETTERS[2+letN],adj=0)
  
      # PLOT INFECTIONS
      
      
      plot(per.year[1,],pch=19,cex=0.9,col=rgb(0.25,0.25,0.25),ylim=c(0,45),xlim=c(0,n_part+1),xlab="participant",ylab="estimated number of infections", xaxs="i", yaxs="i")
      for(kk in 1:n_part){ # Iterate across test years
        lines(c(kk,kk),c(per.year[2,kk],per.year[3,kk]),col=rgb(0.25,0.25,0.25))
      }
      
      ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
      #hist(ind.infN,breaks=seq(0,max(ind.infN)+1,2),col =rgb(0.8,0.8,0.8),xlab="number of infections",prob=TRUE,xlim=c(0,50),ylim=c(0,0.1),xaxs="i",yaxs="i", main=NULL,ylab="density")
      #abline(v=median(ind.infN),col="red",lty=2,lwd=1.5); abline(v=mean(ind.infN),col="red",lwd=1.5)
      title(main=LETTERS[3+letN],adj=0)
    
    } # END OF LOOP
    
   dev.copy(pdf,paste("plot_simulations/attack",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=8,height=3)
    #dev.copy(pdf,paste("plot_simulations/attack",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=8,height=6,useDingbats = F)
    dev.off()
    
    

  # - - - - - - - - - - - -
  # Plot antigenic map
  if(plotmap==T){
    load("datasets/spline_fn.RData") # load spline function for map
    ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
    
    par(mfrow=c(1,1))
    par(mar = c(4,4,1,1))
    
    vals1=predict(am.spl,scalemap(inf_years,inf_years))
    map.sample=length(map.tabCollect)

    MTx=c(332,372)
    MTy=c(245,261)
    plot(vals1,type="l",xlim=MTx,ylim=MTy,col="white",lwd=2,xlab="strain dimension 1", xaxs="i", yaxs="i",
         ylab="strain dimension 2")
    
    points(ag.coord$AG_y,ag.coord$AG_x,col='grey',pch=19,cex=0.8)
    lines(vals1,col="blue")
    points(vals1$x,vals1$y,col="blue",pch=19,cex=0.5)
    points(vals1$x[1],vals1$y[1],col="blue",pch=19)
    points(vals1$x[45],vals1$y[45],col="blue",pch=19)
    text(x=335.9,y=255.2,labels="1968",col="blue")
    text(x=370.2,y=248.4,labels="2012",col="blue")
    
    dev.copy(pdf,paste("plot_simulations/antigenic_map",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=6,height=4)
    dev.off()
    
  }

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compare multiple MCMC outputs for vector of different cross-sectional years fitted
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

plot.compare<-function(simDat=F,loadseed=1,define.year.vec=c(2007:2012)){
  
  loadseed="1" #"1_w12"
  n.test=length(define.year.vec)
  names1=c("test","mu","tau1","tau2","disp_k","sigma","muShort","error","infections")
  store.val=array(NA,dim=c(3,length(names1),n.test),dimnames=list(NULL,names1,NULL))
  
  labels.Y=if(length(define.year.vec)>6){c(define.year.vec[1:6],"All")}else{define.year.vec}
  
  for(kk in 1:n.test){
  
    load(paste("posterior_sero_runs/outputR_f",define.year.vec[kk],"_s",loadseed,".RData",sep=""))
    # Store median and 95% CrI
    lik.tot=rowSums(likelihoodtab); maxlik=max(lik.tot); runsPOST=length(lik.tot[lik.tot!=-Inf]); runs1=ceiling(0.25*runsPOST)
    hist.sample=length(historytabCollect[,1])/n_part; ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
    store.val[,,kk]=cbind(rep(kk,3),
                     c.nume(as.data.frame(thetatab)$mu[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$tau1[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$tau2[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$disp_k[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$sigma[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$muShort[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$error[runs1:runsPOST]),
                     c.nume(ind.infN))
  }
  # - - - - - - - - - - - - - - - 
  # Plot comparison of parameters
  par(mfrow=c(2,4))
  
  range.p=rbind(c(0,0),c(0,5),c(0,0.5),c(0,0.5),c(0,5),c(0,0.5),c(0,10),c(0,0.1),c(0,30))# define parameter ranges for plots
  range.p[5,]=c(0,max(store.val[1,5,]))
  for(jj in 2:length(names1)){   # Iterate across parameters
    colA=rgb(0,0,0.8)
    plot(c(1:n.test),c(1:n.test),pch=19,col=rgb(1,1,1),ylim=range.p[jj,],xaxt="n",xlab="test year",ylab="estimate",main=names1[jj])
    axis(1,at=c(1:n.test),labels=labels.Y)
    grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  
      for(kk in 1:n.test){ # Iterate across test years
        points(kk,store.val[1,jj,kk],pch=19,col=colA)
        lines(c(kk,kk),c(store.val[2,jj,kk],store.val[3,jj,kk]),col=colA)
      }
  }

  dev.copy(pdf,paste("plot_simulations/posterior_compare",paste(define.year.vec,"_",collapse="",sep=""),".pdf",sep=""),width=12,height=8)
  dev.off()
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run multi-chain diagnostics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


plot.multi.chain.posteriors<-function(simDat=F,flu.type="H3HN",loadpick=c(1:4),burnCut=0.25,year_test=c(2007:2012),plotmap=F,fr.lim=F,linearFn=F,runsPOST=NULL){
  
  # simDat=F;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; loadpick=c(1:4); burnCut=0.25; linearFn=T
  # simDat=T;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; loadpick=c(1:4); burnCut=0.25; loadseed=1; linearFn=T; runsPOST=NULL
  
  storeMu = NULL
  storeMu2 = NULL
  storeLik = NULL
  storeSigma = NULL
  storeSigma2 = NULL
  storeError = NULL
  storeTau2 = NULL
  storeWane = NULL
  
  if(simDat==T){load(paste("R_datasets/Simulated_data_SIM_",1,".RData",sep=""))
    theta.true=theta.sim.out
  }
  
  col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  # Orange, blue, green, pink
  
  for(loadseed in loadpick){
  
    if(simDat==F){loadseedA=paste(loadseed,"_",flu.type,sep="")}else{loadseedA=paste("SIM_",loadseed,sep="")}
    load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseedA,"_lin",linearFn,".RData",sep=""))
    

    # Define lengths and sizes of inputs
    lik.tot=rowSums(likelihoodtab)
    if(is.null(runsPOST)){
      runsPOST = min(8e5, length(lik.tot[lik.tot!=-Inf]) )
    }
    maxlik=max(lik.tot[1:runsPOST])
    #plot(as.data.frame(thetatab)$mu[runs1:runsPOST],type="l",ylab="mu")
    
    # - - - - - - - 
    # Calculate ESS by burn-in
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
  
  # Plot MCMC chains
  if(flu.type=="H3FS" ){ par(mfrow=c(2,3)) }else{ par(mfrow=c(2,4)) }
  
  par(xpd=F)
  par(mar = c(3.5,3.5,1,1))
  par(mgp=c(1.8,0.6,0))
  colA=rgb(0.8,0.8,0.8)
  Ctrue="black";Wtrue=2
  
  maxlik=max(storeLik)
  plot(1:length(storeLik[1,]),storeLik[1,],type="l",col=col.list[[1]],xlab="iteration",ylab="likelihood",ylim=c(maxlik-ifelse(simDat==F,1000,500),maxlik+100) ); 
  for(ii in 2:length(loadpick)){ lines(storeLik[ii,],type="l",col=col.list[[ii]] , ylim=c(maxlik-ifelse(simDat==F,1000,500),maxlik+100)) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(maxlik-ifelse(simDat==F,1000,500),maxlik+100),col="gray",lty=2)

  if(!(flu.type=="H3FS" )){ 
    plot(storeWane[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(omega),ylim=c(0,1.2)); for(ii in 2:length(loadpick)){ lines(storeWane[ii,],type="l",col=col.list[[ii]]) }
    lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,1.2),col="gray",lty=2)
    if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["wane"]],theta.true[["wane"]]),col=Ctrue,lwd=Wtrue) }
  }
  plot(storeMu[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(mu[1]),ylim=c(0,4)); for(ii in 2:length(loadpick)){ lines(storeMu[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,4),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST),c(theta.true[["mu"]],theta.true[["mu"]]),col=Ctrue,lwd=Wtrue) }
  
  if(!(flu.type=="H3FS")){ 
    plot(storeMu2[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(mu[2]),ylim=c(0,4)); for(ii in 2:length(loadpick)){ lines(storeMu2[ii,],type="l",col=col.list[[ii]]) }
    lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,4),col="gray",lty=2)
    if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["muShort"]],theta.true[["muShort"]]),col=Ctrue,lwd=Wtrue) }
  }
  
  plot(storeSigma[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(sigma[1]),ylim=c(0,0.4)); for(ii in 2:length(loadpick)){ lines(storeSigma[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,0.4),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["sigma"]],theta.true[["sigma"]]),col=Ctrue,lwd=Wtrue) }
  
  if(!(flu.type=="H3FS" )){ 
    plot(storeSigma2[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(sigma[2]),ylim=c(0,0.4)); for(ii in 2:length(loadpick)){ lines(storeSigma2[ii,],type="l",col=col.list[[ii]]) }
    lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,0.4),col="gray",lty=2)
    if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["sigma2"]],theta.true[["sigma2"]]),col=Ctrue,lwd=Wtrue) }
  }
  
  plot(storeError[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(epsilon),ylim=c(0,2.5)); for(ii in 2:length(loadpick)){ lines(storeError[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,2.5),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["error"]],theta.true[["error"]]),col=Ctrue,lwd=Wtrue) }
  
  plot(storeTau2[1,],type="l",col=col.list[[1]],xlab="iteration",ylab=expression(tau),ylim=c(0,0.15)); for(ii in 2:length(loadpick)){ lines(storeTau2[ii,],type="l",col=col.list[[ii]]) }
  lines(c(burnCut*ltheta/switch1,burnCut*ltheta/switch1),c(0,0.15),col="gray",lty=2)
  if(simDat==T){ lines(c(-1000,runsPOST/switch1),c(theta.true[["tau2"]],theta.true[["tau2"]]),col=Ctrue,lwd=Wtrue) }
  
  dev.copy(png,paste("plot_simulations/MCMC_chains",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,"_lin",linearFn,".png",sep=""),units="cm",width=20,height=10,res=300)
  dev.off()
  
  
}
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot expected titres using sampled posterior estimates against true titres 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# plot.posterior.titres(loadseed="SIM",simDat=T,define.year=c(2007:2012))

plot.posterior.titres<-function(loadseed=1,year_test=c(2007:2012),flu.type,simDat=F,btstrap=5,plotRes=F,linearFn=F){
  
  # simDat=F; flu.type="H1" ; year_test= c(2009:2011); btstrap=5; 
  # simDat=T;loadseed=1;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; plot.corr=F; linearFn=T;btstrap=5
  
  if(simDat==F){
    if(flu.type=="H3HN"){load("R_datasets/HaNam_data.RData")}
    if(flu.type=="H3FS"){load("R_datasets/FluScapeH3_data.RData")}
    if(flu.type=="B"){load("R_datasets/Fluscape_data.RData")}
    if(flu.type=="H1"){load("R_datasets/HK_data.RData")}
    #loadseed=1 # DEBUG
    loadseed=paste(loadseed,"_",flu.type,sep="")
  }else{
    load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))
    test.list=test.listSim
    hist.true=historytabSim
    theta.true=theta.sim.out
  }
  
  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,"_lin",linearFn,".RData",sep="")) # Note that this includes test.listPost
  
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
    pickA=sample(c(runs1:runsPOST),1)
    
    pickAhist=ceiling(pickA/hist.switch)+1 # check which history this specifies
    
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
    
  for(pickyr in 1:n.test){
    
    # Note here that inf_years and strain_years are loads from main_model.R
    # Output expected titre - could include Poisson measurement uncertainty here?
    #print("Need to add antigenic.map.in to SIMULATION") ?
    simulate_data(test.yr[pickyr],historytabPost=hist.sample,
                  inf_years,
                  strain_years,
                  n_part,thetastar=theta.max,p.inf=0.1,
                  #pmask=c("sigma2"), # For old fitted data, need to specify that sigma2 wasn't fitted
                  linD=linearFn,
                  am.spline=am.spl)
    
    load("R_datasets/Simulated_dataPost_1.RData")
    
    # Mask infections after test year
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
  if(plotRes==T){
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
  
  par(mfrow=c(1,1))
  par(mar = c(5,5,1,1))
  breaks0=seq(round(min(c(compTab,compTabNULL)))-0.5,round(max(c(compTab,compTabNULL)))+1.5,1)
  hist(compTab , breaks = breaks0 , freq=F, col =rgb(0.8,0.8,0.8),main=NULL,xlab="expected titre - observed titre"); #title(LETTERS[1],adj=0)
  #hist(compTabNULL , breaks = breaks0 , freq=F, col =rgb(0.8,0.8,0.8),main=NULL,xlab="expected titre - Pois(expected titre)"); #title(LETTERS[2],adj=0)
  
  stats.store = rbind( c(sum( abs(compTab)<=1 )/length(compTab), sum( abs(compTab)<=2 )/length(compTab), sum( abs(compTab)<=3 )/length(compTab) ), 
                       c(sum( abs(compTabNULL)<=1 )/length(compTabNULL), sum( abs(compTabNULL)<=2 )/length(compTabNULL), sum( abs(compTabNULL)<=3 )/length(compTabNULL)  )
                )
  colnames(stats.store)=c("within 1","within 2","within 3")
  rownames(stats.store)=c("model","null")
  write.csv(stats.store,paste("plot_simulations/titre_compare/Store_residuals_",loadseed,".csv",sep=""))
  
  dev.copy(pdf,paste("plot_simulations/titre_compare/Residuals_",loadseed,"_lin",linearFn,".pdf",sep=""),width=5,height=4)
  dev.off()
  
  
  # - - - - - - - - - - - - - - - - - - - 
  # Calculate likelihood of data

  
  }
  
  # - - - - - - - - - - - - - - - - - - - 
  # Plot figures from MCMC posteriors

    colN=if(n.test>1){n.test}else{5}
    loopN=if(n.test>1){5}else{25}
    par(mfrow=c(5,colN)); 
    par(mar = c(2,2,1.5,0.2));
    par(mgp=c(2,0.5,0))
    # Mask infections after test year
    
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
        #dev.copy(pdf,paste("plot_simulations/titre_compare/sim",loadseed,"_",ii0,"P.pdf",sep=""),width=10,height=8)
        dev.copy(png,paste("plot_simulations/titre_compare/sim",loadseed,"_",ii0,"P.png",sep=""),units="cm",width=25,height=15,res=150)
        dev.off()
      }
    
    }  # end loop over participants
      
  } # end loop over test years

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# FIGURE PLOT Plot expected titres using sampled posterior estimates against true titres 

plot.posterior.titres.select<-function(loadseed=1,year_test=c(2007:2012),flu.type,simDat=F,btstrap=5,part_pick=c(31,57,25),year_pick=c(2008:2010),linearFn=T){
  
  # year_pick=c(2008:2010);part_pick=c(15,31,57)
  
  # year_pick=c(2008:2010);part_pick=c(31,57,25) ; simDat = F;  year_test=c(2007:2012); btstrap=50; loadseed = 1; flu.type="H3HN"; linearFn=T
  
  # year_pick=c(2009);part_pick=c(31,57,25,1,2,3) ; simDat = F;  year_test=c(2009); btstrap=50; loadseed = 1; flu.type="H3FS"; linearFn=T
  
  
  if(simDat==F){
    if(flu.type=="H3HN"){load("R_datasets/HaNam_data.RData") }
    if(flu.type=="H3FS"){load("R_datasets/FluscapeH3_data.RData") }
    if(flu.type=="H1"){load("R_datasets/HK_data.RData")}
    loadseed=1 # DEBUG
    loadseed=paste(loadseed,"_",flu.type,sep="")
  }else{
    load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))
    test.list=test.listSim
    hist.true=historytabSim
  }
  
  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,"_lin",linearFn,".RData",sep="")) # Note that this includes test.listPost
  
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.2*runsPOST)
  hist.switch= runsPOST/ (length(historytabCollect[,1])/n_part - 1)
  
  # Set up matrices to store -- need btstrap >1
  n.strains=length(strain_years) # this loads from main_model.R
  n.testMatch=match(year_pick,test.yr)
  n.test=length(year_pick)
  n.inf=length(inf_years)
  n.ppart=length(part_pick)
  
  store.mcmc.test.data=array(NA, dim=c(btstrap,n.ppart,n.strains,n.test,2)) # Store expected titres for each test year
  store.mcmc.hist.data=array(NA, dim=c(btstrap,n.ppart,n.inf,n.test)) # Store history for each test year
  
  # - - - - - - - - - - - -
  # Sample from MCMC runs to get stored matrices of expected titre and estimated infection years
  
  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1)
    pickAhist=ceiling(pickA/hist.switch)+1 # check which history this specifies -- **BASED ON 20 step resampling
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
    
    for(pickyr in 1:n.test){
      
      # Note here that inf_years and strain_years are loads from main_model.R
      # Output expected titre - could include Poisson measurement uncertainty here?
      simulate_data(year_pick[pickyr],historytabPost=hist.sample,
                    inf_years,
                    strain_years,
                    n_part,thetastar=theta.max,p.inf=0.1,
                    #pmask=c("sigma2"), # For old fitted data, need to specify that sigma2 wasn't fitted
                    linD=linearFn)
      
      load("R_datasets/Simulated_dataPost_1.RData")
      
      # Mask infections after test year
      for(ii0 in 1:n.ppart){
        
        sim.titre=test.listSim[[ part_pick[ii0] ]][[1]] # sort sample years - drawn from simulated data above
        hist.sampleB=hist.sample;  hist.sampleB[,as.numeric(colnames(hist.sample))> year_pick[pickyr] ]=0 # don't show infections after test year
        store.mcmc.test.data[sampk,ii0,,pickyr,1]= min(inf_years)-1+sort(sim.titre["sample.index",]) # Sampled strain years
        s.titre=sim.titre["titredat",order(sim.titre["sample.index",])]
        store.mcmc.test.data[sampk,ii0,,pickyr,2]=s.titre # Sampled expected titre
        #store.mcmc.test.data[sampk,ii0,,pickyr,2]=rpois(length(s.titre),lambda=s.titre) # Sampled expected titre - include Poisson noise
        store.mcmc.hist.data[sampk,ii0,,pickyr]=hist.sampleB[ part_pick[ii0] ,] # Sampled history
        
      }  # end loop over participants
      
    } # end loop over test years
  } # end bootstrap loop
  
  # - - - - - - - - - - - - 
  # Plot figures from MCMC posteriors
  
  if(flu.type=="H3FS"){par(mfrow=c(4,3))}else{par(mfrow=c(4,3))} #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  par(mgp=c(1.8,0.6,0))
  ymax=8.2
  titleN = ifelse(flu.type=="H3FS",10,1)
  
  yrange1=c(-0.2,9)
  
  for(ii0 in 1:n.ppart){

    for(pickyr in 1:n.test){
    
      simtitreX=store.mcmc.test.data[,ii0,,pickyr,1]
      simtitreY=store.mcmc.test.data[,ii0,,pickyr,2]
      hist.sample=store.mcmc.hist.data[,ii0,,pickyr] # for participant ii0 in year pickyr

      #if(ii0 ==3 & pickyr == 1 ){par(mar = c(5,5,1,1))}
      #par(mar = c(2,2,1,1)) #B L T R

      par(mar = c(2,3,1,1))
      plot(inf_years,8*hist.sample[1,],type="l",yaxs="i",ylim=yrange1,col='white',xlab=ifelse(ii0==3,"",""),ylab="log titre",main=year_pick[pickyr],xlim=c(1966,2011))
      axis(side = 1, at = inf_years, labels = FALSE, tck = -0.01)
      
      # Sample from infection history
      for(jj in 1:n.inf){
        lines(min(inf_years)-1+c(jj,jj),c(yrange1[1]+0.05,yrange1[2]-0.05),col=rgb(0,0,0,sum(hist.sample[,jj])/btstrap),lwd=2) # Plot estimated infections
      }

      # Calculate credible interval for expected titres
      medP=apply(simtitreY,2,function(x){median(x)}); ciP195=apply(simtitreY,2,function(x){quantile(x,0.025)}); ciP295=apply(simtitreY,2,function(x){quantile(x,0.975)});
      ciP1=apply(simtitreY,2,function(x){quantile(x,0.25)}); ciP2=apply(simtitreY,2,function(x){quantile(x,0.75)})
      #alpha50 = qchisq(0.5, 1); alpha95 = qchisq(0.95, 1); medP=apply(simtitreY,2,function(x){mean(x)})
      #ciP1=(medP + alpha50/2) - sqrt(alpha50)*sqrt(medP + alpha50/4); ciP2=(medP + alpha50/2) + sqrt(alpha50)*sqrt(medP + alpha50/4) # Poisson CI
      #ciP195=(medP + alpha95/2) - sqrt(alpha95)*sqrt(medP + alpha95/4); ciP295=(medP + alpha95/2) + sqrt(alpha95)*sqrt(medP + alpha95/4) # Poisson CI
      
      polygon(c(simtitreX[1,],rev(simtitreX[1,])),c(ciP195,rev(ciP295)),lty=0,col=rgb(0,0.3,1,0.1))
      polygon(c(simtitreX[1,],rev(simtitreX[1,])),c(ciP1,rev(ciP2)),lty=0,col=rgb(0,0.3,1,0.2))
      lines(simtitreX[1,],medP,pch=1,col='blue')
      points(simtitreX[1,],medP,pch=19,cex=0.5,col='blue')
      
      # Plot true titres - check select correct values
      points(min(inf_years)-1+test.list[[ part_pick[ii0] ]][[ n.testMatch[pickyr] ]][4,],test.list[[ part_pick[ii0] ]][[ n.testMatch[pickyr] ]][2,],pch=1,col='red')
      
      text(x=1965,y=8,labels=LETTERS[titleN],adj=0,font=2,cex=1.5)
      titleN=titleN+1
      
    }  # end loop over participants
  } # end loop over test years
  
  dev.copy(pdf,paste("plot_simulations/titre_compare/FIGURE_titre_plot",loadseed,".pdf",sep=""),width=6.5,height=7,useDingbats=F)
  dev.off()
  
  
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot simulated data

plot.sim.data<-function(){

  loadseed="SIM"
  load(paste("R_datasets/Simulated_data_",loadseed,"_1.RData",sep=""))
  
  #load("R_datasets/Simulated_dataPost_1.RData")
  
  #Compare model fits using posterior infection history (historytabPost) and parameters
  picktest=c(2007:2012)
  n.test=length(picktest)
  
  simulate_data(test_years=picktest,historytabPost=historytabSim,
                inf_years,
                strain_years,
                n_part=npartM,thetastar=thetaSim,p.inf=0.1) #theta.max
  
  load("R_datasets/Simulated_dataPost_1.RData")
  par(mfrow=c(2,5))
  par(mar = c(5,5,1,1))
  for(pickyr in 1:n.test){
    for(ii0 in 1:n_part){
      plot(8*historytabSim[ii0,],type="l",ylim=c(0,9),col='white')
      #for(jj in 1:length(lenhis)){
      #  lines(c(jj,jj),c(0,9*historytabSim[ii0,jj]),col='red')
      #}
      histSim1=historytabSim[ii0,]
      histSim1[inf_years>picktest[pickyr]]=0 # don't show infections after test year
      lenhis=rep(0,length(histSim1))
      
      for(jj in 1:length(lenhis)){
        lines(c(jj,jj),c(0,9*histSim1[jj]),col='blue')
      }
      lines(test.listSim[[ii0]][[pickyr]][4,],test.listSim[[ii0]][[pickyr]][2,],col=rgb(0.8,0.8,0.8),lwd=2) # Plot estimated infections
      points(test.listSim[[ii0]][[pickyr]][4,],test.listSim[[ii0]][[pickyr]][2,],pch=19)
      if(ii0 %% 10==0){
        dev.copy(pdf,paste("plot_simulations/simPlot",ii0,"P_",pickyr,".pdf",sep=""),width=12,height=6)
        dev.off()
      }
    }
  }

}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Rewind history and run with flat incidence - DEPRECATED

run.titre.time<-function(loadseed=1,year_test=c(2007:2012),flu.type="H3HN",simDat=F,btstrap=5,n_partSim=2,simTest.year=c(1968:2010)){
  
  if(simDat==F){
    if(flu.type=="H3HN"){load("R_datasets/HaNam_data.RData")}
    if(flu.type=="B"){load("R_datasets/Fluscape_data.RData")}
    if(flu.type=="H1"){load("R_datasets/HK_data.RData")}
    loadseed=paste(loadseed,"_",flu.type,sep="")
  }else{
    load(paste("R_datasets/Simulated_data_",loadseed,"_1.RData",sep=""))
  }
  
  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep="")) # Note that this includes test.listPost
  
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.2*runsPOST)
  
  # Set up matrices to store -- need btstrap >1
  strain_years=inf_years # look at strains from every year
  n.strains=length(strain_years) # this loads from main_model.R
  n.inf=length(inf_years)
  hist.sample0=rep(c(1,rep(0,29)),100)[1:n.inf] # CURRENTLY JUST FOR ONE PARTICIPANT
  simTest.year=sort(c(inf_years[hist.sample0==1],inf_years[hist.sample0==1]+1)) # infection year and one year after
  n.test=length(simTest.year)

  store.mcmc.test.data=array(NA, dim=c(btstrap,n_partSim,n.strains,n.test,2)) # Store expected titres for each test year
  store.mcmc.hist.data=array(NA, dim=c(btstrap,n_partSim,n.inf,n.test)) # Store history for each test year
  
  # - - - - - - - - - - - -
  # Sample from MCMC runs to get stored matrices of expected titre and estimated infection years
  
  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1)
    pickAhist=ceiling(pickA/20)+1 # check which history this specifies
    hist.sample=rbind(hist.sample0,hist.sample0)
    theta.max=as.data.frame(thetatab)[pickA,]
    
    for(pickyr in 1:n.test){ # ITERATE OVER TIME HERE
      
      # Note here that inf_years and strain_years are loads from main_model
      #hist.sample0[ inf_years<simTest.year[pickyr] ] # only take years up to test year -- already included in simulation function!
      simulate_data(simTest.year[pickyr],historytabPost=hist.sample,
                    inf_years,
                    strain_years,
                    n_partSim,thetastar=theta.max,p.inf=0.1,
                    #pmask=c("sigma2"), # For old fitted data, need to specify that sigma2 wasn't fitted
                    linD=F)
      
      load("R_datasets/Simulated_dataPost_1.RData")
      
      # Mask infections after test year
      for(ii0 in 1){
        
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
  
  # - - - - - - - - - - - - 
  # Plot development of titres
  
  par(mfrow=c(1,1)); par(mar = c(5,5,1,1))
  
  plot(inf_years,8*hist.sample[1,],type="l",ylim=c(0,9),col='white',xlab="year",ylab="titre")
  
  for(pickyr in 1:n.test){

    # Mask infections after test year
    
    for(ii0 in 1){
      simtitreX=store.mcmc.test.data[,ii0,,pickyr,1]
      simtitreY=store.mcmc.test.data[,ii0,,pickyr,2]
      hist.sample=store.mcmc.hist.data[,ii0,,pickyr] # for participant ii0 in year pickyr
      
      # Sample from infection history
      for(ksamp in 1:btstrap){
        for(jj in 1:n.inf){
          lines(min(inf_years)-1+c(jj,jj),c(-1,12*hist.sample[ksamp,jj]-1),col=rgb(0.8,0.8,0.8,0.01),lwd=2) # Plot estimated infections
        }
      }
      
      # Calculate credible interval for expected titres
      medP=apply(simtitreY,2,function(x){median(x)})
      ciP1=apply(simtitreY,2,function(x){quantile(x,0.025)})
      ciP2=apply(simtitreY,2,function(x){quantile(x,0.975)})
      polygon(c(simtitreX[1,],rev(simtitreX[1,])),c(ciP1,rev(ciP2)),lty=0,col=rgb(0,0.3,1,0.2))
      lines(simtitreX[1,],medP,pch=1,col='blue')
      points(simtitreX[1,],medP,pch=19,cex=0.5,col='blue')
      
      if(ii0 %% 10==0){
        dev.copy(pdf,paste("plot_simulations/sim",ii0,"P_",pickyr,".pdf",sep=""),width=12,height=6)
        dev.off()
      }
      
    }  # end loop over participants
    
  } # end loop over test years
  
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot different aspects of immune response for paper figures

plot.antibody.changes<-function(loadseed=1,year_test=c(2007:2012),flu.type="H3HN",simDat=F,btstrap=5,n_partSim=2,simTest.year=c(1968:2010)){
  
  # btstrap=50 ; n_partSim=2 ; simTest.year=c(1968:2010)
  load("R_datasets/HaNam_data.RData")
  loadseed="1_H3" #paste(loadseed,"_",flu.type,sep="")

  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep="")) # Note that this includes test.listPost

  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.2*runsPOST)
  thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
  
  sampleX = 10
  year_x = seq(0,3,0.01) ;   year_x2 = seq(0,3,0.1) ; cross_x = seq(0,30,0.1) ;  cross_x2 = seq(0,30,1)
  Dmed = 3; Dlong = 10; Tmed = 0.5; Tlon = 3 # Define values for plots
  
  # 1D plots
  # Store bootstrap runs
  storetitreS = NULL; storetitreM = NULL; storetitreL = NULL;
  storetitreCrossS = NULL ; storetitreCrossM = NULL; storetitreCrossL = NULL
  
  for(sampk in 1:btstrap){
    # Sample from MCMC posterior to get trajectory
    pickA=sample(c(runs1:runsPOST),1)
    theta.max=as.data.frame(thetatab)[pickA,]
    
    storetitreS = rbind(storetitreS, sapply( theta.max$mu + theta.max$muShort*exp(- year_x * theta.max$wane), function(x){min(8, x)}) )
    #storetitreM = rbind(storetitreM, sapply( theta.max$mu*exp(-theta.max$sigma*Dmed) + theta.max$muShort*exp(- year_x * theta.max$wane)*exp(-theta.max$sigma2*Dmed), function(x){min(8, x)}) )
    storetitreL = rbind(storetitreL, sapply( theta.max$mu*exp(-theta.max$sigma*Dlong) + theta.max$muShort*exp(- year_x * theta.max$wane)*exp(-theta.max$sigma2*Dlong), function(x){min(8, x)}) )
    
    storetitreCrossS = rbind(storetitreCrossS, sapply( theta.max$mu*exp(-theta.max$sigma*abs(cross_x)) + theta.max$muShort*exp(-theta.max$sigma2*abs(cross_x)), function(x){min(8, x)}) )
    #storetitreCrossM = rbind(storetitreCrossM, sapply( theta.max$mu*exp(-theta.max$sigma*abs(cross_x)) + theta.max$muShort*exp(-theta.max$sigma2*abs(cross_x))*exp(- Tmed * theta.max$wane), function(x){min(8, x)}) )
    storetitreCrossL = rbind(storetitreCrossL, sapply( theta.max$mu*exp(-theta.max$sigma*abs(cross_x)) + theta.max$muShort*exp(-theta.max$sigma2*abs(cross_x))*exp(- Tlon * theta.max$wane), function(x){min(8, x)}) )
    
  }
  
  titreL = apply(storetitreL,2,function(x){c.nume(x)})
  #titreM = apply(storetitreM,2,function(x){c.nume(x)})
  titreS = apply(storetitreS,2,function(x){c.nume(x)})
  titreCrossL = apply(storetitreCrossL,2,function(x){c.nume(x)})
  #titreCrossM = apply(storetitreCrossM,2,function(x){c.nume(x)})
  titreCrossS = apply(storetitreCrossS,2,function(x){c.nume(x)})
  
  col1P=rgb(0,0,0.8) # rgb(0,0.5,0)
  col1=rgb(0,0,0.8,0.005) # rgb(0,0.5,0,0.03)
  col1a=rgb(0,0,0.8,0.1) # rgb(0,0.5,0,0.03)
  col2P=rgb(0,0,1)
  col2=rgb(0,0,1,0.2)
  col3P=rgb(1,0,0.5)
  col3=rgb(1,0,0.2,0.2)
  
  pois.gen <- function(x,rep){rpois(rep,lambda = x)} # Function to generate random samples
  
  # - - - - - - - - - - - - - - - - - - - - - - 
  # PLOT PARAMETERS FIGURE 3
  # 1. Plot titre waning
  
  par(mfrow=c(1,2)); par(mar = c(3,3,1,1))
  par(mgp=c(1.8,0.6,0))
  
  plot(year_x,titreS[1,],ylim=c(-0.1,8.5),type="l",ylab="log titre",xlab="years since infection",col=NULL,xaxs="i",yaxs="i")
  polygon(c(year_x,rev(year_x)),c(titreS[2,],rev(titreS[3,])),col=col1a,lty=0)
  #polygon(c(year_x,rev(year_x)),c(titreM[2,],rev(titreM[3,])),col=col2,lty=0)
  #polygon(c(year_x,rev(year_x)),c(titreL[2,],rev(titreL[3,])),col=col3,lty=0)
  
  for(ii in 1:length(year_x2)){
    titreSsim = pois.gen(storetitreS[,(ii-1)*sampleX+1],btstrap) # Simulate from Poisson - sample fewer than previous
    tstore=data.frame(table(titreSsim),stringsAsFactors = F)
    t.alpha=tstore$Freq/btstrap
    points(rep(year_x2[ii],length(tstore$titreSsim)),as.numeric(levels(tstore$titreSsim))[tstore$titreSsim],col=rgb(0,0,0.8,t.alpha),pch=19,cex=0.5)
  }

  lines(year_x,titreS[1,],col=col1P ,lwd = 2)  
  #lines(year_x,titreM[1,],col=col2P)
  #lines(year_x,titreL[1,],col=col3P)
  
  title(main=LETTERS[1],adj=0)
  
  # 2. Plot cross reaction
  par(mar = c(3,2,1,1))
  
  plot(cross_x,titreCrossS[1,],ylim=c(-0.1,8.5),type="l",ylab="",xlab="antigenic distance",col=NULL,xaxs="i",yaxs="i")
  polygon(c(cross_x,rev(cross_x)),c(titreCrossS[2,],rev(titreCrossS[3,])),col=col1a,lty=0)
  #polygon(c(cross_x,rev(cross_x)),c(titreCrossM[2,],rev(titreCrossM[3,])),col=col2,lty=0)
  polygon(c(cross_x,rev(cross_x)),c(titreCrossL[2,],rev(titreCrossL[3,])),col=col3,lty=0)

  for(ii in 1:length(cross_x2)){
    titreSsim = pois.gen(storetitreCrossS[,(ii-1)*sampleX+1],btstrap) # Simulate from Poisson - sample fewer than previous
    tstore=data.frame(table(titreSsim),stringsAsFactors = F) ; t.alpha=tstore$Freq/btstrap
    points(rep(cross_x2[ii],length(tstore$titreSsim)),as.numeric(levels(tstore$titreSsim))[tstore$titreSsim],col=rgb(0,0,0.8,t.alpha),pch=19,cex=0.5)
    
    titreSsim = pois.gen(storetitreCrossL[,(ii-1)*sampleX+1],btstrap)
    tstore=data.frame(table(titreSsim),stringsAsFactors = F) ; t.alpha=tstore$Freq/btstrap
    points(rep(cross_x2[ii],length(tstore$titreSsim)),as.numeric(levels(tstore$titreSsim))[tstore$titreSsim],col=rgb(1,0,0.5,t.alpha),pch=19,cex=0.5)
  }
  
  lines(cross_x,titreCrossS[1,],col=col1P,lwd = 2)
  #lines(cross_x,titreCrossM[1,],col=col2P)
  lines(cross_x,titreCrossL[1,],col=col3P)
  
  title(main=LETTERS[2],adj=0)
  
  dev.copy(pdf,paste("plot_simulations/parameters",loadseed,".pdf",sep=""),width=7,height=4)
  dev.off()
  
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Rewind history and reconstruct antibody landscape for Figure 2

run.historical.landscapes<-function(loadseed=1,year_test=c(2007:2012),linearFn=F,flu.type="H3HN",simDat=F,btstrap=5,n_partSim=2,simTest.year=c(1968:2010),d.step=0.5,ymax=5){
    
    # btstrap=50 ; n_partSim=2 ; simTest.year=c(1968:2010) ; d.step = 0.5 ; flu.type="H3HN"; year_test=c(2007:2012); loadseed = 1; linearFn=T
    load("R_datasets/HaNam_data.RData")
    loadseed=paste(loadseed,"_",flu.type,sep="")
    
    load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,"_lin",linearFn,".RData",sep="")) # Note that this includes test.listPost
    
    lik.tot=rowSums(likelihoodtab)
    runsPOST=length(lik.tot[lik.tot!=-Inf])
    runs1=ceiling(0.2*runsPOST)
    thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
    
    # Read in antigenic map data and specify test strains
    
    ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
    strain_names=ag.coord$viruses
    
    x.range <- seq(floor(min(ag.coord$AG_x))-2,ceiling(max(ag.coord$AG_x))+1,d.step)
    y.range <- seq(floor(min(ag.coord$AG_y))-2,ceiling(max(ag.coord$AG_y)),d.step)
    points.j <- expand.grid(x.range,y.range) # Define list of points to evaluate
    
    xx=scalemap(inf_years,inf_years)
    yy=predict(am.spl,xx)$y
    
    # Set up matrices to store -- need btstrap >1
    strain_years=inf_years # look at strains from every year
    n.strains=length(strain_years) # this loads from main_model.R
    n.inf=length(inf_years)
    
    pick_years = c(1,21) #c(1,31)
    hist.sample0=rep(c(rep(0,pick_years[1]-1),1,rep(0,pick_years[2]-1)),100)[1:n.inf] # CURRENTLY JUST FOR ONE PARTICIPANT
    simTest.year=sort(c(inf_years[hist.sample0==1],inf_years[hist.sample0==1]+1)) # infection year and one year after
    n.test=length(simTest.year)

    
    # Set up sample sapce
    
    sampleX = 10

    # 1D plots
    # Store bootstrap runs
    storetitreS1 <- rep(NA,nrow=length(points.j[,1]))
    storetitreL1 <- rep(NA,length(points.j[,1]))
    storetitreS2 <- rep(NA,length(points.j[,1]))
    storetitreL2 <- rep(NA,length(points.j[,1]))
    storetitreS3 <- rep(NA,length(points.j[,1]))
    #storetitreL3 <- rep(NA,length(points.j[,1]))
    
    right.censor <- function(x){min(8, x)}
    pickA <- (lik.tot==max(lik.tot)) # pick MAP estimate #sample(c(runs1:runsPOST),1)
    theta.max <- head(as.data.frame(thetatab)[pickA,],1)
    
    # Iterate across whole sample sapce
    for(sampk in 1:length(points.j[,1])){
      year_x=pick_years[2]-pick_years[1]
      
      # Need to add AGS
      
      cross_x1 = sqrt( sum((as.numeric(points.j[sampk,]) - c(yy[pick_years[1]],xx[pick_years[1]]) )^2 ) ) # calculate 2D distance - note xx,yy switch
      cross_x2 = sqrt( sum((as.numeric(points.j[sampk,]) - c(yy[pick_years[2]],xx[pick_years[2]]) )^2 ) ) # calculate 2D distance
      
      
      # Need to edit the linear function here 
      if(linearFn==F){
      
        storetitreS1[sampk] = (theta.max$mu*exp(-theta.max$sigma*cross_x1) + theta.max$muShort*exp(- 0 * theta.max$wane)*exp(-theta.max$sigma2*cross_x1) ) %>% right.censor()
        storetitreL1[sampk] = ( theta.max$mu*exp(-theta.max$sigma*cross_x1) + theta.max$muShort*exp(- (1) * theta.max$wane)*exp(-theta.max$sigma2*cross_x1) )%>% right.censor()
        
        storetitreS2[sampk] = ( (1+theta.max$tau1)*(theta.max$mu*exp(-theta.max$sigma*cross_x1) + 
                        theta.max$muShort*exp(- (year_x) * theta.max$wane)*exp(-theta.max$sigma2*cross_x1) ) + 
                          exp(- theta.max$tau2)*(theta.max$mu*exp(-theta.max$sigma*cross_x2) + 
                        theta.max$muShort*exp(- 0 * theta.max$wane)*exp(-theta.max$sigma2*cross_x2) )
                        ) %>% right.censor()
        
        storetitreL2[sampk] = ( (1+theta.max$tau1)*(theta.max$mu*exp(-theta.max$sigma*cross_x1) + 
                        theta.max$muShort*exp(- (year_x+1 ) * theta.max$wane)*exp(-theta.max$sigma2*cross_x1) ) + # Note the year_x +1 here
                        exp(- theta.max$tau2)*(theta.max$mu*exp(-theta.max$sigma*cross_x2) + 
                        theta.max$muShort*exp(- (1) * theta.max$wane)*exp(-theta.max$sigma2*cross_x2) )
                        ) %>% right.censor()
      
      }else{
        
        storetitreS1[sampk] = (theta.max$mu*max(0,1-theta.max$sigma*cross_x1) + theta.max$muShort*max(0,1- 0 * theta.max$wane)*max(0,1-theta.max$sigma2*cross_x1) ) %>% right.censor()
        storetitreL1[sampk] = ( theta.max$mu*max(0,1-theta.max$sigma*cross_x1) + theta.max$muShort*max(0,1 - (1) * theta.max$wane)*max(0,1-theta.max$sigma2*cross_x1) )%>% right.censor()
        
        storetitreS2[sampk] = ( (1+theta.max$tau1)*(theta.max$mu*max(0,1-theta.max$sigma*cross_x1) + 
                                                      theta.max$muShort*max(0,1- (year_x) * theta.max$wane)*(1-theta.max$sigma2*cross_x1) ) + 
                                  max(0,1- theta.max$tau2)*(theta.max$mu*max(0,1-theta.max$sigma*cross_x2) + 
                                                           theta.max$muShort*max(0,1- 0 * theta.max$wane)*(1-theta.max$sigma2*cross_x2) )
        ) %>% right.censor()
        
        storetitreL2[sampk] = ( (1+theta.max$tau1)*(theta.max$mu*max(0,1-theta.max$sigma*cross_x1) + 
                                                      theta.max$muShort*max(0,1- (year_x+1 ) * theta.max$wane)*max(0,1-theta.max$sigma2*cross_x1) ) + # Note the year_x +1 here
                                  max(0,1- theta.max$tau2)*(theta.max$mu*max(0,1-theta.max$sigma*cross_x2) + 
                                                           theta.max$muShort*max(0,1- (1) * theta.max$wane)*max(0,1-theta.max$sigma2*cross_x2) )
        ) %>% right.censor()
        
      }

    }
    

    
    
    # PLOT FIGURES

    col1P=rgb(0,0,0.8) # rgb(0,0.5,0)
    
    pred_matrixS1 <- matrix(storetitreS1,byrow=F,nrow=length(x.range))
    pred_matrixL1 <- matrix(storetitreL1,byrow=F,nrow=length(x.range))
    pred_matrixS2 <- matrix(storetitreS2,byrow=F,nrow=length(x.range))
    pred_matrixL2 <- matrix(storetitreL2,byrow=F,nrow=length(x.range))
    
    # Calculate titres along summary path
    predOutput <- function(matrixA){apply(cbind(xx,yy),1,function(zz){ ydist=abs(y.range-zz[1]); xdist=abs(x.range-zz[2]); matrixA[xdist==min(xdist),ydist==min(ydist)] })}

    # Define list of isolate years
    sY <- strain_years_convert()
    
    # Plot panels and antigenic summary paths
    layout(matrix(c(rep(1,3),rep(2,3),rep(1,3),rep(2,3),rep(1,3),rep(2,3),
                    rep(3,3),rep(4,3),
                    rep(5,3),rep(6,3),rep(5,3),rep(6,3),rep(5,3),rep(6,3),
                    rep(7,3),rep(8,3)
                    ), 8,6, byrow=T) )
    
    par(mgp=c(1.8,0.6,0))
    tranP <- 1; sizP <- 0.9 ; grY <- 0.6
    
    par(mar = c(4,4,2,2))
    image2D(z = t(pred_matrixS1), x = y.range, y = x.range, xlab="antigenic dimension 1", ylab="antigenic dimension 2", zlim = c(0, ymax),
            main=inf_years[pick_years[1] ],col=rev(ramp.col (col = c("blue",rgb(0.4,0.6,1),"white"), n = 100, alpha = 1))) #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")
    points(ag.coord[sY <= (pick_years[1] + 1968 -1),c("AG_y","AG_x")],pch=19,col=rgb(0,0,0,tranP))
    points(ag.coord[sY > (pick_years[1] + 1968 -1),c("AG_y","AG_x")],pch=19,col=rgb(grY,grY,grY,tranP))
    lines(xx,yy,col="black",lty=2)
    
    points(x=xx[pick_years[1]],y=yy[pick_years[1]],cex=1.5, col="red", lwd=2)
    title(main=LETTERS[1],adj=0)

    par(mar = c(4,4,2,2))
    image2D(z = t(pred_matrixL1), x = y.range, y = x.range, xlab="antigenic dimension 1", ylab="antigenic dimension 2", zlim = c(0, ymax),
            main=inf_years[pick_years[1] ]+1,col=rev(ramp.col (col = c("blue",rgb(0.4,0.6,1),"white"), n = 100, alpha = 1))) #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")
    points(ag.coord[sY <= (pick_years[1] + 1968 ),c("AG_y","AG_x")],pch=19,col=rgb(0,0,0,tranP))
    points(ag.coord[sY > (pick_years[1] + 1968 ),c("AG_y","AG_x")],pch=19,col=rgb(grY,grY,grY,tranP))
    lines(xx,yy,col="black",lty=2)
    
    points(x=xx[pick_years[1]],y=yy[pick_years[1]],cex=1.5, col="red", lwd=2)
    title(main=LETTERS[3],adj=0)
    
    infcol="red"
    
    # Project along summary path
    par(mar = c(2,4,1,2))
    plot(inf_years,predOutput(pred_matrixS1),type="l",col="blue",ylim=c(0,ymax),ylab="log titre",xlab="",yaxs="i", lwd=2) # Need to include uncertainty?
    lines(c(inf_years[pick_years[1] ],inf_years[pick_years[1] ]),c(0,ymax),lwd=2,col=infcol)
    title(main=LETTERS[2],adj=0)
    
    plot(inf_years,predOutput(pred_matrixL1),type="l",col="blue",ylim=c(0,ymax),ylab="log titre",xlab="",yaxs="i", lwd=2) # Need to include uncertainty?
    lines(c(inf_years[pick_years[1] ],inf_years[pick_years[1] ]),c(0,ymax),lwd=2,col=infcol)
    title(main=LETTERS[4],adj=0)
    
    par(mar = c(4,4,2,2))
    image2D(z = t(pred_matrixS2), x = y.range, y = x.range, xlab="antigenic dimension 1", ylab="antigenic dimension 2", zlim = c(0, ymax),
            main=inf_years[pick_years[2] ],col=rev(ramp.col (col = c("blue",rgb(0.4,0.6,1),"white"), n = 100, alpha = 1))) #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")
    points(ag.coord[sY <= (pick_years[2] + 1968 -1),c("AG_y","AG_x")],pch=19,col=rgb(0,0,0,tranP))
    points(ag.coord[sY > (pick_years[2] + 1968 -1),c("AG_y","AG_x")],pch=19,col=rgb(grY,grY,grY,tranP))
    lines(xx,yy,col="black",lty=2)
    
    points(x=c(xx[pick_years[1]],xx[pick_years[2]]),y=c(yy[pick_years[1]],yy[pick_years[2]]),cex=1.5, col="red", lwd=2)
    title(main=LETTERS[5],adj=0)
    
    par(mar = c(4,4,2,2))
    image2D(z = t(pred_matrixL2), x = y.range, y = x.range, xlab="antigenic dimension 1", ylab="antigenic dimension 2", zlim = c(0, ymax),
            main=inf_years[pick_years[2] ]+1,col=rev(ramp.col (col = c("blue",rgb(0.4,0.6,1),"white"), n = 100, alpha = 1))) #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")
    points(ag.coord[sY <= (pick_years[2] + 1968 ),c("AG_y","AG_x")],pch=19,col=rgb(0,0,0,tranP))
    points(ag.coord[sY > (pick_years[2] + 1968 ),c("AG_y","AG_x")],pch=19,col=rgb(grY,grY,grY,tranP))
    lines(xx,yy,col="black",lty=2)
    
    points(x=c(xx[pick_years[1]],xx[pick_years[2]]),y=c(yy[pick_years[1]],yy[pick_years[2]]),cex=1.5, col="red", lwd=2)
    title(main=LETTERS[7],adj=0)
    
    # Project along summary path
    par(mar = c(2,4,1,2))
    plot(inf_years,predOutput(pred_matrixS2),type="l",col="blue",ylim=c(0,ymax),ylab="log titre",xlab="",yaxs="i", lwd=2) # Need to include uncertainty?
    lines(c(inf_years[pick_years[1] ],inf_years[pick_years[1] ]),c(0,ymax),lwd=2,col=infcol)
    lines(c(inf_years[pick_years[2] ],inf_years[pick_years[2] ]),c(0,ymax),lwd=2,col=infcol)
    title(main=LETTERS[6],adj=0)
    plot(inf_years,predOutput(pred_matrixL2),type="l",col="blue",ylim=c(0,ymax),ylab="log titre",xlab="",yaxs="i", lwd=2) # Need to include uncertainty?
    lines(c(inf_years[pick_years[1] ],inf_years[pick_years[1] ]),c(0,ymax),lwd=2,col=infcol)
    lines(c(inf_years[pick_years[2] ],inf_years[pick_years[2] ]),c(0,ymax),lwd=2,col=infcol)
    title(main=LETTERS[8],adj=0)
    
    
    dev.copy(pdf,paste("plot_simulations/simulate_new_response/map_space",loadseed,".pdf",sep=""),width=8,height=7,useDingbats=F)
    dev.off()

    
}

# Plot H3N2 Vietnam data

plot_h3_reports <- function(){
  
  par(mfrow=c(3,1))
  par(mgp=c(1.8,0.6,0))
  par(mar = c(3,3,1,1))
  
  test.yr=c(2007:2012)
  
  collect.list = data.frame(rbind(c(2006,"2005-12-31"),c(2007,"2007-12-31"),c(2008,"2008-12-31"),
                                  c(2009,"2009-06-30"),c(2010,"2010-04-30"), c(2011,"2011-07-31"),c(2012,"2012-05-31") ),stringsAsFactors = F) # Dat of isolation
  names(collect.list)=c("year","sample")
  collect.list$year=as.numeric(collect.list$year)
  collect.list$sample=as.Date(collect.list$sample)
  
  flu.isolates=data.frame(read.csv("datasets/Vietnam_H3.csv",stringsAsFactors = F)) # Data from http://www.who.int/influenza/gisrs_laboratory/flunet/en/
  flu.isolates=flu.isolates[flu.isolates$Start_date >= as.Date(collect.list[2,2]), ]
  flu.isolates$Start_date=as.Date(flu.isolates$Start_date)
  flu.isolates$A_H3[is.na(flu.isolates$A_H3)]=0 # Set blank spaces to zero
  
  # Count samples within region
  isolatetab=NULL
  for(ii in 2:length(test.yr) ){ # iterate from 2 as don't know initial period of risk
    isolatetab=c(isolatetab,
                 sum(flu.isolates[flu.isolates$Start_date > collect.list[ii,"sample"] & flu.isolates$Start_date < collect.list[ii+1,"sample"] ,"A_H3"])
    )
  }

  xLims=c(as.Date("2007-12-30"),as.Date("2012-06-01"))

  plot(flu.isolates$Start_date,flu.isolates$total_flu_positive,type="l",xlab="year",ylab="Total isolates detected",ylim=c(0,100),xlim=xLims);title(LETTERS[1],adj=0)
  
  
  plot(flu.isolates$Start_date,flu.isolates$A_H3,type="l",xlab="year",ylab="H3 isolates detected",ylim=c(0,100),xlim=xLims);title(LETTERS[2],adj=0)
  for(ii in 1:length(collect.list$year)){
    lines(c(collect.list[ii,"sample"],collect.list[ii,"sample"]),c(0,100),col="red")
  }
  
  # Plot from 3 as want to line up end points
  plot(collect.list$sample[3:7], isolatetab,ylim=c(0,800),pch=19,col="red",xlab="year",ylab="H3 isolates detected",xlim=xLims);title(LETTERS[3],adj=0)
  

  dev.copy(pdf,paste("plot_simulations/H3Data.pdf",sep=""),width=6,height=6,useDingbats=F)
  dev.off()
  
}


# Plot H3N2 China data

plot_h3_china_reports <- function(){
  
  par(mfrow=c(2,1))
  par(mgp=c(1.8,0.6,0))
  par(mar = c(4,4,1,1))
  
  test.yr=c(2004:2009)
  
  #collect.list = data.frame(rbind(c(2006,"2005-12-31"),c(2007,"2007-12-31"),c(2008,"2008-12-31"),
  #                                c(2009,"2009-06-30"),c(2010,"2010-04-30"), c(2011,"2011-07-31"),c(2012,"2012-05-31") ),stringsAsFactors = F) # Dat of isolation
  collect.list = data.frame(rbind(c(2004,"2004-12-31"),c(2005,"2005-12-31"),c(2006,"2006-12-31"),c(2007,"2007-12-31"),c(2008,"2008-12-31"),
                                  c(2009,"2009-12-31"),c(2010,"2010-12-31"), c(2011,"2011-12-31"),c(2012,"2012-12-31") ),stringsAsFactors = F) # Dat of isolation
  
  names(collect.list)=c("year","sample")
  collect.list$year=as.numeric(collect.list$year)
  collect.list$sample=as.Date(collect.list$sample)
  
  flu.isolates=data.frame(read.csv("datasets/China_H3.csv",stringsAsFactors = F)) # Data from http://www.who.int/influenza/gisrs_laboratory/flunet/en/
  flu.isolates$Start_date=as.Date(flu.isolates$Start_date)
  flu.isolates$A_H3[is.na(flu.isolates$A_H3)]=0 # Set blank spaces to zero
  
  # Count samples within region
  isolatetab=NULL
  for(ii in 2:length(test.yr) ){
    isolatetab=c(isolatetab,
                 sum(flu.isolates[flu.isolates$Start_date > collect.list[ii,"sample"] & flu.isolates$Start_date < collect.list[ii+1,"sample"] ,"A_H3"])
    )
  }
  
  
  xLims=c(as.Date("2004-11-01"),as.Date("2009-06-01"))
  
 # plot(flu.isolates$Start_date,flu.isolates$total_flu_positive,type="l",xlab="Date",ylab="Total isolates detected",ylim=c(0,1000),xlim=xLims);title(LETTERS[1],adj=0)
  
  
  plot(flu.isolates$Start_date,flu.isolates$A_H3,type="l",xlab="Date",ylab="H3 isolates detected",ylim=c(0,2000),xlim=xLims);title(LETTERS[2],adj=0)
  for(ii in 1:length(collect.list$year)){
    lines(c(collect.list[ii,"sample"],collect.list[ii,"sample"]),c(0,10000),col="red")
  }
  
  plot(collect.list$sample[3:7], isolatetab,ylim=c(0,20000),pch=19,col="red",xlab="year",ylab="H3 isolates detected",xlim=xLims);title(LETTERS[3],adj=0)
  
  
  dev.copy(pdf,paste("plot_simulations/H3Data.pdf",sep=""),width=5,height=7,useDingbats=F)
  dev.off()
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run multi-chain diagnostics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


plot.multi.true.vs.estimated<-function(simDat=T,flu.type="H3HN",loadpick=c(1:4),burnCut=0.25,year_test=c(2007:2012),plotmap=F,fr.lim=F,linearFn=T,runsPOST=NULL){
  
  # simDat=T;year_test=c(2007:2012);plotmap=F;fr.lim=T;flu.type="H3HN"; loadpick=c(1:10); burnCut=0.25; loadseed=1; linearFn=T; runsPOST=NULL

  
  col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  # Orange, blue, green, pink
  par(mfrow=c(1,1),mar = c(3,3,1,1),mgp=c(1.8,0.6,0))
  plot(0.1,0.1,pch=19,col="white",xlim=c(0,0.3),ylim=c(0,0.3),xlab="true attack rate",ylab="estimated attack rate", xaxs="i", yaxs="i")
  lines(c(0,1),c(0,1),col='grey')
  
  vals.blank=NULL
  
  for(loadseed in loadpick){
    
    if(simDat==F){loadseedA=paste(loadseed,"_",flu.type,sep="")}else{loadseedA=paste("SIM_",loadseed,sep="")}
    load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseedA,"_lin",linearFn,".RData",sep=""))
    
    hist.sample=length(historytabCollect[,1])/n_part # need this sample value because table is stacked
    
    ind.infN=rowSums(historytabCollect[((round(0.2*hist.sample)*n_part)+1):(hist.sample*n_part),])
    
    
    yob.data=cbind(rep(1,n_part),rep(1,n_part)) # Import age distribution
    n.alive=n_part+0*inf_years
    
    attack=colSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])/(length(ind.infN)*(n.alive/length(yob.data[,1]))) #scale by proportion alive
    attackCI=NULL
    for(jj in 1:length(inf_years)){
      htest <- binom.test(round(n.alive*attack)[jj], n.alive[jj], p = 1,conf.level=0.95)
      meanA=attack[jj]
      conf1=htest$conf.int[1]
      conf2=htest$conf.int[2]
      attackCI=rbind(attackCI,c(meanA,conf1,conf2))
    }
    attackCI=data.frame(attackCI)
    names(attackCI)=c("mean","CI1","CI2")
    load(paste("R_datasets/Simulated_data_",loadseedA,".RData",sep=""))
    attack.yr=colSums(historytabSim)/n_part
    #attack.yr = read.csv(paste("datasets/sim_attackS",loadseed,".csv",sep=""))[,2] # TRUE VALUES

    # - - - - - - - - - - - - - - - - - - 
    # Calculate and plot four fold rise in data
    sconverttab = NULL
    
    for(kk in 2:(length(test.yr)-1) ){ # Only valid for 2008-2011 (no test strains for 2012)
      pyear2=0
      pyear4=0
      nyear=0
      for(ii in 1:n_part){
        t.part1=test.listSim[[ii]][[kk-1]]
        t.part2=test.listSim[[ii]][[kk]]
        
        # Check to match test strains
        matchd1d2 = t.part2[3,]==test.yr[kk]
        
        if(length(matchd1d2) > 0){
          
          diffT = t.part2[2,matchd1d2] - t.part1[2,matchd1d2] # Compare titres
          nyear = nyear +1
          if(median(diffT) >= 2){pyear4 = pyear4 + 1}
          if(median(diffT) >= 1){pyear2 = pyear2 + 1}
        }
        
      }
      sconverttab=rbind(sconverttab, c(pyear4/nyear,pyear2/nyear))
    }
    
    pick_r=match(test_years[2:(length(test.yr)-1)],inf_years)

    for(kk in pick_r){ # Iterate across test years
      points(attack.yr[kk],attackCI$mean[kk],pch=19,col=rgb(1,0,0),cex=0.5)
      lines(c(attack.yr[kk],attack.yr[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=rgb(1,0,0))
    }
    
    points(attack.yr[pick_r], sconverttab[,2],pch=1,cex=1.2,col=rgb(0,0,0))
    points(attack.yr[pick_r], sconverttab[,1],pch=19,cex=1.2,col=rgb(0,0,0))
    
    vals.blank = rbind(vals.blank, cbind(attack.yr[pick_r],attackCI$mean[pick_r],sconverttab) )
    
  }
  
  dev.copy(pdf,paste("plot_simulations/True_vs_rise",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=5,height=4)
  dev.off()
  
  vals.blank0 = vals.blank %>% data.frame()
  names(vals.blank0)=c("true","est","rise4","rise2")
  
  # Plot residuals
  
  breaksN=seq(-0.31,0.31,0.02)
  hist(vals.blank0$est-vals.blank0$true,breaks = breaksN,col=rgb(1,0,0,0.5),border="grey")
  hist(vals.blank0$rise4-vals.blank0$true, breaks = breaksN,add=T,col=rgb(0,0,0,0.2),border=NULL)
  hist(vals.blank0$rise2-vals.blank0$true, breaks = breaksN,add=T,col=rgb(0,1,0,0.2),border=NULL)
  
  par(mar=c(3,3,1,1),mgp=c(2,0.7,0))
  
  plot(density(vals.blank0$est-vals.blank0$true),col="white",frame=T,xaxs="i",yaxs="i",ylab="density",xlab="simulation residual",main="",xlim=c(-0.2,0.2),ylim=c(0,100))
  lines(c(0,0),c(0,100),col="grey")
  lines(density(vals.blank0$est-vals.blank0$true,adjust=20),col="red",lwd=2)
  lines(density(vals.blank0$rise4-vals.blank0$true,adjust=1),col="black",lty=1,lwd=2)
  lines(density(vals.blank0$rise2-vals.blank0$true,adjust=1),col="black",lty=2,lwd=2)
  
  
  dev.copy(pdf,paste("plot_simulations/TrueDens",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""), width=3,height=2)
  dev.off()
  
  
  
}





# Convert map ID tags to strain years

strain_years_convert <- function(){
  
  ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
  
  # Convert antigenic coords into cluster centroids
  strain_years=as.numeric(sapply(ag.coord$viruses,function(x){
    a1=max(which(strsplit(x, "")[[1]]=="/"))
    lstr=nchar(x)
    yr1=substr(x, a1+1, lstr)
    
    if(nchar(yr1)>4){yr1=substr(yr1, 1, 4)}
    year=yr1
    if(nchar(yr1)==2 & as.numeric(yr1)>15){year=paste("19",yr1,sep="")}
    if(nchar(yr1)==2 & as.numeric(yr1)<15){year=paste("20",yr1,sep="")}
    year
  }
  ))
  
}