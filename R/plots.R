#' Generates all output plots
#'
#' Given an output file location, locations of output files from the MCMC procedure and some control parameters, generates and saves infection history plots, posterior plots (trace and density), attack rates and total number of infections.
#' @param outputDir the file location to save output to
#' @param adaptive_period specifies how many iterations to discard as burning
#' @param chainFile the full file location of the MCMC chain output for theta, the parameter vector
#' @param infHistFile the full file location of the MCMC chain output for infection histories
#' @param titreDat the data frame of titre data used in the MCMC
#' @param antigenicMap the unmelted antigenic map specifying the relationship between strains
#' @param parTab the full parameter control table as used in the MCMC procedure
#' @param ages a data frame with columns for individual and date of birth (titled DOB), where the date of birth is in the same time resolution as the model (ie. DOB in years would be 1991, in months would be 1991*12)
#' @param nIndiv specifies how many random individuals to sample for the infection history plot
#' @param nSamp how many samples from the posterior to take when generating the plot
#' @param filename the actual filename is appended to this part. For example, "test" would give "test_plots.svg", "test_ar.svg" etc
#' @param thin Optional parameter to thin the infection history chain for speed
#' @return nothing
#' @export
generate_all_plots <- function(outputDir, adaptive_period, chainFile, infHistFile,
                               titreDat, antigenicMap, parTab, ages, nIndiv=10,nSamp=1000,
                               filename="test",thin=100){

    ## Get all possible strains
    strainIsolationTimes <- unique(antigenicMap$inf_years)

    ## Read in MCMC chains and discard burnin
    chain <- read.csv(chainFile)
    chain <- chain[chain$sampno > adaptive_period,]
    infectionHistories <- data.table::fread(infHistFile,data.table=FALSE)
    infectionHistories <- infectionHistories[infectionHistories$sampno > adaptive_period,]
    infectionHistories <- infectionHistories[complete.cases(infectionHistories),]

    ## Thin infection history chain 
    sampnos <- unique(infectionHistories$sampno)
    sampnos <- sampnos[seq(1,length(sampnos),by=thin)]

    ## Take subset of individuals
    infectionHistories <- infectionHistories[infectionHistories$sampno %in% sampnos,]
    individuals <- sample(unique(titreDat$individual),nIndiv)

    ## Infection history plot
    p1 <- plot_infection_histories(chain, infectionHistories,titreDat, individuals, antigenicMap,
                                   ages, parTab, nSamp)
    ## Posterior distribution
    posteriors <- plot_posteriors(chain,parTab,TRUE,TRUE,TRUE,TRUE,filename)
    
    xs <- min(strainIsolationTimes):max(strainIsolationTimes)

    ## Attack rates for each year of circulation
    arP <- plot_attack_rates(infectionHistories, titreDat,ages,xs)

    ## Distribution of total number of infections by individual
    infP <- plot_number_infections(infectionHistories)

    ## Save to SVG
    to.svg(print(p1),paste0(filename,"_plots.svg"))
    to.svg(print(arP),paste0(filename,"_ar.svg"))
    to.svg(print(infP),paste0(filename,"_infP.svg"))
}


#' @export
generate_quantiles <- function(x, sigF=3, qs=c(0.025,0.5,0.975),asText=TRUE){
    res <- signif(quantile(x, qs),sigF)
    if(asText){
        res <- paste(res[2]," (",res[1],"-",res[3],")",sep="")
    }
    return(res)
}


#' Generate titre credible intervals
#'
#' Generates credible intervals on titres and infection histories from an MCMC chain output.
#' @param chain the full MCMC chain to generate titre trajectories from
#' @param infectionHistories the MCMC chain for infection histories
#' @param dat the data frame of titre data
#' @param individuals the subset of individuals to generate credible intervals for
#' @param antigenicMap the unmelted antigenic map
#' @param ages the data frame of ages for each individual, with columns for individual and DOB (date of birth)
#' @param parTab the table controlling the parameters in the MCMC chain
#' @param nsamp number of samples to take from posterior
#' @return a list with the titre predictions (95% credible intervals, median and multivariate posterior mode) and the probabilities of infection for each individual in each epoch
#' @export
get_titre_predictions <- function(chain, infectionHistories, dat,
                                  individuals, antigenicMap,
                                  ages, parTab,
                                  nsamp=100){

    ## Take subset of individuals
    dat <- dat[dat$individual %in% individuals,]
    infectionHistories <- infectionHistories[infectionHistories$individual %in% individuals,]
    ages <- ages[ages$individual %in% individuals,]

    ## Format the antigenic map to solve the model
    antigenicMap <- antigenicMap[1:ncol(infectionHistories[,!(colnames(infectionHistories) %in% c("sampno","individual"))]),]
    strainIsolationTimes <- unique(antigenicMap$inf_years)
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))
    
    nstrain <- length(unique(antigenicMap$inf_years))

    ## Empty data structures to save output to
    allRes <- NULL
    infHistDens <- NULL
    refTime <- max(dat$samples,na.rm=TRUE)

    ## Need to align the iterations of the two MCMC chains
    ## and choose some random samples
    tmpSamp <- sample(intersect(unique(infectionHistories$sampno),unique(chain$sampno)),nsamp)

    ## See the function in posteriors.R
    f <- create_model_func(parTab,dat,antigenicMap)
    predicted_titres <- matrix(nrow=nrow(dat),ncol=nsamp)

    ## For each sample, take values for theta and infection histories and simulate titres
    for(i in 1:nsamp){
        pars <- get_index_pars(chain, which(chain$sampno == tmpSamp[i]))
        tmpInfHist <- as.matrix(infectionHistories[infectionHistories$sampno == tmpSamp[i], 1:(ncol(infectionHistories)-2)])
        predicted_titres[,i] <- f(pars, tmpInfHist)
    }

    ## Get 95% credible interval and means
    dat2 <- t(apply(predicted_titres,1, function(x) quantile(x, c(0.025,0.5,0.975))))

    tmpInfChain <- infectionHistories[infectionHistories$sampno %in% tmpSamp,]

    ## Get infection history density for each individual and each epoch
    infHistDens <- plyr::ddply(tmpInfChain, ~individual, function(x) colSums(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
    colnames(infHistDens) <- c("individual",strainIsolationTimes)
    infHistDens <- reshape2::melt(infHistDens, id.vars="individual")
    infHistDens$variable <- as.integer(as.character(infHistDens$variable))
    infHistFinal <- NULL

    ## For each individual, get density for the probability that an epoch was an infection time
    ## The point of the following loop is to mask the densities where infection epochs were either
    ## before an individual was born or after the time that a blood sample was taken
    for(indiv in unique(infHistDens$individual)){
        sampleTimes <- unique(dat[dat$individual==indiv,"samples"])
        tmp <- NULL
        age <- ages[ages$individual == indiv, "DOB"]
        for(samp in sampleTimes){
            indivInfHist <- infHistDens[infHistDens$individual == indiv,]
            indivInfHist[indivInfHist$variable > samp, "value"] <- 0
            indivInfHist <- cbind(indivInfHist, "samples"=samp)
            tmp <- rbind(tmp, indivInfHist)
        }
        infHistFinal <- rbind(infHistFinal, tmp)
        
    }

    ## Find multivariate posterior mode estimate from the chain
    samps <- intersect(unique(infectionHistories$sampno), unique(chain$sampno))
    tmpChain <- chain[chain$sampno %in% samps,]
    tmpInfChain1 <- infectionHistories[infectionHistories$sampno %in% samps,]

    bestPars <- get_best_pars(tmpChain)
    bestI <- tmpChain$sampno[which.max(tmpChain$lnlike)]
    bestInf <- as.matrix(tmpInfChain1[tmpInfChain1$sampno == bestI,1:(ncol(tmpInfChain1)-2)])

    ## Generate trajectory for best parameters
    bestTraj <- f(bestPars, bestInf)
    dat2 <- as.data.frame(dat2)
    colnames(dat2) <- c("lower","median","upper")
    dat2$max <- bestTraj
    dat2 <- cbind(dat,dat2)
    
    return(list("predictions"=dat2,"histories"=infHistFinal))
}

#' Plots infection histories and titre model fits
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed titres and model predicted trajectories
#' @param chain the full MCMC chain to generate titre trajectories from
#' @param infectionHistories the MCMC chain for infection histories
#' @param dat the data frame of titre data
#' @param individuals the subset of individuals to generate credible intervals for
#' @param antigenicMap the unmelted antigenic map
#' @param ages the data frame of ages for each individual, with columns for individual and DOB (date of birth)
#' @param parTab the table controlling the parameters in the MCMC chain
#' @param nsamp number of samples to take from posterior
#' @return a ggplot2 object
#' @export
plot_infection_histories <- function(chain, infectionHistories, dat,
                                     individuals, antigenicMap,ages,parTab,
                                     nsamp=100){
    tmp <- get_titre_predictions(chain, infectionHistories,dat, individuals, antigenicMap,ages, parTab, nsamp)
    dens <- tmp[[1]]
    infHist <- tmp[[2]]
    p <- ggplot(dens) + 
        geom_line(aes(x=virus,y=max),col="blue") + 
        geom_ribbon(aes(x=virus,ymin=lower,ymax=upper),alpha=0.25,fill="blue") + 
        geom_vline(data=infHist,aes(xintercept=variable,alpha=value))+
        geom_point(data=dens,aes(x=virus,y=titre),col="red",size=0.5) + 
        facet_grid(individual~samples) + 
        theme_bw() +
        theme(axis.text.x=element_text(angle=45,hjust=1,size=8))+
        scale_alpha(limits=c(0,1),range=c(0,1))+
        xlab("Year") +
        ylab("Titre") +
        scale_y_continuous(limits=c(0,8),breaks=seq(0,8,by=2))
    p

}

#' @export
plot_posteriors <- function(chain, parTab,
                            calculateESS=TRUE,
                            plotCorr=TRUE,
                            savePlots=FALSE,
                            plotMCMC=TRUE,
                            saveLoc=""){

   
    ## Find best parameters and use only free parameters
    best_pars <- get_best_pars(chain)

    ## Combined chain
    freeChain <- chain[,which(parTab$fixed == 0)+1]
    thinFreeChain <- freeChain[sample(1:nrow(freeChain),1000,replace=T),]

    
    ## Plot correlations
    corP <- NULL
    if(plotCorr){
        corP <- GGally::ggpairs(thinFreeChain) + theme_bw()
        
        if(savePlots){
            to.svg(print(corP), paste0(saveLoc,"corPlot.svg"))
            to.svg(coda::autocorr.plot(freeChain),paste0(saveLoc,"autocorr.svg"))
            #to.pdf(print(corP), paste0(saveLoc,"corPlot.pdf"))
            #to.pdf(coda::autocorr.plot(freeChain),paste0(saveLoc,"autocorr.pdf"))
        }
    }

    ## Get quantiles and create table of results
    thinFreeChain$sigma1drop <- thinFreeChain$mu*thinFreeChain$sigma1
    thinFreeChain$sigma2drop <- thinFreeChain$mu_short*thinFreeChain$sigma2
    thinFreeChain$wane_yr <- thinFreeChain$mu_short*thinFreeChain$wane
    thinFreeChain$errorCorrect1 <- pnorm(3,mean=1.5,sd=thinFreeChain$error)-pnorm(2,mean=2.5,sd=thinFreeChain$error)
    thinFreeChain$errorCorrect2 <- pnorm(4,mean=1.5,sd=thinFreeChain$error)-pnorm(1,mean=2.5,sd=thinFreeChain$error)

    ## Calculate ESS for all free parameters
    all_ess <- NULL
    if(calculateESS){
        all_ess <- coda::effectiveSize(freeChain)
        all_ess <- c(all_ess, rep(NA,5))
    }

    densP <- traceP <- NULL
    if(plotMCMC){
        densP <- bayesplot::mcmc_hist(freeChain) + theme()
        traceP <- bayesplot::mcmc_trace(freeChain) + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
        if(savePlots){
            to.svg(print(densP), paste0(saveLoc,"densities.svg"))
            to.svg(print(traceP), paste0(saveLoc,"traces.svg"))
            #to.pdf(print(densP), paste0(saveLoc,"densities.pdf"))
            #to.pdf(print(traceP), paste0(saveLoc,"traces.pdf"))
        }
    }

        
    results <- apply(thinFreeChain, 2, function(x) generate_quantiles(x))
    allResults <- data.frame(parameter=c(parTab[which(parTab$fixed==0),"names"],
                                         "sigma1drop","sigma2drop","wane_yr","errorCorrect1","errorCorrect2"),
                             estimate=results,
                             ESS=all_ess)
    if(savePlots) write.table(allResults,paste0(saveLoc,"estimates.csv"), row.names=FALSE,sep=",")
    return(list(results=allResults,corPlot=corP,densP=densP,traceP=traceP))    
}

#' Plot historical attack rates
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories
#' @param infectionHistories the MCMC chain for infection histories
#' @param dat the data frame of titre data
#' @param ages the data frame of ages for each individual, with columns for individual and DOB (date of birth)
#' @param yearRange vector of the first and last epoch of potential circulation
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates <- function(infectionHistories, dat, ages, yearRange){
    ## Maximum number of strains
    nstrain <- ncol(infectionHistories[,!(infectionHistories$colnames %in% c("sampno","individual"))])

    ## Find inferred total number of infections from the MCMC output
    tmp <- plyr::ddply(infectionHistories,~sampno,function(x) colSums(x[,1:(ncol(infectionHistories)-2)]))

    ## Scale by number of individuals that were alive in each epoch
    ## and generate quantiles
    n_alive <- sapply(yearRange, function(x) nrow(ages[ages$DOB <= x,]) )
    quantiles <- apply(tmp[,2:ncol(tmp)],2, function(x) quantile(x,c(0.025,0.5,0.975)))
    quantiles <- quantiles/n_alive
    quantiles <- as.data.frame(t(quantiles))
    colnames(quantiles) <- c("lower","median","upper")
    quantiles$year <- yearRange
    quantiles$taken <- quantiles$year %in% unique(dat$samples)

    ## Colour depending on whether or not titres were taken in each year
    quantiles$taken <- ifelse(quantiles$taken,"Yes","No")
    
    p <- ggplot(quantiles) + 
        geom_pointrange(aes(x=year,y=median,ymin=lower,ymax=upper,col=taken)) +
        scale_y_continuous(limits=c(0,1),expand=c(0,0)) +  
        theme_bw() +
        #theme(text=element_text(family="Arial")) +
        ylab("Estimated attack rate") +
        xlab("Year")
    return(p)   
    
}

#' @export
plot_number_infections <- function(infChain){
  n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
  n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
  #n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
  indivHist <- plyr::ddply(n_inf_chain,.(individual),function(x) quantile(x$V1, c(0.025,0.5,0.975)))
    #indivHist <- plyr::ddply(infectionHistories,.(individual),function(x) quantile(rowSums(x),c(0.025,0.5,0.975)))
    colnames(indivHist) <- c("individual","lower","median","upper")
    indivHist <- indivHist[order(indivHist$median),]
    indivHist$individual <- 1:nrow(indivHist)
    p <- ggplot(indivHist) + 
        geom_pointrange(aes(x=individual,y=median,ymin=lower,ymax=upper),
                        size=0.1,shape=21,fatten=0.1) +
        scale_y_continuous(limits=c(0,ncol(infChain)-2),expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0),breaks=seq(0,1100,by=100)) +
        xlab("Individual") +
        ylab("Estimated number of infections") +
        theme_bw()# +
  #theme(text=element_text(family="Arial"))
    return(p)
}


#' @export
plot_data <- function(data, infectionHistories, strainIsolationTimes, n_samps, startInf=NULL){
    indivs <- unique(data$individual)
    infHist <- as.data.frame(cbind(indivs, infectionHistories))
    colnames(infHist) <- c("individual",strainIsolationTimes)
    meltedInfHist <- reshape2::melt(infHist, id.vars="individual")
    meltedInfHist$variable <- as.numeric(as.character(meltedInfHist$variable))
    meltedInfHist <- meltedInfHist[meltedInfHist$value > 0,]
    samps <- sample(unique(data$individual), n_samps)    
    
    p1 <- ggplot(data[data$individual %in% samps,]) +
        geom_line(aes(x=as.integer(virus),y=titre)) +
        #geom_vline(aes(xintercept=samples), col="red",linetype="dashed") +
        geom_vline(data=meltedInfHist[meltedInfHist$individual %in% samps,], aes(xintercept=variable), col="red",linetype="dashed")+
        theme_bw()
    
    if(!is.null(startInf)){
        startInfHist <- as.data.frame(cbind(indivs, startInf))
        colnames(startInfHist) <- c("individual",strainIsolationTimes)
        meltedStartHist <- reshape2::melt(startInfHist, id.vars="individual")
        meltedStartHist$variable <- as.numeric(as.character(meltedStartHist$variable))
        meltedStartHist <- meltedStartHist[meltedStartHist$value > 0,]
        p1 <- p1 + geom_vline(data=meltedStartHist[meltedStartHist$individual %in% samps,], aes(xintercept=variable), col="blue",linetype="dashed")
    }
    p1 <- p1 +
        facet_grid(individual~samples)
    return(p1)
}


#' @export
generate_cumulative_inf_plots <- function(infChainFile, burnin, nIndiv=10, realInfHist=NULL){
    infChain <- data.table::fread(infChainFile,data.table=FALSE)
    infChain <- infChain[infChain$sampno >= burnin,]
    tmp <- plyr::ddply(infChain, ~individual, function(x) colSums(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
    tmp <- reshape2::melt(tmp, id.vars="individual")
    tmp$variable <- as.integer(tmp$variable)

    indivs <- 
    if(!is.null(realInfHist)){
        infHist1 <- as.data.frame(realInfHist)
        infHist1$individual <- 1:n_indiv
        infHist1 <- infHist1[infHist1$individual %in% 1:10,]
   
        colnames(infHist1) <- c(strainIsolationTimes,"individual")
        infHist1 <- reshape2::melt(infHist1, id.vars="individual")
        infHist1 <- infHist1[infHist1$value == 1,]
    }

    ## Generate data to plot the MCMC starting infection history
    firstSamp <- infChain[infChain$sampno == 1,!(colnames(infChain) == "sampno")]
    firstSamp <- reshape2::melt(firstSamp, id.vars="individual")
    firstSamp <- firstSamp[firstSamp$value == 1,]

    ## Generate lower, upper and median cumulative infection histories from the
    ## MCMC chain
    infChain <- infChain[infChain$sampno > mcmcPars["adaptive_period"],]
    histProfiles <- t(apply(infChain, 1, function(x) cumsum(x[1:(ncol(infChain)-2)])))
    histProfiles <- as.data.frame(cbind(histProfiles, "individual"=infChain[,c("individual")]))

    histProfiles_lower <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.025))))
    histProfiles_lower <- reshape2::melt(histProfiles_lower, id.vars="individual")
    colnames(histProfiles_lower) <- c("individual","variable","lower")

    histProfiles_upper <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.975))))
    histProfiles_upper <- reshape2::melt(histProfiles_upper, id.vars="individual")
    colnames(histProfiles_upper) <- c("individual","variable","upper")

    histProfiles_median <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.5))))
    histProfiles_median <- reshape2::melt(histProfiles_median, id.vars="individual")
    colnames(histProfiles_median) <- c("individual","variable","median")

    ## Merge these quantiles into a data frame for plotting
    quantHist <- merge(histProfiles_lower, histProfiles_upper, by=c("individual","variable"))
    quantHist <- merge(quantHist, histProfiles_median, by=c("individual","variable"))

    ## If available, process the real infection history matrix for plotting
    if(!is.null(realInfHist)){
        realHistProfiles <- as.data.frame(t(apply(realInfHist, 1, cumsum)))
        colnames(realHistProfiles) <- strainIsolationTimes
        
        realHistProfiles$individual <- 1:n_indiv
        realHistProfiles <- realHistProfiles[realHistProfiles$individual %in% 1:10,]
        realHistProfiles <- reshape2::melt(realHistProfiles,id.vars="individual")
    }

    ## Process starting point from MCMC chain
    startHistProfiles <- as.data.frame(t(apply(startInf, 1, cumsum)))
    colnames(startHistProfiles) <- strainIsolationTimes
    startHistProfiles$individual <- 1:n_indiv
    startHistProfiles <- startHistProfiles[startHistProfiles$individual %in% 1:10,]
    startHistProfiles <- reshape2::melt(startHistProfiles,id.vars="individual")

    p1 <- ggplot(quantHist[quantHist$individual %in% seq(1,10,by=1),]) + 
        geom_line(aes(x=as.integer(variable),y=median)) + 
        geom_ribbon(aes(x=as.integer(variable),ymin=lower,ymax=upper), alpha=0.2)
    if(!is.null(realInfHist)){
        p1 <- p1 +         
            geom_line(data=realHistProfiles[realHistProfiles$individual %in% seq(1,10,by=1),],aes(x=as.integer(variable),y=value), col="blue")
    }
    p1 <- p1 + 
        geom_line(data=startHistProfiles[startHistProfiles$individual %in% seq(1,10,by=1),],aes(x=as.integer(variable),y=value), col="red") +
        facet_wrap(~individual, scales="free_y")
    return(p1)
}
