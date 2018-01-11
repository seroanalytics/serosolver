#' @export
generate_all_plots <- function(outputDir, adaptive_period, chainFile, infHistFile,
                               titreDat, antigenicMap, parTab, ages, nIndiv=10,nSamp=1000,
                               fileName="test"){
    strainIsolationTimes <- unique(antigenicMap$inf_years)
    chain <- read.csv(chainFile)
    infectionHistories <- data.table::fread(infHistFile,data.table=FALSE)
    individuals <- sample(unique(titreDat$individual),nIndiv)
    p1 <- plot_infection_histories(chain, infectionHistories,titreDat, individuals, antigenicMap, ages, nSamp)
    posteriors <- plot_posteriors(chain,parTab,TRUE,TRUE,TRUE,TRUE,fileName)
    xs <- min(strainIsolationTimes):max(strainIsolationTimes)
    arP <- plot_attack_rates(infectionHistories, titreDat,ages,xs)
    infP <- plot_number_infections(infectionHistories)
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


#' Need to mask by age and sample time
#' @export
get_titre_predictions <- function(chain, infectionHistories, dat,
                                  individuals, antigenicMap,
                                  ages,
                                  nsamp=100){
    antigenicMap <- antigenicMap[1:ncol(infectionHistories[,!(colnames(infectionHistories) %in% c("sampno","individual"))]),]
    strainIsolationTimes <- unique(antigenicMap$inf_years)
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))
    nstrain <- length(unique(antigenicMap$inf_years))
    allRes <- NULL
    infHistDens <- NULL
    refTime <- max(dat$samples,na.rm=TRUE)
    
    for(indiv in individuals){
        tmpInfHist <- infectionHistories[infectionHistories$individual==indiv,]
        tmpDat <- dat[dat$individual==indiv,]
        tmpDat <- tmpDat[complete.cases(tmpDat),]
        age <- as.numeric(ages[ages$individual == indiv,"age"])
        birthYear <- refTime - age
        tmpSamp <- sample(intersect(unique(tmpInfHist$sampno),unique(chain$sampno)),nsamp)
        for(sampleTime in unique(tmpDat$samples)){
            predicted_titres <- matrix(nrow=nsamp,ncol=length(strainIsolationTimes))
            tmpIHist <- tmpInfHist
            tmpIHist[,which(strainIsolationTimes > sampleTime)] <- 0
            tmpIHist[,which(strainIsolationTimes < birthYear)] <- 0
            infCounts <- apply(tmpIHist[,1:nstrain],2,function(x) table(factor(x, levels=c(0,1))))
            infSD <- apply(tmpIHist[,1:nstrain],2,sd)
            y <- (infCounts[2,]/colSums(infCounts))
            histDens <- data.frame(year=strainIsolationTimes,density=y)
            histDens$year <- as.numeric(as.character(histDens$year))
            
            for(i in 1:nsamp){
                sampno <-tmpSamp[i]
                index <- which(chain$sampno == sampno)
                infectionHistory <- tmpIHist[tmpIHist$sampno == sampno,1:nstrain]
                pars <- get_index_pars(chain,index)
                antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
                antigenicMapLong[antigenicMapLong < 0] <- 0
                antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
                antigenicMapShort[antigenicMapShort < 0] <- 0
                predicted_titres[i,] <- infection_model_indiv(pars, as.numeric(infectionHistory),
                                                              sampleTime,
                                                              strainIsolationTimes,
                                                              antigenicMapLong,antigenicMapShort)
            }
            dat2 <- t(apply(predicted_titres,2, function(x) quantile(x, c(0.025,0.5,0.975))))
            dat2 <- as.data.frame(cbind(year=strainIsolationTimes,dat2,
                                        max=infection_model_indiv(bestpars, as.numeric(infectionHistory),
                                                    sampleTime,strainIsolationTimes,
                                                    antigenicMapLong,antigenicMapShort)))
            colnames(dat2) <- c("year","lower","median","upper","max")
            dat2$samples <- sampleTime
            dat2$individual <- indiv
            allRes <- rbind(allRes, dat2)
            
            histDens$samples <- sampleTime
            histDens$individual <- indiv
            infHistDens <- rbind(infHistDens,histDens)
        }
    }
    return(list("predictions"=allRes,"histories"=infHistDens))
}

#' @export
plot_infection_histories <- function(chain, infectionHistories, dat,
                                     individuals, antigenicMap,ages,
                                     nsamp=100){
    tmp <- get_titre_predictions(chain, infectionHistories,dat, individuals, antigenicMap,ages, nsamp)
    dens <- tmp[[1]]
    infHist <- tmp[[2]]
    p <- ggplot(dens) + 
        geom_line(aes(x=year,y=median),col="blue") + 
        geom_ribbon(aes(x=year,ymin=lower,ymax=upper),alpha=0.25,fill="blue") + 
        geom_vline(data=infHist,aes(xintercept=year,alpha=density))+
        geom_point(data=dat[dat$individual %in% individuals,],aes(x=virus,y=titre),col="red",size=0.5) + 
        facet_grid(individual~samples) + 
        theme_bw() +
        theme(text=element_text(family="Arial"), axis.text.x=element_text(angle=45,hjust=1,family="Arial",size=8))+
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
        densP <- bayesplot::mcmc_hist(freeChain) + theme(text=element_text(family="Arial"))
        traceP <- bayesplot::mcmc_trace(freeChain) + theme(text=element_text(family="Arial"),
                                                           axis.text.x=element_text(angle=45,hjust=1,size=8,
                                                                                    family="Arial"))
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


#' @export
plot_attack_rates <- function(infectionHistories, dat, ages, yearRange){
    nstrain <- ncol(infectionHistories[,!(infectionHistories$colnames %in% c("sampno","individual"))])
    tmp <- plyr::ddply(infectionHistories,~sampno,function(x) colSums(x[,1:47]))
    ages$birthYear <- max(dat$samples,na.rm=TRUE) - ages$age
    n_alive <- sapply(yearRange, function(x) nrow(ages[ages$birthYear <= x,]) )
    quantiles <- apply(tmp[,2:ncol(tmp)],2, function(x) quantile(x,c(0.025,0.5,0.975)))
    quantiles <- quantiles/n_alive
    quantiles <- as.data.frame(t(quantiles))
    colnames(quantiles) <- c("lower","median","upper")
    quantiles$year <- yearRange
    quantiles$taken <- quantiles$year %in% unique(dat$samples)
    quantiles$taken <- ifelse(quantiles$taken,"Yes","No")
    
    p <- ggplot(quantiles) + 
        geom_pointrange(aes(x=year,y=median,ymin=lower,ymax=upper,col=taken)) +
        scale_y_continuous(limits=c(0,1),expand=c(0,0)) +  
        theme_bw() +
        theme(text=element_text(family="Arial")) +
        ylab("Estimated attack rate") +
        xlab("Year")
    return(p)   
    
}

#' @export
plot_number_infections <- function(infectionHistories){
    indivHist <- plyr::ddply(infectionHistories,.(individual),function(x) quantile(rowSums(x),c(0.025,0.5,0.975)))
    colnames(indivHist) <- c("individual","lower","median","upper")
    indivHist <- indivHist[order(indivHist$median),]
    indivHist$individual <- 1:nrow(indivHist)
    p <- ggplot(indivHist) + 
        geom_pointrange(aes(x=individual,y=median,ymin=lower,ymax=upper),
                        size=0.1,shape=21,fatten=0.1) +
        scale_y_continuous(limits=c(0,40),expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0),breaks=seq(0,1100,by=100)) +
        xlab("Individual") +
        ylab("Estimated number of infections") +
        theme_bw() +
        theme(text=element_text(family="Arial"))
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
        geom_point(aes(x=as.integer(virus),y=titre)) +
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
