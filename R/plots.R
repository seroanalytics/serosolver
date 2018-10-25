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
    infectionHistories <- data.table::fread(infHistFile)
    #infectionHistories <- infectionHistories[infectionHistories$sampno > adaptive_period,]
    infectionHistories <- subset(infectionHistories, sampno > adaptive_period)
    infectionHistories <- infectionHistories[complete.cases(infectionHistories),]
    ## Attack rates for each year of circulation
    xs <- min(strainIsolationTimes):max(strainIsolationTimes)
    arP <- plot_attack_rates(infectionHistories, titreDat,ages,xs)
    
    ## Thin infection history chain 
    sampnos <- unique(infectionHistories$sampno)
    sampnos <- sampnos[seq(1,length(sampnos),by=thin)]

    ## Take subset of individuals
                                        #infectionHistories <- infectionHistories[infectionHistories$sampno %in% sampnos,]
    infectionHistories <- subset(infectionHistories, sampno %in% sampnos)
    individuals <- sample(unique(titreDat$individual),nIndiv)
    infectionHistories <- subset(infectionHistories, i %in% individuals)

    ## Infection history plot
    p1 <- plot_infection_histories(chain, infectionHistories,titreDat, individuals, antigenicMap,
                                   ages, parTab, nSamp)
    ## Posterior distribution
    posteriors <- plot_posteriors(chain,parTab,TRUE,FALSE,TRUE,TRUE,filename)
    
    xs <- min(strainIsolationTimes):max(strainIsolationTimes)

   

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
#' @param addResiduals if true, returns an extra output summarising residuals between the model prediction and data
#' @param mu_indices vector of integers. for random effects on boosting parameter, mu. If random mus are included in the parameter table, this vector specifies which mu to use for each circulation year. For example, if years 1970-1976 have unique boosting, then mu_indices should be c(1,2,3,4,5,6). If every 3 year block shares has a unique boosting parameter, then this should be c(1,1,1,2,2,2)
#' @param measurement_indices default NULL
#' @param for_res_plot TRUE/FALSE value. If using the output of this for plotting of residuals, returns the actual data points rather than summary statistics
#' @return a list with the titre predictions (95% credible intervals, median and multivariate posterior mode) and the probabilities of infection for each individual in each epoch
#' @export
get_titre_predictions <- function(chain, infectionHistories, titreDat,
                                  individuals, antigenicMap,
                                  parTab,
                                  nsamp=100, addResiduals=FALSE,
                                  mu_indices=NULL,
                                  measurement_indices=NULL,
                                  for_res_plot=FALSE){
    ## Need to align the iterations of the two MCMC chains
    ## and choose some random samples
    samps <- intersect(unique(infectionHistories$sampno), unique(chain$sampno))
    chain <- chain[chain$sampno %in% samps,]
    infectionHistories <- infectionHistories[infectionHistories$sampno %in% samps,]

    ## Take subset of individuals
    titreDat <- titreDat[titreDat$individual %in% individuals,]
    infectionHistories <- infectionHistories[infectionHistories$i %in% individuals,]

    titreDat$individual <- match(titreDat$individual, individuals)
    infectionHistories$i <- match(infectionHistories$i, individuals)

    ages <- unique(titreDat[,c("individual","DOB")])    

    ## Format the antigenic map to solve the model
    strainIsolationTimes <- unique(antigenicMap$inf_years)
    nstrain <- length(strainIsolationTimes)
    n_indiv <- length(individuals)
    
    ## Empty data structures to save output to
    allRes <- NULL
    infHistDens <- NULL

    refTime <- max(titreDat$samples,na.rm=TRUE)
 
    
    tmpSamp <- sample(samps, nsamp)

    ## See the function in posteriors.R
    f <- create_posterior_func(parTab,titreDat,antigenicMap, 100,mu_indices=mu_indices,measurement_indices=measurement_indices)
    
    predicted_titres <- residuals <- matrix(nrow=nrow(titreDat),ncol=nsamp)
    samp_record <- numeric(nsamp)
    ## For each sample, take values for theta and infection histories and simulate titres

    for(i in 1:nsamp){
        index <- tmpSamp[i]
        pars <- get_index_pars(chain, which(chain$sampno == index))
        tmpInfHist <- infectionHistories[infectionHistories$sampno == index,]
        tmpInfHist <- as.matrix(Matrix::sparseMatrix(i=tmpInfHist$i, j=tmpInfHist$j, x=tmpInfHist$x,dims=c(n_indiv, nstrain)))
        predicted_titres[,i] <- f(pars, tmpInfHist)
        residuals[,i] <- titreDat$titre - floor(predicted_titres[,i])
        samp_record[i] <- index
    }

    colnames(predicted_titres) <- tmpSamp
    if(for_res_plot) return(list(residuals, samp_record, titreDat,predicted_titres))
    residuals <- cbind(titreDat, residuals)
    ## Get 95% credible interval and means
    dat2 <- t(apply(predicted_titres,1, function(x) quantile(x, c(0.025,0.5,0.975))))
    residuals <- t(apply(residuals,1, function(x) quantile(x, c(0.025,0.5,0.975))))
    residuals <- cbind(titreDat, residuals)
    ## Find multivariate posterior mode estimate from the chain
    bestPars <- get_best_pars(chain)
    bestI <- chain$sampno[which.max(chain$lnlike)]
    bestInf <- infectionHistories[infectionHistories$sampno == bestI,]
    bestInf <- as.matrix(Matrix::sparseMatrix(i=bestInf$i, j=bestInf$j, x=bestInf$x,dims=c(n_indiv, nstrain)))
    ## Generate trajectory for best parameters
    bestTraj <- f(bestPars, bestInf)
    bestResiduals <- titreDat$titre - floor(bestTraj)
    bestResiduals <- cbind(titreDat, bestResiduals, "sampno"=bestI)
    dat2 <- as.data.frame(dat2)
    colnames(dat2) <- c("lower","median","upper")
    dat2$max <- bestTraj
    dat2[dat2 < 0] <- 0
    dat2 <- cbind(titreDat,dat2)
    
    tmpInfChain <- data.table(subset(infectionHistories, sampno %in% tmpSamp))
    ## Get infection history density for each individual and each epoch
    data.table::setkey(tmpInfChain, "i","j")
    infHistDens <- tmpInfChain[,list(V1=sum(x)/length(tmpSamp)),by=key(tmpInfChain)]
    infHistDens$j <- strainIsolationTimes[infHistDens$j]
    colnames(infHistDens) <- c("individual","variable","value")
    infHistFinal <- NULL

    ## For each individual, get density for the probability that an epoch was an infection time
    ## The point of the following loop is to mask the densities where infection epochs were either
    ## before an individual was born or after the time that a blood sample was taken
    for(indiv in unique(infHistDens$individual)){
        sampleTimes <- unique(titreDat[titreDat$individual==indiv,"samples"])
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
    dat2$individual <- individuals[dat2$individual]
    infHistFinal$individual <- individuals[infHistFinal$individual]
    if(addResiduals){
        result <- list("predictions"=dat2,"histories"=infHistFinal, "residuals"=residuals,"bestRes"=bestResiduals)
    } else {
        result <- list("predictions"=dat2,"histories"=infHistFinal)
    }
    return(result)
}

#' Plots infection histories and titre model fits
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed titres and model predicted trajectories
#' @param chain the full MCMC chain to generate titre trajectories from
#' @param infectionHistories the MCMC chain for infection histories
#' @param titreDat the data frame of titre data
#' @param individuals the subset of individuals to generate credible intervals for
#' @param antigenicMap the unmelted antigenic map
#' @param ages the data frame of ages for each individual, with columns for individual and DOB (date of birth)
#' @param parTab the table controlling the parameters in the MCMC chain
#' @param nsamp number of samples to take from posterior
#' @param mu_indices if random effects on boosting parameter, mu, this specifies which entry in the parameter table corresponds to which year. See \code{\link{run_MCMC}}
#' @param measurement_indices default NULL, optional vector giving the index of `measurement_bias` that each strain uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clus
#' @return a ggplot2 object
#' @export
plot_infection_histories <- function(chain, infectionHistories, titreDat,
                                     individuals, antigenicMap,parTab,
                                     nsamp=100,
                                     mu_indices=NULL,
                                     measurement_indices=NULL){
    ages <- unique(titreDat[,c("individual","DOB")])
    individuals <- individuals[order(individuals)]

    ## Generate titre predictions
    tmp <- get_titre_predictions(chain, infectionHistories,titreDat, individuals,
                                 antigenicMap, parTab, nsamp, FALSE,mu_indices,
                                 measurement_indices)

    ## Use these titre predictions and summary statistics on infection histories
    dens <- tmp[[1]]
    infHist <- tmp[[2]]
    p <- ggplot(dens) + 
        geom_line(aes(x=virus,y=median),col="blue") +
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

#' Plot historical attack rates monthly
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories for monthly. The main difference compared to the normal attack rate plot is that pointrange plots don't make as much sense at a very fine time resolution.
#' @param infectionHistories the MCMC chain for infection histories
#' @param dat the data frame of titre data
#' @param ages the data frame of ages for each individual, with columns for individual and DOB (date of birth)
#' @param yearRange vector of the first and last epoch of potential circulation
#' @param n_alive vector with the number of people alive in each year of circulation. Can be left as NULL, and ages will be used to infer this
#' @param ymax Numeric. the maximum y value to put on the axis
#' @param buckets Integer. How many buckets of time is each year split into? ie. 12 for monthly data, 4 for quarterly etc.
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates_monthly<- function(infectionHistories, dat, ages, yearRange,n_alive=NULL,ymax=0.1, buckets=1){
    ## Find inferred total number of infections from the MCMC output
    ##tmp <- plyr::ddply(infectionHistories,~sampno,function(x) colSums(x[,1:(ncol(infectionHistories)-2)]))
    ## Scale by number of individuals that were alive in each epoch
    ## and generate quantiles
    months <- 1:buckets
    years <- range(floor(yearRange/buckets))
    years <- years[1]:years[2]
    labels <- c(sapply(years, function(x) paste0(months, "/",x)))
    labels1 <- labels[1:length(yearRange)]
    labels1 <- labels1[seq(1,length(labels1),by=buckets)]
    yearBreak <- yearRange[seq(1,length(yearRange),by=buckets)]
    if(is.null(n_alive)){
      if(is.null(ages)) n_alive <- length(unique(dat$individual))
      else n_alive <- sapply(yearRange, function(x) nrow(ages[ages$DOB <= x,]) )
    }
    

    ##quantiles <- apply(tmp[,2:ncol(tmp)],2, function(x) quantile(x,c(0.025,0.5,0.975)))
    data.table::setkey(infectionHistories, "sampno","j")
    tmp <- infectionHistories[,list(V1=sum(x)),by=key(infectionHistories)]
    ##tmp <- ddply(infectionHistories, c("sampno","j"), function(x) sum(x$x))
    quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025,0.5,0.975)))
    colnames(quantiles) <- c("j","lower","median","upper")
    quantiles[c("lower","median","upper")] <- quantiles[c("lower","median","upper")]/n_alive
    ##quantiles <- as.data.frame(t(quantiles))
    quantiles$year <- yearRange[quantiles$j]
    quantiles$taken <- quantiles$year %in% unique(dat$samples)

    ## Colour depending on whether or not titres were taken in each year
    quantiles$taken <- ifelse(quantiles$taken,"Yes","No")
    
    p <- ggplot(quantiles) + 
        geom_ribbon(aes(x=year, ymin=lower,ymax=upper),fill="red",alpha=0.2) +
        geom_line(aes(x=year,y=median),col="red")+
        geom_point(aes(x=year,y=median),col="purple",size=0.5)+
        scale_y_continuous(limits=c(-0.005,ymax),expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0),breaks=yearBreak,labels=labels1)+
        theme_bw() +
        theme(axis.text.x=element_text(angle=45,hjust=1))+
        #theme(text=element_text(family="Arial")) +
        ylab("Estimated monthly per capita incidence") +
        xlab("Date")
    return(p)   
    
}

#' Plot historical attack rates
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories
#' @param infectionHistories the MCMC chain for infection histories
#' @param dat the data frame of titre data
#' @param ages the data frame of ages for each individual, with columns for individual and DOB (date of birth)
#' @param yearRange vector of the first and last epoch of potential circulation
#' @param n_alive vector with the number of people alive in each year of circulation. Can be left as NULL, and ages will be used to infer this
#' @param pointsize Numeric - how big should each point be?
#' @param fatten Numeric - fatten parameter for ggplot pointrange
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates <- function(infectionHistories, dat, ages, yearRange,n_alive=NULL,pointsize=1, fatten=1){
    ## Find inferred total number of infections from the MCMC output
    ##tmp <- plyr::ddply(infectionHistories,~sampno,function(x) colSums(x[,1:(ncol(infectionHistories)-2)]))
    ## Scale by number of individuals that were alive in each epoch
    ## and generate quantiles
    if(is.null(n_alive)){
        if(is.null(ages)) n_alive <- length(unique(dat$individual))
        else n_alive <- sapply(yearRange, function(x) nrow(ages[ages$DOB <= x,]) )
    }

    ##quantiles <- apply(tmp[,2:ncol(tmp)],2, function(x) quantile(x,c(0.025,0.5,0.975)))
    data.table::setkey(infectionHistories, "sampno","j")
    tmp <- infectionHistories[,list(V1=sum(x)),by=key(infectionHistories)]
    ##tmp <- ddply(infectionHistories, c("sampno","j"), function(x) sum(x$x))
    quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025,0.5,0.975)))
    colnames(quantiles) <- c("j","lower","median","upper")
    quantiles[c("lower","median","upper")] <- quantiles[c("lower","median","upper")]/n_alive
    quantiles$year <- yearRange[quantiles$j]
    quantiles$taken <- quantiles$year %in% unique(dat$samples)

    ## Colour depending on whether or not titres were taken in each year
    quantiles$taken <- ifelse(quantiles$taken,"Yes","No")

    p <- ggplot(quantiles) + 
        geom_pointrange(aes(x=year,y=median,ymin=lower,ymax=upper,col=taken),size=pointsize, fatten=fatten) +
        scale_y_continuous(limits=c(-0.1,1),expand=c(0,0)) +  
        theme_bw() +
        ylab("Estimated attack rate") +
        xlab("Year")
    return(p)   
    
}

#' @export
plot_number_infections <- function(infChain){
    ##n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
    ##n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
    n_strain <- max(infChain$j)
    data.table::setkey(infChain, "i","sampno")
    n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
    ##n_inf_chain <- ddply(infChain, c("sampno","i"), function(x) sum(x$x))
    ##n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
    indivHist <- plyr::ddply(n_inf_chain,.(i),function(x) quantile(x$V1, c(0.025,0.5,0.975)))
    ##indivHist <- plyr::ddply(infectionHistories,.(individual),function(x) quantile(rowSums(x),c(0.025,0.5,0.975)))
    colnames(indivHist) <- c("individual","lower","median","upper")
    indivHist <- indivHist[order(indivHist$median),]
    indivHist$individual <- 1:nrow(indivHist)
    p <- ggplot(indivHist) + 
        geom_pointrange(aes(x=individual,y=median,ymin=lower,ymax=upper),
                        size=0.1,shape=21,fatten=0.1) +
        scale_y_continuous(limits=c(0,n_strain),expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0),breaks=seq(0,1100,by=100)) +
        xlab("Individual") +
        ylab("Estimated number of infections") +
        theme_bw()# +
  #theme(text=element_text(family="Arial"))
    return(p)
}

#' Useful plot for looking at simulated data
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
generate_cumulative_inf_plots <- function(infChainFile, burnin, indivs, realInfHist=NULL, startInf=NULL,
                                          strainIsolationTimes, nsamp=100, ages=NULL, numberCol=1){
    infChain <- data.table::fread(infChainFile)
    firstSamp <- infChain[infChain$sampno == min(infChain$sampno),]
    infChain <- infChain[infChain$sampno >= burnin,]

    samps <- sample(unique(infChain$sampno), nsamp)
    infChain <- infChain[infChain$sampno %in% samps,]
    
    infChain1 <- infChain[infChain$i %in% indivs,]
    data.table::setkey(infChain1, "i","j")
    maxSampno <- length(unique(infChain1$sampno))
    densities <- infChain1[,list(V1=sum(x)/maxSampno),by=key(infChain1)]
    densities$j <- as.numeric(strainIsolationTimes[densities$j])
    densities$i <- as.numeric(densities$i)
    allCombos <- data.table(expand.grid(i=indivs,j=strainIsolationTimes))
    allCombos$j <- as.numeric(allCombos$j)
    allCombos$i <- as.numeric(allCombos$i)
    allCombos <- data.table::fsetdiff(allCombos[,c("i","j")], densities[,c("i","j")])
    allCombos$V1 <- 0
    densities <- rbind(allCombos, densities)

    if(!is.null(realInfHist)){
        infHist1 <- as.data.frame(realInfHist)
        infHist1 <- infHist1[indivs,]
        infHist1$individual <- indivs
        
        colnames(infHist1) <- c(strainIsolationTimes,"i")
        infHist1 <- reshape2::melt(infHist1, id.vars="i")
        infHist1$variable <- as.numeric(as.character(infHist1$variable))                                                   
        infHist1 <- infHist1[infHist1$value == 1,]
    }
    density_plot <- ggplot() +
        geom_line(data=densities,aes(x=j,y=V1))
    if(!is.null(realInfHist)){
        density_plot <- density_plot +
            geom_vline(data=infHist1,aes(xintercept=variable),col="red")
    }
    density_plot <- density_plot + 
        facet_wrap(~i) +
        xlab("Year") +
        ylab("Density") +
        theme_bw()
    
    ## Generate lower, upper and median cumulative infection histories from the
    ## MCMC chain
    tmpInfChain <- infChain[infChain$i %in% indivs,]
    histProfiles <- ddply(tmpInfChain, .(i, sampno), function(x){
        empty <- numeric(length(strainIsolationTimes))
        empty[x$j] <- 1
        cumsum(empty)
    })
    
    
    histProfiles <- histProfiles[,colnames(histProfiles) != "sampno"]
    colnames(histProfiles) <- c("i",strainIsolationTimes)
    histProfiles_lower <- ddply(histProfiles, ~i, function(x) apply(x, 2, function(y) quantile(y, c(0.025))))
    histProfiles_lower <- reshape2::melt(histProfiles_lower, id.vars="i")
    colnames(histProfiles_lower) <- c("individual","variable","lower")

    histProfiles_upper <- ddply(histProfiles, ~i, function(x) apply(x, 2, function(y) quantile(y, c(0.975))))
    histProfiles_upper <- reshape2::melt(histProfiles_upper, id.vars="i")
    colnames(histProfiles_upper) <- c("individual","variable","upper")

    histProfiles_median <- ddply(histProfiles, ~i, function(x) apply(x, 2, function(y) quantile(y, c(0.5))))
    histProfiles_median <- reshape2::melt(histProfiles_median, id.vars="i")
    colnames(histProfiles_median) <- c("individual","variable","median")

    ## Merge these quantiles into a data frame for plotting
    quantHist <- merge(histProfiles_lower, histProfiles_upper, by=c("individual","variable"))
    quantHist <- merge(quantHist, histProfiles_median, by=c("individual","variable"))

    ## If available, process the real infection history matrix for plotting
    if(!is.null(realInfHist)){
        realHistProfiles <- as.data.frame(t(apply(realInfHist, 1, cumsum)))
        colnames(realHistProfiles) <- strainIsolationTimes
        
        realHistProfiles <- realHistProfiles[indivs,]
        realHistProfiles$individual <- indivs
        realHistProfiles <- reshape2::melt(realHistProfiles,id.vars="individual")
    }

    ## Process starting point from MCMC chain
    if(!is.null(startInf)){
        startHistProfiles <- as.data.frame(t(apply(startInf, 1, cumsum)))
        colnames(startHistProfiles) <- strainIsolationTimes
        startHistProfiles <- startHistProfiles[indivs,]
        startHistProfiles$individual <- indivs
        startHistProfiles <- reshape2::melt(startHistProfiles,id.vars="individual")
    }
    
    p1 <- ggplot(quantHist[quantHist$individual %in% indivs,]) + 
        geom_line(aes(x=as.integer(as.character(variable)),y=median)) + 
        geom_ribbon(aes(x=as.integer(as.character(variable)),ymin=lower,ymax=upper), alpha=0.2)
       
    if(!is.null(realInfHist)){
        p1 <- p1 +         
            geom_line(data=realHistProfiles[realHistProfiles$individual %in% indivs,],aes(x=as.integer(as.character(variable)),y=value), col="blue")
    }
    if(!is.null(startInf)){
        p1 <- p1 + 
            geom_line(data=startHistProfiles[startHistProfiles$individual %in% indivs,],aes(x=as.integer(as.character(variable)),y=value), col="red")
    }

    if(!is.null(ages)){
        tmpAge <- ages[ages$individual %in% indivs,]
        ageMask <- create_age_mask(tmpAge[,2], strainIsolationTimes)
        ageDat <- data.frame(j=ageMask,individual=indivs[order(indivs)])
        p1 <- p1 + geom_vline(data=ageDat,aes(xintercept=strainIsolationTimes[j]),col="purple",linetype="dashed")
    }
    p1 <- p1 + facet_wrap(~individual, ncol=numberCol) +
        theme_bw() +
        ylab("Cumulative infections") +
        xlab("Circulation time")
    return(list(p1,density_plot))
}
