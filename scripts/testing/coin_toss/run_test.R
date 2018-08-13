run_test <- function(runID,chainNo,
                     theta_proposal, lambda_proposal, Z_proposal, 
                     sampPropn, yearPropn, 
                     swapPropn=0.5, adaptive=TRUE,
                     pars, coin_probs, coin_results, dat, samps,
                     iter, thin, adapt_period, adapt_freq, plot=FALSE){
    n <- ncol(coin_results)
    indivs <- nrow(coin_results)
    
    newModelDir <- paste(getwd(),"/outputs1/",runID)
    if(!dir.exists(newModelDir)) dir.create(newModelDir,recursive=TRUE)
    
    sampPropn_lookup <- c("A"=1,"B"=0.5,"C"=0.1)
    namesA <- names(sampPropn_lookup)[which(sampPropn_lookup == sampPropn)]
    namesB <- names(sampPropn_lookup)[which(sampPropn_lookup == yearPropn)]
    ## Create filenames to save chains

    topWD <- sprintf("%s/outputs1/%s_%s_%s_%s_%s_%s", getwd(),runID, theta_proposal,lambda_proposal,
                              Z_proposal,namesA,namesB)
    arbit_filename <- sprintf("%s/%s_%s_%s_%s_%s_%s_%d", topWD,runID, theta_proposal,lambda_proposal,
                              Z_proposal,namesA,namesB,chainNo)
    if(!dir.exists(topWD)) dir.create(topWD,recursive=TRUE)
    
    fixed <- c(0,0,0)
    fixed_probs <- rep(0,n)
    covMat_theta <- diag(length(fixed[which(fixed==0)]))
    covMat_probs <- diag(length(fixed_probs[which(fixed_probs==0)]))
    
    data_suggested_coins <- colSums(coin_results)/nrow(coin_results)
    
    real_pars <- pars
    real_coin_results <- coin_results
    real_heads <- colSums(real_coin_results)
    real_heads <- c(real_heads, sum(real_heads))
    
    startPars <- pars
    startPars[1] <- runif(1,0,10)
    startPars[2] <- runif(1,0,1)
    startPars[3] <- runif(1,0,5)
    startPars[which(fixed == 1)] <- pars[which(fixed == 1)]
    startProbs <- runif(n,0,1)
    startProbs[which(fixed_probs == 1)] <- coin_probs[which(fixed_probs == 1)]
    start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE,prob = c(0.9,0.1)),nrow=indivs)
    covMat_theta <- diag(length(fixed[which(fixed==0)]))
    covMat_probs <- diag(length(fixed_probs[which(fixed_probs==0)]))
    
    if(adaptive==FALSE){
        adaptive_period <- 0
    } else {
        adaptive_period <- adapt_period
    }
    
    if(theta_proposal == "univariate_theta"){
        step_theta <- rep(0.1,3)
        step_prob <- rep(0.1,n)
    } else {
        step_theta <- 0.1
        step_prob <- 0.1
    }
    message(arbit_filename)
    if(lambda_proposal == "none"){
        res <- run_MCMC_gibbs(startPars, fixed_pars=fixed,
                              coin_results=start_coins,dat=dat,samps,iter=iter,covMat_theta,thin=thin,step_theta,
                              lower_bounds=c(0,0,0),upper_bounds=c(10,1,5),
                              adapt_freq=adapt_freq,adaptive_period=adaptive_period,printF=1000,temp=1,
                              sampPropn=sampPropn,yearPropn=yearPropn,
                              swapPropn=swapPropn,theta_proposal=theta_proposal,Z_proposal=Z_proposal,
                              alpha=1,beta=1)
    } else {
        res <- run_MCMC_group(startPars, startProbs, 
                              fixed, fixed_probs, 
                              start_coins,dat,samps, iter,
                              lower_bounds=c(0,0,0),upper_bounds=c(10,1,5),
                              covMat_theta, covMat_probs, thin=thin,
                              step_theta=step_theta,step_prob=step_prob,
                              adapt_freq=adapt_freq,adaptive_period=adaptive_period,
                              printF=1000,temp=1,
                              sampPropn=sampPropn,yearPropn=yearPropn,swapPropn=swapPropn,
                              theta_proposal=theta_proposal,Z_proposal=Z_proposal)
    }
    no_infections<- extract_number_infections_from_chain(res[[3]], n, TRUE)
    no_infections <- as.data.frame(no_infections[no_infections[,"sampno"] > adaptive_period,])
    chain <- res[[2]]
    if(lambda_proposal == "none") {
        colnames(chain) <- c("sampno","mu","cr","sd")
    } else {
        colnames(chain) <- c("sampno","mu","cr","sd", 1:n)
    }
    chain <- as.data.frame(chain[chain[,"sampno"] > adaptive_period,])
    chain <- merge(chain, no_infections, by="sampno")
    print(head(chain))
    write.table(chain, paste0(arbit_filename,"_chain.csv"),sep=",",row.names=FALSE)
    ess <- effectiveSize(chain)
    if(plot){
        pdf(paste0(arbit_filename,"_plot.pdf"))
        plot(as.mcmc(chain))
        dev.off()
    }
    return(chain)
}
