titre_dependent_boosting_plot <- function(chain, n){
    sampnos <- sample(unique(chain$sampno),n)
    titres <- seq(0,8,by=0.1)
    store <- matrix(nrow=n, ncol=length(titres))
    i <- 1
    for(samp in sampnos){
        pars <- as.numeric(chain[chain$sampno == samp,])
        names(pars) <- colnames(chain)
        mu <- pars["mu"] + pars["mu_short"]

        gradient <- pars["gradient"]
        boost_limit <- pars["boost_limit"]
        boost <- mu*(1-gradient*titres)
        boost[which(titres > boost_limit)] <-  mu*(1-gradient*boost_limit)
        store[i,] <- boost
        i <- i + 1
    }
    range <- apply(store, 2, function(x) quantile(x, c(0.025,0.5,0.975)))
    return(range)
}
