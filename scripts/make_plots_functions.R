
process_infection_times <- function(filename){
  
  titre_dat <- read.csv(paste("data/",filename,"_titre_dat.csv", sep = ""))
  inf_chain  <- data.table::fread(paste("chains/",filename,"_infection_histories.csv",sep=""))
  inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
  
  inf_chain_long <- expand_chain_wrapper(titre_dat,ages,inf_chain,xs=1:length(strain_isolation_times))
  colnames(inf_chain_long) <- c('sampno', 'individual', strain_isolation_times)
  
  #predicted infection times
  predicted <- find_likely_infection(inf_chain_long)
  ind_infections <- predicted$ind_infections

  #true infection times
  inf_hist <- read.csv("data/UK_inf_hist.csv")
  
  #rename for easy plotting
  observed_no <- rowSums(inf_hist)
  predicted_no <- ind_infections[,2]
  
  plot(observed_no, predicted_no,
       xlab = "observed", ylab = "predicted", pch = 19, main = filename,
       ylim = c (0, max(observed_no, predicted_no)) )
  abline(a = 0, b = 1, col = "gray") 
  
  barplot(table(predicted_no,observed_no), legend = TRUE,
          xlab = "observed number of infections",
          args.legend=list(title="predicted"))
  
  
}

find_likely_infection<-function(inf_chain){
  
  #extract ids
  ids <- ages[, 1]

  #storage
  ind_infections<-matrix(nrow=length(ids),ncol=2)
  inf_times_list<-vector("list",length(ids))
  
  for(i in 1:length(ids)){
    #subset by indivivdual 
    infections <- inf_chain[which(inf_chain[,2]==ids[i]),-c(1,2)]
    #how many infections 
    no_infs<-rowSums(infections)
    #median across all posterios samples
    med_no_infs<-median(no_infs)
    #which infections had the highest posterior density
    which_infections<-colSums(infections)/sum(colSums(infections))
    #what was the most likely infection time
    inf_times<-head(order(which_infections,decreasing = T),n=med_no_infs)
    #store results
    inf_times_list[[i]] <- inf_times
    ind_infections[i,]<-c(ids[i],med_no_infs)
  }
  
  return(list(inf_times_list = inf_times_list, ind_infections = ind_infections))
}

make_plots <- function(filename){
  
  titre_dat <- read.csv(paste("data/",filename,"_titre_dat.csv", sep = ""))
  
  chain <- read.csv(paste("chains/",filename,"_chain.csv",sep=""))
  chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
  ## attack rate plot
  inf_chain  <- data.table::fread(paste("chains/",filename,"_infection_histories.csv",sep=""))
  inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
  p<- plot_attack_rates(infection_histories = inf_chain, merge(titre_dat,ages), strain_isolation_times) 
  p <- p + ggtitle(paste(filename)) + geom_point(data = attack_rates_df, aes(y = AR, x = inf_years))
  ggsave(paste(filename,".png",sep=""),p, width = 10, height = 6)
  
  png(paste(filename, "infection_times.png", sep =""),units = "px", res = 200, width = 2000, height = 1000)
  par(mfrow = c(1,2))
  process_infection_times(filename )
  dev.off()
  
  ## antibody parameter plot
  png(paste(filename, "posterior_pars.png") ,units = "px", res = 200, width = 2000, height = 2000)
  par(mfrow=c(3,3))
  #mu
  min_val <- min(chain[, "mu"])
  max_val <- max(chain[, "mu"])
  hist(chain[, "mu"], xlab = expression(paste('long boost (',mu[1],')')), freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  #a
  min_val <- min(chain[, "a"])
  max_val <- max(chain[, "a"])
  hist(chain[, "a"], xlab = "a", freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  
  #mu short
  min_val <- min(chain[, "mu"] / chain[, "a"])
  max_val <- max(chain[, "mu"] / chain[, "a"])
  hist(chain[, "mu"]/chain[, "a"], xlab = expression(paste('short boost (',mu[2],')')), freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  
  #sigma 1
  min_val <- min(chain[, "sigma1"])
  max_val <- max(chain[, "sigma1"])
  hist(chain[, "sigma1"], xlab = expression(paste('long cross (',sigma[1],')')), freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  #b
  min_val <- min(chain[, "b"])
  max_val <- max(chain[, "b"])
  hist(chain[, "b"], xlab = "b", freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  
  
  #sigma 2
  min_val <- min(chain[, "sigma1"] * chain[, "b"])
  max_val <- max(chain[, "sigma1"] * chain[, "b"])
  hist(chain[, "sigma1"] * chain[, "b"], xlab = expression(paste('short cross (',sigma[2],')')), freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  
  #waning
  # min_val <- min(chain[, "wane"])
  # max_val <- max(chain[, "wane"])
  # hist(chain[, "wane"], xlab = expression(paste('waning (',omega,')')), freq = F, main = "",
  #      xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  # 
  #error
  min_val <- min(chain[, "error"])
  max_val <- max(chain[, "error"])
  hist(chain[, "error"], xlab = expression(paste('error (',epsilon,')')), freq = F, main = "",
       xlim = c(min_val, max_val), breaks=seq(min_val,max_val,l=20),col=rgb(0,0.5,0.5,alpha=0.5))
  
  plot(chain$likelihood, type = "l")
  dev.off()
  

  
}



expand_chain_wrapper<-function(titreDat,ages,infChain,xs){
  # Extract numebr of indivs
  n_indiv<-dim(ages)[1]
  # Split vector by 10
  id_seq<-seq(1,n_indiv,by=5)
  if(tail(id_seq,n=1)!=n_indiv) id_seq<-c(id_seq,n_indiv)
  # How many strains
  n_strain<-length(xs)
  # Create and store empty table
  tmp_table <- matrix(NA, nrow=0,ncol=n_strain+2)
  ## Write starting infection histories
  write.table(tmp_table, "infChain_long.csv", row.names=FALSE, col.names=TRUE, sep=",",append=FALSE)
  # for each subset of ids, expand infection history file and write it 
  for(k in 1:(length(id_seq)-1)){
    infChain_sub<-infChain[which(infChain$i>=id_seq[k]&infChain$i<id_seq[k+1]),]
    
    historyTab<-data.frame(expand_summary_infChain(infChain_sub,i_vec=id_seq[k]:(id_seq[k+1]-1),j_vec=1:n_strain))
    if(any(is.na(historyTab))) print(dim(titreDat))
    write.table(historyTab, file="infChain_long.csv",
                col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)

  }
  #final one
  infChain_sub<-infChain[which(infChain$i>=id_seq[k]&infChain$i<=id_seq[k+1]),]
  
  historyTab<-data.frame(expand_summary_infChain(infChain_sub,j_vec=1:n_strain,i_vec=id_seq[k]:(id_seq[k+1])))
  
  if(any(is.na(historyTab))) print(dim(titreDat))
  write.table(historyTab, file="infChain_long.csv",
              col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
  
  
  #then read in infection history file
  infChain <- read.csv("infChain_long.csv")
  colnames(infChain)[1:2]<-c('sampno', 'individual')
  
  #return attack rate 
  return(infChain)
}