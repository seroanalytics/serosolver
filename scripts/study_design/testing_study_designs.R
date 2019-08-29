code.dir <- "~/Documents/Github/serosolver"
setwd(code.dir)
devtools::load_all()

## Specify paramters controlling the MCMC procedure
mcmc_pars <- c("iterations"=500000,"adaptive_period"=500000, "burnin"=50000,
               "hist_sample_prob"=0.5,"thin"=50,"thin_hist"=200,
               "swap_propn"=0.5,"hist_switch_prob"=0.2,"year_swap_propn"=0.5)


## Antigenic map for cross reactivity parameters
antigenic_map <- read.csv(file.path(code.dir,"data/fonville_map_approx.csv"),
                         stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenic_map, buckets=1)

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968,]
strain_isolation_times <- unique(fit_dat$inf_years)

## Read in parameter table
par_tab <- read.csv(file.path(code.dir,"inputs/par_tab_study_design.csv"))


library(doParallel)
library(foreach)

#create cluster with desired number of cores
cl<-makeCluster(4)

#register cluster
registerDoParallel(cl)

# ###mcmc chains
# for(f in 1:3){
#   foreach(i=1:3)%do%{
#     for(type in c("L","CS")){
#       filename <- paste0(i,"study",f,"_b_", type, sep="")
#       
#       titre_dat <- read.csv(paste("data/study_design/",filename,"dat.csv",sep = ""))
#       ages <- read.csv(paste("data/study_design/",filename,"ages.csv",sep = ""))
#       
#      tryCatch( res <- run_MCMC(par_tab = par_tab, titre_dat = merge(titre_dat,ages, by = "individual"), 
#                       antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
#                       mvr_pars = NULL, start_inf_hist = NULL, filename=paste("chains/", filename,sep=""),
#                       CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = NULL,
#                       version = 2,  
#                       fast_version = TRUE),
#       error = function(e) write.table(filename, paste(filename,".txt",sep="")))
#     }
#   }
# }



stopCluster(cl)

# variable attack rate

filename <- "AR_change_L"

titre_dat <- read.csv(paste("data/study_design/",filename,"dat.csv",sep = ""))
ages <- read.csv(paste("data/study_design/",filename,"ages.csv",sep = ""))

res <- run_MCMC(par_tab = par_tab, titre_dat = merge(titre_dat,ages, by = "individual"), 
                antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                mvr_pars = NULL, start_inf_hist = NULL, filename=paste("chains/", filename,sep=""),
                CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = NULL,
                version = 2,  
                fast_version = TRUE)


filename <- "AR_change_b_L"

titre_dat <- read.csv(paste("data/study_design/",filename,"dat.csv",sep = ""))
ages <- read.csv(paste("data/study_design/",filename,"ages.csv",sep = ""))

res <- run_MCMC(par_tab = par_tab, titre_dat = merge(titre_dat,ages, by = "individual"), 
                antigenic_map = fit_dat, mcmc_pars = mcmc_pars,
                mvr_pars = NULL, start_inf_hist = NULL, filename=paste("chains/", filename,sep=""),
                CREATE_POSTERIOR_FUNC = create_posterior_func_fast, CREATE_PRIOR_FUNC = NULL,
                version = 2,  
                fast_version = TRUE)


# main figures

which_pars <- c("mu", "mu_short", "wane","sigma1","sigma2","error")
true_val <- par_tab$values[par_tab$names%in%which_pars]

res_master <- vector("list", 3)

for(f in 1:3){
  res_list <- vector("list", 6)
  k <- 1
  for(type in c("CS", "L")){
    for(i in 1:3){
      filename <- paste0(i, "study", f, "_", type, sep="")
      print(filename)
      
      res_list[[k]] <- vector("list", 3)
      titre_dat <- read.csv(paste("data/study_design/",filename,"dat.csv",sep = ""))
      ages <- read.csv(paste("data/study_design/",filename,"ages.csv",sep = ""))
      n_indiv <- dim(ages)[1]

      inf_chain <- data.table::fread(paste("chains/",filename,"_infection_histories.csv",sep=""))
      inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
      
      chain <- read.csv(paste("chains/", filename,"_chain.csv",sep=""),stringsAsFactors=FALSE)
      chain <- chain[chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]
      chain <- chain[,c('mu', 'mu_short','wane','sigma1','sigma2','error','likelihood')]
  
      DOBs <- ages$DOB
      age_mask <- create_age_mask(DOBs, strain_isolation_times)
      strain_mask <- create_strain_mask(titre_dat,strain_isolation_times)
      masks <- data.frame(cbind(age_mask, strain_mask))
      n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x) nrow(masks[masks$age_mask <=x & masks$strain_mask >= x,]))  

      p <- plot_attack_rates(infection_histories = inf_chain, merge(titre_dat, ages), strain_isolation_times, by_val = 10)
      pg <- ggplot_build(p)
      AR <- pg$data[[1]]
      

      mat <- matrix(nrow = dim(chain)[1], ncol=6)

      for(j in 1:6){
        chain_adj <-chain[,which(colnames(chain)==which_pars[j])]
        chain_adj <-(chain_adj - true_val[j]) / true_val[j]
        mat[,j] <- chain_adj
      }
      
      res_list[[k]][[1]] <- AR
      res_list[[k]][[2]] <- chain
      res_list[[k]][[3]] <- mat
      
      
      k <- k + 1
    }
  }
  res_master[[f]] <- res_list
}

png("Figure2.png", width = 3300, height = 2000, res = 300, units = "px")
age_names <- c("ages 2 - 4", "ages 40 - 75", "ages 2 - 75")
m <- matrix(c(1, 1, 1, 2,3,4,5,6,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.1,0.4,0.4))
par(mar = c(0,1,1,1), cex.axis=1.5, cex.main=1.5,cex.lab=1.5)
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", c("1 strain/year", "2 strain/year","4 strain/year","Cross-sectional","Longtitudinal"),
       col = c(1, "orange", "red", "black" ,"black"), lty = c(1, 1, 1, NA, NA), 
       pch = c(NA, NA, NA, 19, 17), bty = "n", ncol = 1, horiz = TRUE ,cex = 1.5)

for(f in 1:3){
  par(mar = c(5, 4, 4, 2) + 0.1)
  res <- res_master[[f]]
  
  x_lim <- c(2011.5, 2015.5)
  
  adj <- 0.2
  increment <- 0.05
  current_res <- res[[1]][[1]]
  plotCI(current_res$x - adj - increment, current_res$y, ui=current_res$ymax ,li=current_res$ymin, 
         ylab = "Attack rate", xlab = "Virus year" , ylim = c(0, 1),
         main=paste(age_names[f]), xaxt = "n", pch = 19, xlim = x_lim, col = "black")
  current_res<-res[[2]][[1]]
  plotCI(current_res$x - adj, current_res$y, ui=current_res$ymax ,li=current_res$ymin, 
       pch = 19, add = T, col="orange")
  current_res<-res[[3]][[1]]
  plotCI(current_res$x -  adj + increment, current_res$y, ui=current_res$ymax ,li=current_res$ymin, 
         pch = 19, add = T, col="red")
  axis(1, at = seq(2012, 2015, by = 1), labels = c("t-3", "t-2","t-1", "t"))
  
  current_res<-res[[4]][[1]]
  plotCI(current_res$x  + increment*2 , current_res$y, ui=current_res$ymax ,li=current_res$ymin, 
         pch = 17, add = T, col="black")
  current_res<-res[[5]][[1]]
  plotCI(current_res$x  + increment*3, current_res$y, ui=current_res$ymax ,li=current_res$ymin, 
         pch = 17, add = T, col="orange")
  current_res<-res[[6]][[1]]
  plotCI(current_res$x  + increment*4, current_res$y, ui=current_res$ymax ,li=current_res$ymin, 
         pch = 17, add = T, col="red")
  abline(h = 0.15, col = "gray", lty = "dashed")
  axis(1, at = seq(2012, 2015, by = 1), labels = c("t-3", "t-2","t-1", "t"))
}


res <- res_master[[1]]
current_res <- res[[6]][[1]]
plotCI(current_res$x-0.12, current_res$y, ui = current_res$ymax, li=current_res$ymin,
       ylab = "Attack rate", xlab = "Virus year", ylim = c(0, 0.6),
       main = paste("Longitudinal, ", age_names[1]), xaxt = "n", 
       pch=17, xlim = x_lim, col = "red")
abline(h=0.15, col = "gray")
axis(1, at = seq(2012 ,2015, by = 1),labels = c("t-3", "t-2","t-1", "t"))
for(f in 2:3){
  res <- res_master[[f]]
  current_res <- res[[6]][[1]]
  plotCI(current_res$x-0.12, current_res$y, ui = current_res$ymax, li=current_res$ymin,
         ylab = "Attack rate", xlab = "Virus year", ylim = c(0, 0.6),
         main = paste("Longitudinal, ", age_names[f]), xaxt = "n", 
         pch=17, xlim = x_lim, col = "red")
  
  abline(h = 0.15, col = "gray", lty = "dashed")
  axis(1, at = seq(2012 ,2015, by = 1),labels = c("t-3", "t-2","t-1", "t")) 
}
dev.off()

#Figure 3
col_vec <- c("black", "orange", "red", "black", "orange", "red")
type_vec <-c (1, 1, 1, 2, 2, 2)
pch_vec <- c(rep(19, 3), rep(17, 3))

png("Figure3.png", width = 3000, height = 1600, res = 300, units = "px")
m <- matrix(c(1, 1, 2,3,4,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.1,0.4,0.4))
par(mar = c(0,1,1,1), cex.axis=1.5, cex.main=1.5,cex.lab=1.5)
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", c("1 strain/year", "2 strain/year","4 strain/year","Cross-sectional","Longtitudinal"),
       col = c(1, "orange", "red", "black" ,"black"), lty = c(1, 1, 1, NA, NA), 
       pch = c(NA, NA, NA, 19, 17), bty = "n", ncol = 1, horiz = TRUE , cex = 1.5)

par(mar = c(3, 4, 4, 2) + 0.1)

titles <- c( expression(paste('wane (',omega,'), ages 2 - 4')),
             expression(paste('long cross (',sigma[1],'), ages 2 - 4')),
             expression(paste('wane (',omega,'), ages 40 - 75')),
             expression(paste('long cross (',sigma[1],'), ages 40 - 75')))
k <- 1
for(f in 1:2){
  res <- res_master[[f]]
  for(which_par in c("wane", "sigma1")){
    which_col <- which(which_pars == which_par)
    tmp_mat<- matrix(ncol = 3, nrow = 6)
    for(i in 1:6){
      res_chain <- res[[i]][[3]]
      tmp_mat[i, ] <- quantile(res_chain[,which_col], c(0.025,0.5,0.975))
    }
    
    plotCI(tmp_mat[, 2],ui = tmp_mat[, 3],li = tmp_mat[, 1],
           ylab="Relative error",xlab = "", xaxt = "n", col = col_vec,
           pch = pch_vec, main = titles[k])
    abline(h=0, col = "gray", lty = "dashed")
    
    abline(v = 3.5)
    k <- k +1
  }
}
dev.off()


# Figure 4
titles <- c( 'ages 2 - 4','ages 40 - 75')

labels_vec <- c(expression(paste('long boost, ',mu[1])),
                expression(paste('short boost, ',mu[2])),
                expression(paste('wane, ',omega)),
                expression(paste('long cross, ',sigma[1])),
                expression(paste('short cross, ',sigma[2])),
                expression(paste('error, ',epsilon)))
short_indicies <- c(2, 3, 5)
long_indicies <- c(1, 4)

png("Figure4a.png", width = 3000, height = 3000, res = 300, units = "px")
par(mfrow = c(2, 2), cex = 1.1,mar = c(5, 4, 4, 2) + 0.1)
tmp_list <- vector("list", 2)
for(f in 1:2){
  res <- res_master[[f]]
  res_chain <- res[[6]][[3]]
  tmp_mat <-  sapply(1:6, function(j) quantile(res_chain[, j], c(0.025, 0.5, 0.975))) 
  tmp_list[[f]] <- tmp_mat
  
  plotCI(tmp_mat[2,short_indicies ],ui = tmp_mat[3,short_indicies],li = tmp_mat[1, short_indicies],
         ylab="Relative error",xlab = "", col = "red",
         pch = 17, xaxt="n", main = titles[f], xlim = c(0, 4))
  text(seq(1, length(short_indicies),by = 1),  par("usr")[3] - 0.15, srt = 30, adj = 1,labels = labels_vec[short_indicies], xpd = TRUE)
  abline(h=0, col = "gray",lty = "dashed")
  
  plotCI(c(1, 2), tmp_mat[2,long_indicies ],ui = tmp_mat[3,long_indicies],li = tmp_mat[1, long_indicies],
         ylab="Relative error",xlab = "", col = "red",
         pch = 17, xaxt="n", main = titles[f], xlim = c(0, 3))
  text(seq(1, length(long_indicies),by = 1),  par("usr")[3] - 0.05, srt = 30, adj = 1,labels = labels_vec[long_indicies], xpd = TRUE)
  abline(h=0, col = "gray", lty = "dashed")
  
}
dev.off()

png("Figure4b.png", width = 3000, height = 3000, res = 300, units = "px")
par(mfrow = c(2, 2), cex = 1.1,mar = c(5, 4, 4, 2) + 0.1)
tmp_list <- vector("list", 2)
for(f in 1:2){
  res <- res_master[[f]]
  res_chain <- res[[3]][[3]]
  tmp_mat <-  sapply(1:6, function(j) quantile(res_chain[, j], c(0.025, 0.5, 0.975))) 
  tmp_list[[f]] <- tmp_mat
  
  plotCI(tmp_mat[2,short_indicies ],ui = tmp_mat[3,short_indicies],li = tmp_mat[1, short_indicies],
         ylab="Relative error",xlab = "", col = "red",
         pch = 19, xaxt="n", main = titles[f], xlim = c(0, 4))
  text(seq(1, length(short_indicies),by = 1),  par("usr")[3] - 0.15, srt = 30, adj = 1,labels = labels_vec[short_indicies], xpd = TRUE)
  abline(h=0, col = "gray",lty = "dashed")
  
  plotCI(c(1, 2), tmp_mat[2,long_indicies ],ui = tmp_mat[3,long_indicies],li = tmp_mat[1, long_indicies],
         ylab="Relative error",xlab = "", col = "red",
         pch = 19, xaxt="n", main = titles[f], xlim = c(0,3))
  text(seq(1, length(long_indicies),by = 1),  par("usr")[3] - 0.2, srt = 30, adj = 1,labels = labels_vec[long_indicies], xpd = TRUE)
  abline(h=0, col = "gray",lty = "dashed")
  
}
dev.off()

# Figure 5
filename_variable <- "AR_change_L"

titre_dat <- read.csv(paste("data/study_design/",filename_variable,"dat.csv",sep = ""))
ages <- read.csv(paste("data/study_design/",filename_variable,"ages.csv",sep = ""))

n_indiv <- dim(ages)[1]

inf_chain <- data.table::fread(paste("chains/",filename_variable,"_infection_histories.csv",sep=""))
inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]

DOBs <- ages$DOB
age_mask <- create_age_mask(DOBs, strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat,strain_isolation_times)
masks <- data.frame(cbind(age_mask, strain_mask))
n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x) nrow(masks[masks$age_mask <=x & masks$strain_mask >= x,]))  

p <- plot_attack_rates(infection_histories = inf_chain, merge(titre_dat, ages), strain_isolation_times, by_val = 10)
pg <- ggplot_build(p)
AR_variable <- pg$data[[1]]
 
viruses <- unique(titre_dat$samples)
#find based on threshold
prev_variable <- sapply(1:length(viruses), function(x){
  dat <- titre_dat[which(titre_dat$sample==viruses[x]&titre_dat$virus==viruses[x]),]
  infected <- as.numeric(dat$titre >= 3)
  tab.inf <- table(infected)
  tab.inf[2] / sum(tab.inf)
  return(tab.inf[2] / sum(tab.inf))})


####prev from constant 
#all age scenario 4
filename_constant <-  paste0( 3, "study", 3, "_L", sep = "")
titre_dat <- read.csv(paste("data/study_design/",filename_constant,"dat.csv",sep = ""))
ages <- read.csv(paste("data/study_design/",filename_constant,"ages.csv",sep = ""))

n_indiv <- dim(ages)[1]

inf_chain <- data.table::fread(paste("chains/",filename_constant,"_infection_histories.csv",sep=""))
inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"] + mcmc_pars["burnin"]),]

DOBs <- ages$DOB
age_mask <- create_age_mask(DOBs, strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat,strain_isolation_times)
masks <- data.frame(cbind(age_mask, strain_mask))
n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x) nrow(masks[masks$age_mask <=x & masks$strain_mask >= x,]))  

p <- plot_attack_rates(infection_histories = inf_chain, merge(titre_dat, ages), strain_isolation_times, by_val = 10)
pg <- ggplot_build(p)
AR_constant <- pg$data[[1]]

viruses <- unique(titre_dat$samples)
#find based on threshold
prev_constant <- sapply(1:length(viruses), function(x){
  dat <- titre_dat[which(titre_dat$sample==viruses[x]&titre_dat$virus==viruses[x]),]
  infected <- as.numeric(dat$titre >= 3)
  tab.inf <- table(infected)
  tab.inf[2] / sum(tab.inf)
  return(tab.inf[2] / sum(tab.inf))})

####plot together 
png('AR_compare.png',units='px',res=300,width=4000,height=2000)
par(mfrow=c(1,2), cex.axis=1.5, cex.main=1.5,cex.lab=1.5)
plotCI(AR_constant$x, AR_constant$y, ui=AR_constant$ymax ,li=AR_constant$ymin, 
       ylab = "Attack rate", xlab = "Virus year" , ylim = c(0, 0.6),
       main="Constant attack rate", xaxt = "n", pch = 19, xlim = c(2012, 2015), col = "black")
points(c(2014, 2015), prev_constant, col = "red", pch = 15)
lines(seq(2012, 2015, by =1 ), rep(0.15, 4), pch = 17, col = "gray", lty = "dashed")
axis(1, at = seq(2012 ,2015, by = 1),labels = c("t-3", "t-2","t-1", "t"))

plotCI(AR_variable$x, AR_variable$y, ui=AR_variable$ymax ,li=AR_variable$ymin, 
       ylab = "Attack rate", xlab = "Virus year" , ylim = c(0, 0.6),
       main="Variable attack rate", xaxt = "n", pch = 19, xlim = c(2012, 2015), col = "black")
points(c(2014, 2015), prev_variable, col = "red", pch = 15)
attack_rates_variable <- rep(c(0,0.15),length(strain_isolation_times) / 2)
lines(seq(2012, 2015, by =1 ), tail(attack_rates_variable, n = 4), pch = 17, col = "gray", lty = "dashed")
axis(1, at = seq(2012 ,2015, by = 1),labels = c("t-3", "t-2","t-1", "t"))
legend('topright',c('True attack rate','Model attack rate','Titre based attack rate'),
       pch=c(NA,19,15),col=c('gray',1,2),bty='n',
       lty = c("dashed", NA, NA))
dev.off()

## Supp figures
# Other antibody parameters

png("FigureS1.png", width = 3000, height = 1600, res = 300, units = "px")
m <- matrix(c(1, 1, 2,3,4,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.1,0.4,0.4))
par(mar = c(0,1,1,1), cex.axis=1.5, cex.main=1.5,cex.lab=1.5)
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", c("1 strain/year", "2 strain/year","4 strain/year","Cross-sectional","Longtitudinal"),
       col = c(1, "orange", "red", "black" ,"black"), lty = c(1, 1, 1, NA, NA), 
       pch = c(NA, NA, NA, 19, 17), bty = "n", ncol = 1, horiz = TRUE , cex = 1.5)

par(mar = c(3, 4, 4, 2) + 0.1)

titles <- c( expression(paste('short boost (',mu[2],'), ages 2 - 4')),
             expression(paste('short cross (',sigma[2],'), ages 2 - 4')),
             expression(paste('short boost (',mu[2],'), ages 40 - 75')),
             expression(paste('short cross (',sigma[2],'), ages 40 - 75')))
k <- 1
for(f in 1:2){
  res <- res_master[[f]]
  for(which_par in c("mu_short", "sigma2")){
    which_col <- which(which_pars == which_par)
    tmp_mat<- matrix(ncol = 3, nrow = 6)
    for(i in 1:6){
      res_chain <- res[[i]][[3]]
      tmp_mat[i, ] <- quantile(res_chain[,which_col], c(0.025,0.5,0.975))
    }
    
    plotCI(tmp_mat[, 2],ui = tmp_mat[, 3],li = tmp_mat[, 1],
           ylab="Relative error",xlab = "", xaxt = "n", col = col_vec,
           pch = pch_vec, main = titles[k])
    abline(h=0, col = "gray", lty = "dashed")
    k <- k +1
  }
}
dev.off()

png("FigureS2.png", width = 3000, height = 1600, res = 300, units = "px")
m <- matrix(c(1, 1, 2,3,4,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.1,0.4,0.4))
par(mar = c(0,1,1,1), cex.axis=1.5, cex.main=1.5,cex.lab=1.5)
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", c("1 strain/year", "2 strain/year","4 strain/year","Cross-sectional","Longtitudinal"),
       col = c(1, "orange", "red", "black" ,"black"), lty = c(1, 1, 1, NA, NA), 
       pch = c(NA, NA, NA, 19, 17), bty = "n", ncol = 1, horiz = TRUE , cex = 1.5)

par(mar = c(3, 4, 4, 2) + 0.1)

titles <- c( expression(paste('long boost (',mu[1],'), ages 2 - 4')),
             expression(paste('error (',epsilon,'), ages 2 - 4')),
             expression(paste('long boost (',mu[1],'), ages 40 - 75')),
             expression(paste('error (',epsilon,'), ages 40 - 75')))
k <- 1
for(f in 1:2){
  res <- res_master[[f]]
  for(which_par in c("mu", "error")){
    which_col <- which(which_pars == which_par)
    tmp_mat<- matrix(ncol = 3, nrow = 6)
    for(i in 1:6){
      res_chain <- res[[i]][[3]]
      tmp_mat[i, ] <- quantile(res_chain[,which_col], c(0.025,0.5,0.975))
    }
    
    plotCI(tmp_mat[, 2],ui = tmp_mat[, 3],li = tmp_mat[, 1],
           ylab="Relative error",xlab = "", xaxt = "n", col = col_vec,
           pch = pch_vec, main = titles[k])
    abline(h=0, col = "gray", lty = "dashed")
    abline(v = 3.5)
    
    k <- k +1
  }
}
dev.off()


