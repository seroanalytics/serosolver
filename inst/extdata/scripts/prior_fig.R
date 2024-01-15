library(ggpubr)
library(viridis)
library(cowplot)

strain_isolation_times <- 1968:2009
n_indiv <- 100

## Prior V2
## Plot 1A
data_prior_2 <- overall_data[overall_data$variable %in% c("total","i.1","j.1") & 
                               overall_data$metric != "Total across nm" &
                               overall_data$model == "Beta prior on times",]
theory_prior_2 <- theory[theory$metric != "Total across nm" &
                           theory$model == "Beta prior on times",]
tmp_prior_2 <- tmp[tmp$model == "Beta prior on times" & tmp$metric == "Total across m",]
p1a <- ggplot(data_prior_2[data_prior_2$metric == "Total across m",]) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,fill="grey50") +
  geom_line(data=theory_prior_2[theory_prior_2$metric == "Total across m",], aes(x=x,y=y),col="red") +
  #geom_blank(data=tmp_prior_2,aes(xmax=upper_bound/n_indiv,xmin=lower_bound/n_indiv)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Probability mass") +
  xlab(expression("Number infected in time period"~italic(j)))+
  theme_pubr()+
  theme(text=element_text(size=6))


## Plot 1B
tmp_prior_2 <- tmp[tmp$model == "Beta prior on times" & tmp$metric == "Total across n",]
p1b <- ggplot(data_prior_2[data_prior_2$metric == "Total across n",]) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,fill="grey50") +
  geom_line(data=theory_prior_2[theory_prior_2$metric == "Total across n",], aes(x=x,y=y),col="red") +
  geom_blank(data=tmp_prior_2,aes(xmax=upper_bound,xmin=lower_bound)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Probability mass") +
  xlab("Number of lifetime infections")+
  theme_pubr()+
  theme(text=element_text(size=6))

all_chain_gibbs <- load_mcmc_chains(location="~/Drive/Influenza/serosolver/infection_history_prior_methods_results/chain_no_lik/gibbs",
                                    thin=1,burnin=10000,
                                    par_tab=par_tab,
                                    unfixed=TRUE,convert_mcmc=TRUE)
inf_chain_gibbs <- all_chain_gibbs$inf_chain
## Generate lower, upper and median cumulative infection histories from the
## MCMC chain
tmp_inf_chain <- inf_chain_gibbs[inf_chain_gibbs$i ==1, ]
hist_profiles <- ddply(tmp_inf_chain, .(i, samp_no, chain_no), function(x) {
  empty <- numeric(length(strain_isolation_times))
  empty[x[x$x == 1, "j"]] <- 1
  cumsum(empty)
})


hist_profiles <- hist_profiles[, colnames(hist_profiles) != "samp_no"]
colnames(hist_profiles) <- c("i", "chain_no", strain_isolation_times)
hist_profiles_tmp <- melt(hist_profiles, id.vars=c("i","chain_no"))
res <- ddply(hist_profiles_tmp, .(i, chain_no, variable), function(x) quantile(x$value,c(0.025,0.25,0.5,0.75,0.975)))
colnames(res) <- c("individual", "chain_no", "variable", "lower","lower_mid","median","upper_mid","upper")
res$variable <- as.numeric(res$variable) - 1
p1c <- ggplot(res) +
  geom_ribbon(aes(x=variable, ymin=lower,ymax=upper),fill="gray90") +
  geom_ribbon(aes(x=variable, ymin=lower_mid,ymax=upper_mid),fill="gray60") +
  geom_line(aes(x=variable,y=median)) +
  coord_cartesian(ylim=c(0,n_times),xlim=c(0,n_times)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Number of time periods alive") +
  ylab("Cumulative infections") +
  theme_pubr()+
  theme(text=element_text(size=6))

title1 <- ggdraw() + draw_label("Prior 1/2, uniform/beta prior on attack rates (α=β=1)",fontface="bold",x=0.01,hjust=0,size=8)
row1 <- plot_grid(title1, plot_grid(p1a, p1b, p1c,nrow=1,labels="AUTO",label_size=8),nrow=2,rel_heights = c(0.1,1))


## Prior V3
## Plot 2A
data_prior_3 <- overall_data[overall_data$variable %in% c("total","i.1","j.1") & 
                               overall_data$metric != "Total across nm" &
                               overall_data$model == "Beta prior on individuals",]
theory_prior_3 <- theory[theory$metric != "Total across nm" &
                           theory$model == "Beta prior on individuals",]
tmp_prior_3 <- tmp[tmp$model == "Beta prior on individuals" & tmp$metric == "Total across m",]
p2a <- ggplot(data_prior_3[data_prior_3$metric == "Total across m",]) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,fill="grey50") +
  geom_line(data=theory_prior_3[theory_prior_3$metric == "Total across m",], aes(x=x,y=y),col="red") +
  #geom_blank(data=tmp_prior_2,aes(xmax=upper_bound/n_indiv,xmin=lower_bound/n_indiv)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Probability mass") +
  xlab(expression("Number infected in time period"~italic(j)))+
  theme_pubr() +
  theme(text=element_text(size=6))


## Plot 2B
tmp_prior_3 <- tmp[tmp$model == "Beta prior on individuals" & tmp$metric == "Total across n",]
p2b <- ggplot(data_prior_3[data_prior_3$metric == "Total across n",]) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,fill="grey50") +
  geom_line(data=theory_prior_3[theory_prior_3$metric == "Total across n",], aes(x=x,y=y),col="red") +
  geom_blank(data=tmp_prior_3,aes(xmax=upper_bound,xmin=lower_bound)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Probability mass") +
  xlab("Number of lifetime infections")+
  theme_pubr()+
  theme(text=element_text(size=6))

all_chain_bb <- load_mcmc_chains(location="~/Drive/Influenza/serosolver/infection_history_prior_methods_results/chain_no_lik/beta_binomial/",
                                    thin=1,burnin=10000,
                                    par_tab=par_tab,
                                    unfixed=TRUE,convert_mcmc=TRUE)
all_chain_bb <- all_chain_bb$inf_chain
## Generate lower, upper and median cumulative infection histories from the
## MCMC chain
tmp_inf_chain <- all_chain_bb[all_chain_bb$i ==1, ]
hist_profiles <- ddply(tmp_inf_chain, .(i, samp_no, chain_no), function(x) {
  empty <- numeric(length(strain_isolation_times))
  empty[x[x$x == 1, "j"]] <- 1
  cumsum(empty)
})

hist_profiles <- hist_profiles[, colnames(hist_profiles) != "samp_no"]
colnames(hist_profiles) <- c("i", "chain_no", strain_isolation_times)
hist_profiles_tmp <- melt(hist_profiles, id.vars=c("i","chain_no"))
res <- ddply(hist_profiles_tmp, .(i, chain_no, variable), function(x) quantile(x$value,c(0.025,0.25,0.5,0.75,0.975)))
colnames(res) <- c("individual", "chain_no", "variable", "lower","lower_mid","median","upper_mid","upper")
res$variable <- as.numeric(res$variable) - 1
p2c <- ggplot(res) +
  geom_ribbon(aes(x=variable, ymin=lower,ymax=upper),fill="gray90") +
  geom_ribbon(aes(x=variable, ymin=lower_mid,ymax=upper_mid),fill="gray60") +
  geom_line(aes(x=variable,y=median)) +
  coord_cartesian(ylim=c(0,n_times),xlim=c(0,n_times)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Number of time periods alive") +
  ylab("Cumulative infections") +
  theme_pubr()+
  theme(text=element_text(size=6))

title1 <- ggdraw() + draw_label("Prior 3, beta-binomial prior on number of lifetime infections (α=β=1)",
                                fontface="bold",x=0.01,hjust=0,size=8)
row2 <- plot_grid(title1, plot_grid(p2a, p2b, p2c,nrow=1,labels=c("D","E","F"),label_size=8),nrow=2,rel_heights = c(0.1,1))


## Prior V4
## Plot 3A
data_prior_4 <- overall_data[overall_data$variable %in% c("total","i.1","j.1") & 
                               overall_data$metric != "Total across nm" &
                               overall_data$model == "Beta prior on total",]
theory_prior_4 <- theory[theory$metric != "Total across nm" &
                           theory$model == "Beta prior on total",]
tmp_prior_4 <- tmp[tmp$model == "Beta prior on total" & tmp$metric == "Total across m",]
p3a <- ggplot(data_prior_4[data_prior_4$metric == "Total across m",]) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,fill="grey50") +
  geom_line(data=theory_prior_4[theory_prior_4$metric == "Total across m",], aes(x=x,y=y),col="red") +
  #geom_blank(data=tmp_prior_2,aes(xmax=upper_bound/n_indiv,xmin=lower_bound/n_indiv)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Probability mass") +
  xlab(expression("Number infected in time period"~italic(j)))+
  theme_pubr()+
  theme(text=element_text(size=6))


## Plot 3B
tmp_prior_4 <- tmp[tmp$model == "Beta prior on total" & tmp$metric == "Total across n",]
p3b <- ggplot(data_prior_4[data_prior_4$metric == "Total across n",]) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,fill="grey50") +
  geom_line(data=theory_prior_4[theory_prior_4$metric == "Total across n",], aes(x=x,y=y),col="red") +
  geom_blank(data=tmp_prior_4,aes(xmax=upper_bound,xmin=lower_bound)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Probability mass") +
  xlab("Number of lifetime infections")+
  theme_pubr()+
  theme(text=element_text(size=6))

all_chain_overall <- load_mcmc_chains(location="~/Drive/Influenza/serosolver/infection_history_prior_methods_results/chain_no_lik/overall//",
                                 thin=1,burnin=10000,
                                 par_tab=par_tab,
                                 unfixed=TRUE,convert_mcmc=TRUE)
all_chain_overall <- all_chain_overall$inf_chain
## Generate lower, upper and median cumulative infection histories from the
## MCMC chain
tmp_inf_chain <- all_chain_overall[all_chain_overall$i ==1, ]
hist_profiles <- ddply(tmp_inf_chain, .(i, samp_no, chain_no), function(x) {
  empty <- numeric(length(strain_isolation_times))
  empty[x[x$x == 1, "j"]] <- 1
  cumsum(empty)
})

hist_profiles <- hist_profiles[, colnames(hist_profiles) != "samp_no"]
colnames(hist_profiles) <- c("i", "chain_no", strain_isolation_times)
hist_profiles_tmp <- melt(hist_profiles, id.vars=c("i","chain_no"))
res <- ddply(hist_profiles_tmp, .(i, chain_no, variable), function(x) quantile(x$value,c(0.025,0.25,0.5,0.75,0.975)))
colnames(res) <- c("individual", "chain_no", "variable", "lower","lower_mid","median","upper_mid","upper")
res$variable <- as.numeric(res$variable) - 1
p3c <- ggplot(res) +
  geom_ribbon(aes(x=variable, ymin=lower,ymax=upper),fill="gray90") +
  geom_ribbon(aes(x=variable,ymin=lower_mid,ymax=upper_mid),fill="grey60") +
  geom_line(aes(x=variable,y=median)) +
  coord_cartesian(ylim=c(0,n_times),xlim=c(0,n_times)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Number of time periods alive") +
  ylab("Cumulative infections") +
  theme_pubr()+
  theme(text=element_text(size=6))

title1 <- ggdraw() + draw_label("Prior 4, beta on any infection (α=β=1)",
                                fontface="bold",x=0.01,hjust=0,size=8)
row3 <- plot_grid(title1, plot_grid(p3a, p3b, p3c,nrow=1,labels=c("G","H","I"),label_size=8),nrow=2,rel_heights = c(0.1,1))


overall_p <- plot_grid(row1,row2,row3,nrow=3)
print(overall_p)

cairo_pdf("~/Drive/Influenza/serosolver/methods_paper/PLOS Comp Biol/Figures/Fig3.pdf",width=5.2, height=4.5)
overall_p
dev.off()
