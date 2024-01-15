inpputDat2 <- vacDat_h1n1
########################
## Plot AR by age group in unvac
########################
### Plot inferred attack rates
# arguments
year_range <- min(strain_isolation_times):max(strain_isolation_times)
# plot_attack_rates_monthly function altered to add additional line for simulated attack rates (to compare against inferred attach rates)
months <- 1:buckets
years <- range(floor(year_range/buckets))
years <- years[1]:years[2]
labels <- c(sapply(years, function(x) paste0(months, "/",x)))
labels1 <- labels[1:length(year_range)]
labels1 <- labels1[seq(1,length(labels1),by=buckets)]
year_break <- year_range[seq(1,length(year_range),by=buckets)]

n_indiv <- length(unique(titre_dat2$individual))
# determine n_alive for each age group
n_alive <- c(rep(n_indiv,length(year_range)))

age_ids <- unique(inputDat2[,c("individual","Age_1")])

DOBs1 <- unique(titre_dat2[,c("individual","DOB")])[,2]
ageMask <- create_age_mask(DOBs1, strain_isolation_times)
# Create strain mask
strainMask <- create_strain_mask(titre_dat2,strain_isolation_times)
masks <- data.frame(cbind(ageMask, strainMask))
# Number of people that were born before each year and have had a sample taken since that year happened
n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x) nrow(masks[masks$ageMask <=x & masks$strainMask >= x,]))    

titre_dat_grouped <- merge(titre_dat2, age_ids)
titre_dat_grouped$group <- as.numeric(titre_dat_grouped$Age_1)

# inferred attack rates
inf_chain <- all_chains_vac$inf_chain
inf_chain <- pad_inf_chain(inf_chain)

inf_chain1 <- inf_chain[inf_chain$i%in%which(age_ids$Age_1=='<19'),]
data.table::setkey(inf_chain1, "j","samp_no","chain_no")
tmp1 <- inf_chain1[,list(V1=sum(x)),by=key(inf_chain1)]
quantiles1 <- ddply(tmp1, ~j, function(x) quantile(x$V1, c(0.025,0.1,0.5,0.9,0.975)))
quantiles1$age_group <- '<19'

inf_chain2 <- inf_chain[inf_chain$i%in%which(age_ids$Age_1=='19-64'),]
data.table::setkey(inf_chain2, "j","samp_no","chain_no")
tmp2 <- inf_chain2[,list(V1=sum(x)),by=key(inf_chain2)]
quantiles2 <- ddply(tmp2, ~j, function(x) quantile(x$V1, c(0.025,0.1,0.5,0.9,0.975)))
quantiles2$age_group <- '19-64'

inf_chain3 <- inf_chain[inf_chain$i%in%which(age_ids$Age_1=='>64'),]
data.table::setkey(inf_chain3, "j","samp_no","chain_no")
tmp3 <- inf_chain3[,list(V1=sum(x)),by=key(inf_chain3)]
quantiles3 <- ddply(tmp3, ~j, function(x) quantile(x$V1, c(0.025,0.1,0.5,0.9,0.975)))
quantiles3$age_group <- '>64'

quantiles <- rbind(quantiles1,quantiles2,quantiles3)
colnames(quantiles) <- c("j","lower","lower_50","median","upper_50","upper","age_group")
#quantiles[c("lower","median","upper")] <- quantiles[c("lower","median","upper")]/n_alive[1]
quantiles$year <- year_range[quantiles$j]
quantiles$age_group <- factor(quantiles$age_group,levels=c("<19","19-64",">64"))

n_alive_group <- get_n_alive_group(titre_dat_grouped, strain_isolation_times,melt_dat = TRUE)
n_alive_group$age_group <- levels(age_ids$Age_1)[n_alive_group$group]
n_alive_group$age_group <- factor(n_alive_group$age_group)
quantiles <- merge(quantiles, n_alive_group,by=c("j","age_group"))
quantiles[c("lower","lower_50","median","upper_50","upper")] <- quantiles[c("lower","lower_50","median","upper_50","upper")]/quantiles$n_alive

# Labels
i <- 1:(length(strain_isolation_times) + 1)
labels <- as.Date("01/01/2009", format="%d/%m/%Y") + ((i-1)*365/4 + 365/8)
labels1 <- as.Date("01/01/2009", format="%d/%m/%Y") + ((i-2)*365/4)
labels2 <- c("2009-Q1","2009-Q2","2009-Q3","2009-Q4",
             "2010-Q1","2010-Q2","2010-Q3","2010-Q4",
             "2011-Q1","2011-Q2","2011-Q3","2011-Q4")

### plot
colnames(quantiles)[which(colnames(quantiles) == "age_group")] <- "Age group:"
ymax <- 1
p_age <- ggplot(quantiles) + 
  geom_ribbon(aes(x=year, ymin=lower,ymax=upper,fill=`Age group:`),alpha=0.25) +   
  geom_ribbon(aes(x=year, ymin=lower_50,ymax=upper_50,fill=`Age group:`),alpha=0.5) +   
  geom_line(aes(x=year,y=median,colour=`Age group:`),size=0.75)+ 
  geom_point(aes(x=year,y=median,colour=`Age group:`),size=0.75)+ 
  scale_y_continuous(expand=c(0,0),limits=c(-0.005,ymax),breaks=seq(0,ymax,by=0.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = strain_isolation_times1, labels = labels2,
                     limits=c(min(strain_isolation_times-0.1),max(strain_isolation_times+0.1))) +
  scale_fill_manual(values=colour_fills_age) +
  scale_color_manual(values=colour_fills_age) +
  theme_pubr() +
  theme()+
  ylab("Estimated per capita\n incidence (per quarter)") +
  xlab("Time of infection")+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.key=element_rect(color=NA),
        legend.background = element_blank(),
        legend.title = element_text(size=8,family="sans"),
        legend.text=element_text(size=8,family="sans"),
        legend.position = c(0.8, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_blank(),
        legend.direction = "horizontal",
        legend.margin = margin(6, 6, 6, 6),
        axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=margin(l=10,r=5,t=5)
  )
p_age
