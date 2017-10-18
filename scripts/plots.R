#devtools::load_all()

#infectionHistories1 <- data.table::fread("test1_infectionHistories.csv")
#infectionHistories <- infectionHistories1[infectionHistories1$sampno %in% 
    #                                        unique(infectionHistories1$sampno)[seq(1,length(unique(infectionHistories1$sampno)),by=100)],]
#infectionHistories <- infectionHistories[infectionHistories$sampno == chain[which.max(chain$lnlike),"sampno"]-1,1:47]
i <- 100

dat <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)

tmpDat <- dat[dat$individual==i,]
unique(tmpDat$samples)
tmpDat <- tmpDat[tmpDat$samples == unique(tmpDat$samples)[2],]
tmpDat <- tmpDat[complete.cases(tmpDat),]
tmp <- infectionHistories1[infectionHistories1$individual == i,1:47]
infCounts <- apply(tmp[,1:47],2,function(x) table(factor(x, levels=c(0,1))))
infSD <- apply(tmp[,1:47],2,sd)
y <- (infCounts[2,]/colSums(infCounts))
#print(infSD)
x <- infectionHistories[i,]
wow <- data.frame(real=x,estimated=y)
omg <- data.frame(year=rownames(wow),density=wow$estimated)
omg$year <- as.numeric(as.character(omg$year))

antigenicMap <- read.csv("data/antigenic_map.csv",stringsAsFactors=FALSE)
strainIsolationTimes <- unique(antigenicMap$inf_years)
antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")]))

bestpars <- zikaProj::get_best_pars(chain)
antigenicMapLong <- 1-pars["sigma1"]*antigenicMapMelted
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort <- 1-pars["sigma2"]*antigenicMapMelted
antigenicMapShort[antigenicMapShort < 0] <- 0

infectionHistory <- infectionHistories[infectionHistories$sampno == 89500,1:47]
predicted_titres <- matrix(nrow=1000,ncol=length(strainIsolationTimes))


tmp <- as.data.frame(infectionHistories1[infectionHistories1$individual == i,])
tmpSamp <- sample(intersect(unique(tmp$sampno),unique(chain$sampno)-1),1000)

for(i in 1:1000){
  sampno <-tmpSamp[i]
  index <- which(chain$sampno == sampno+1)
  
  infectionHistory <- tmp[tmp$sampno == sampno,1:47]
  pars <- zikaProj::get_index_pars(chain,index)
  predicted_titres[i,] <- infection_model_indiv(pars, as.numeric(infectionHistory), unique(tmpDat$samples)[1],
                                       strainIsolationTimes,antigenicMapLong,antigenicMapShort)
}
dat2 <- t(apply(predicted_titres,2, function(x) quantile(x, c(0.025,0.5,0.975))))
dat2 <- cbind(year=strainIsolationTimes,dat2, max=infection_model_indiv(bestpars, as.numeric(infectionHistory), unique(tmpDat$samples)[1],
                                  strainIsolationTimes,antigenicMapLong,antigenicMapShort))
dat2 <- as.data.frame(dat2)
colnames(dat2) <- c("year","lower","median","upper","max")

p100 <- ggplot(tmpDat) + geom_point(aes(x=virus,y=titre),col="red",size=1) + geom_vline(data=omg,aes(xintercept=year,alpha=density)) + 
  geom_line(data=dat2,aes(x=year,y=max), col="blue") +
  geom_ribbon(data=dat2,aes(ymin=lower,ymax=upper,x=year),alpha=0.5,fill="blue") +
  scale_y_continuous(limits=c(0,8)) +
   xlab("Year") +
  ylab("Recorded titre") +
    theme_bw() + 
  theme(panel.grid=element_blank())+
  scale_alpha(limits=c(0,1))


png("test.png",height=7,width=8,unit="in",res=300)
plot_grid(p100,p250,p350,ncol=1)
dev.off()
