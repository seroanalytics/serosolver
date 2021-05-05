library(ggplot2)
#############
## Looking at antigenic coordinates of strains that have circulated
## over the years
#############
setwd("~/Drive/Influenza/serosolver/Measurement_error")
devtools::load_all("~/Documents/Fluscape/serosolver")
#devtools::test("~/Documents/Fluscape/serosolver",show_report=TRUE)


extract_and_clean <- function(coords, colI=2,sep_char='[/]'){
  viruses <- sapply(coords[,colI], function(x) strsplit(x,sep_char))
  years <- as.numeric(unlist(lapply(viruses, function(x) x[length(x)])))
  years[years < 10] <- paste0(200,years[years<10])
  years <- as.numeric(years)
  years[years < 50] <- paste0(20,years[years<50])
  years <- as.numeric(years)
  years[years < 100 & years > 50] <- paste0(19,years[years < 100 & years > 50])
  years <- as.numeric(years)
  return(years)
}

###########
## Extract coordinates form Fonville and clean
###########
fonville_coords <- read.csv("data/fonville_coords.csv",stringsAsFactors=FALSE)
years_fonville <- extract_and_clean(fonville_coords,1)
fonville_coords <- data.frame(set="fonville",year=years_fonville, x=fonville_coords$AG.coordinate.1,y=fonville_coords$AG.coordinate.2)
#p1 <- ggplot(fonville_coords) + geom_point(aes(x=y,y=x,col=as.factor(year)))

##########
## Extract coords from Smith and clean
##########
smith_coords <- read.table("smith2004/smith_coords.txt",stringsAsFactors=FALSE)
years_smith <- extract_and_clean(smith_coords,2)
smith_coords <- data.frame(set="smith",year=years_smith,x=smith_coords[,3],y=smith_coords[,4])
smith_coords$y <- -smith_coords$y
#p2 <- ggplot(smith_coords) + geom_point(aes(x=y,y=x,col=as.factor(year)))

coords <- rbind(fonville_coords, smith_coords)

clusters <- read.csv("clusters.csv")

coords <- merge(clusters, coords)
coords1 <- coords[coords$set == "smith",]
coords1 <- coords1[,c("x","y","year")]
colnames(coords1) <- c("Y","X","Strain")

averages <- plyr::ddply(coords1, ~Strain, function(x) c(sum(x$Y)/nrow(x), sum(x$X)/nrow(x) ))
colnames(averages) <- c("Strain","Y","X")
fit_dat <- generate_antigenic_map_flexible(coords1, spar=0.7,year_min=1968,year_max=2010)

plot(fit_dat[,c("x_coord","y_coord")],type="l")

titre_dat <- read.csv("~/Documents/Fluscape/serosolver/tests/testdata/fluscape_sim_annual_dat.csv")
strain_isolation_times <- fit_dat$inf_years
infection_history_mat <- setup_infection_histories_new_2(titre_dat, strain_isolation_times, 5, 2, sample_prob=0.1)


fake_data <- expand.grid(individual=1:25,samples=1968:2010,virus=strain_isolation_times,titre=0,run=1,DOB=1900)
fake_data <- fake_data[order(fake_data$individual, fake_data$run, fake_data$samples, fake_data$virus), ]
par_tab <- read.csv("~/Documents/Fluscape/serosolver/tests/testdata/par_tab_base.csv",stringsAsFactors=FALSE)
par_tab[par_tab$names == "wane","values"] <- 0.5
par_tab[par_tab$names == "tau","values"] <- 0.1
par_tab <- par_tab[par_tab$names != "lambda",]
f <- create_posterior_func_fast(par_tab, fake_data,fit_dat,version=2,solve_likelihood=TRUE,function_type = 3)
inf_hist <- matrix(0, nrow=5,ncol=length(strain_isolation_times))
inf_hist[1,] <- 0
inf_indices <- c(5,12)
inf_hist[1, inf_indices] <- 1
y <- f(par_tab$values, inf_hist)
fake_data$y <- f(par_tab$values, inf_hist)
ggplot(fake_data[fake_data$individual == 1 & fake_data$virus %in% strain_isolation_times[c(5,12,20)] & fake_data$samples < 1985,]) + 
  geom_line(aes(x=samples,y=2*(y^5),col=as.factor(virus)))




p1 <- ggplot(coords[coords$set=="smith",]) + 
  geom_line(data=fit_dat,aes(x=x_coord,y=y_coord),linetype="dashed",size=1,col="gray20")+
  geom_point(aes(x=y,y=x,col=as.factor(year)),size=2) + 
  #geom_point(data=averages,aes(x=X,y=Y,fill=as.factor(Strain)),size=4, shape=21) +
  scale_y_continuous(expand=c(0,0),limits=c(-8,8),breaks=seq(-8,8,by=1), labels=NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(-17,16),breaks=seq(-18,18,by=1), labels=NULL)+
  theme_bw() +
  xlab("Antigenic distance")+
  ylab("Antigenic distance")+
  theme(legend.position="none", panel.grid.major=element_blank(), axis.ticks = element_blank())
p1
  #scale_color_gradient2(low="blue",high="red",mid="green",midpoint=1990)#+
 # facet_wrap(~set,ncol=1,scales="free")

svg("~/Dropbox/serosolver/OrthomyxoConf/antigenic_map.svg",width=4,height=4)
plot(map)
dev.off()



averages <- ddply(smith_coords, ~year, function(x) c(sum(x$x)/nrow(x), sum(x$y)/nrow(x) ))
library(ggrepel)
map <- ggplot(averages) + 
  geom_path(data=averages,aes(x=V2,y=V1),col="gray",alpha=0.5) +
  geom_point(data=averages,aes(x=V2,y=V1, col=as.factor(year)),size=3) + 
  geom_label_repel(data=averages[averages$year %in% seq(1968,2002,by=2),],aes(x=V2,y=V1,label=year),label.size=NA,size=4,fill=NA,family="Arial") +
    #scale_color_brewer(palette="paired") +
  theme_bw()+
  theme(legend.position = "none",
        text=element_text(family="Arial",size=12)) +
  xlab("Antigenic dimension 1") +
  ylab("Antigenic dimension 2")
map


svg("~/Dropbox/serosolver/OrthomyxoConf/antigenic_map.svg",width=4,height=4)
plot(map)
dev.off()

