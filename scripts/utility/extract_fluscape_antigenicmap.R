library(ggplot2)

## Antigenic map fit
buckets <- 1
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fluscape_virus_key <- c("BJ89"=1989, "SC87"=1987, "PH82"=1982, "BR07"=2007, "WU95"=1995, "BK79"=1979, "HK14"=2014, "TX12"=2012, 
                        "PE09"=2010, "BJ92"=1992, "TX77"=1977, "VC09"=2009, "CL04"=2004, "VC98"=1998, "FJ00"=2000, "VC75"=1975, 
                        "MS85"=1985, "FJ02"=2002, "EN72"=1972, "X31"=1969, "HK68"=1968)*buckets

generate_antigenic_map <- function(antigenicDistances, buckets=1){
  ## Following assumptions:
  ## 1. X31 == 1969
  ## 2. PE2009 is like the strain circulating in 2010
  virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
                 "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets

  antigenicDistances$Strain <- virus_key[antigenicDistances$Strain]
  fit <- smooth.spline(antigenicDistances$X,antigenicDistances$Y,spar=0.3)
  x_line <- lm(data = antigenicDistances, X~Strain)
  Strain <- seq(1968*buckets,2015*buckets,by=1)
  x_predict <- predict(x_line,data.frame(Strain))
  y_predict <- predict(fit, x=x_predict)
  fit_dat <- data.frame(x=y_predict$x,y=y_predict$y)
  fit_dat$strain <- Strain
  colnames(fit_dat) <- c("x_coord","y_coord","inf_years")
  return(fit_dat)
}
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets
antigenicMap$Strain <- virus_key[antigenicMap$Strain]
p1 <- ggplot(antigenicMap) + 
  geom_line(data=fit_dat,aes(x=x_coord,y=y_coord), col="red") +
  geom_point(data=antigenicMap,aes(x=X,y=Y)) + 
  geom_label(data=antigenicMap,aes(x=X+4,y=Y+0.25,label=Strain)) +
  theme_bw()
euc_distance <- function(i1, i2, fit_dat){
  return(sqrt((fit_dat[i1,"x_coord"] - fit_dat[i2,"x_coord"])^2 + (fit_dat[i1,"y_coord"] - fit_dat[i2,"y_coord"])^2))
}

#write.table(fit_dat,"data/fluscape_map_quarter.csv",row.names=FALSE,sep=",")
