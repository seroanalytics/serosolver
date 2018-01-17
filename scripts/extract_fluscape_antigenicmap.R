library(ggplot2)

## Antigenic map fit
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
viruses <- unique(antigenicMap$Strain)
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)

## Following assumptions:
## 1. X31 == 1969
## 2. PE2009 is like the strain circulating in 2010

fluscape_virus_key <- c("BJ89"=1989, "SC87"=1987, "PH82"=1982, "BR07"=2007, "WU95"=1995, "BK79"=1979, "HK14"=2014, "TX12"=2012, 
  "PE09"=2010, "BJ92"=1992, "TX77"=1977, "VC09"=2009, "CL04"=2004, "VC98"=1998, "FJ00"=2000, "VC75"=1975, 
  "MS85"=1985, "FJ02"=2002, "EN72"=1972, "X31"=1969, "HK68"=1968)


antigenicMap$Strain <- virus_key[antigenicMap$Strain]
fit <- smooth.spline(antigenicMap$X,antigenicMap$Y,spar=0.3)
x_line <- lm(data = antigenicMap, X~Strain)
Strain <- seq(1968,2015,by=1)
x_predict <- predict(x_line,data.frame(Strain))
y_predict <- predict(fit, x=x_predict)

fit_dat <- data.frame(x=y_predict$x,y=y_predict$y)
fit_dat$strain <- Strain

p1 <- ggplot(antigenicMap) + 
    geom_line(data=fit_dat,aes(x=x,y=y), col="red") +
    geom_point(data=antigenicMap,aes(x=X,y=Y)) + 
    geom_label(data=antigenicMap,aes(x=X+4,y=Y+0.25,label=Strain)) +
    theme_bw()

colnames(fit_dat) <- c("x_coord","y_coord","inf_years")
write.table(fit_dat,"data/fluscape_map_quarter.csv",row.names=FALSE,sep=",")
