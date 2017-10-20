library(tidyr)
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

firstYear <- 1968
lastYear <- 2015
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, 
               "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)

fluscape_virus_key <- c("BJ89"=1989, "SC87"=1987, "PH82"=1982, "BR07"=2007, "WU95"=1995, "BK79"=1979, "HK14"=2014, "TX12"=2012, 
                        "PE09"=2009, "BJ92"=1992, "TX77"=1977, "VC09"=2009, "CL04"=2004, "VC98"=1998, "FJ00"=2000, "VC75"=1975, 
                        "MS85"=1985, "FJ02"=2002, "EN72"=1972, "X31"=1931, "HK68"=1968)


## Extract fluscape data
titreDat <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/HI_titers_paired_R56Pilot.csv",stringsAsFactors = FALSE)
titreDat <- titreDat[,c("Visit","Virus","HI_Titer","Participant_ID")]
titreDat[titreDat$Visit == "v1","Visit"] <- "V1"
titreDat[titreDat$Visit == "20","Visit"] <- "V2"
viruses <- unique(titreDat$Virus)
colnames(titreDat) <- c("visit","virus","titre","individual")
v1 <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Participants_V1.csv",stringsAsFactors = FALSE)
v2 <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Participants_V2.csv",stringsAsFactors = FALSE)
v3 <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Participants_V3.csv",stringsAsFactors = FALSE)
v4 <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Participants_V4.csv",stringsAsFactors = FALSE)

needed_names <- c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_SAMPLE_TIME")
v1 <- v1[,needed_names[!(needed_names %in% c("PART_BIRTH_YEAR","PART_BIRTH_MONTH"))]]
v1$PART_SAMPLE_TIME <- as.Date(v1$PART_SAMPLE_TIME)
colnames(v1)[colnames(v1) == "PART_SAMPLE_TIME"] <- "V1"

v2 <- v2[,needed_names]
v2$PART_SAMPLE_TIME <- as.Date(v2$PART_SAMPLE_TIME)
colnames(v2)[colnames(v2) == "PART_SAMPLE_TIME"] <- "V2"

v3 <- v3[,needed_names]
v3$PART_SAMPLE_TIME <- as.Date(v3$PART_SAMPLE_TIME)
colnames(v3)[colnames(v3) == "PART_SAMPLE_TIME"] <- "V3"

v4 <- v4[,c(needed_names,"PART_BIRTH_YEAR","PART_BIRTH_MONTH")]
v4$PART_SAMPLE_TIME <- as.Date(v4$PART_SAMPLE_TIME)
colnames(v4)[colnames(v4) == "PART_SAMPLE_TIME"] <- "V4"

part_info <- merge(v1, v2, id.vars=c("PARTICIPANT_ID","HH_ID","LOC_ID"),all=TRUE)
part_info <- merge(part_info,v3, id.vars=c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_BIRTH_YEAR","PART_BIRTH_MONTH"),all=TRUE)
part_info <- merge(part_info,v4, id.vars=c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_BIRTH_YEAR","PART_BIRTH_MONTH"),all=TRUE)

PARTICIPANT_ID <- sprintf("L%02dH%02dP%02d",part_info$LOC_ID,part_info$HH_ID,part_info$PARTICIPANT_ID)

maxTime <- max(part_info$V4,na.rm=TRUE)
part_info$PART_BIRTH_YEAR[part_info$PART_BIRTH_YEAR %in% c(999, 888, 0, "\\N")] <- NA
part_info$PART_BIRTH_MONTH[part_info$PART_BIRTH_MONTH %in% c(999, 888, 0, "\\N")] <- NA

birthDates <- apply(part_info, 1, function(x){
  if(!is.na(x["PART_BIRTH_MONTH"]) & !is.na(x["PART_BIRTH_YEAR"])){
   sprintf("%04d-%02d-%02d",as.numeric(x["PART_BIRTH_YEAR"]),as.numeric(x["PART_BIRTH_MONTH"]),1)
  } else {
    NA
}})

birthDates <- as.Date(birthDates) 
part_info$age <- as.numeric(round((maxTime - birthDates)/365))
part_info$DOB <- as.numeric(format(birthDates,"%Y"))
part_info$V1 <-  as.numeric(format(part_info$V1, "%Y"))
part_info$V2 <-  as.numeric(format(part_info$V2, "%Y"))
part_info$V3 <-  as.numeric(format(part_info$V3, "%Y"))
part_info$V4 <-  as.numeric(format(part_info$V4, "%Y"))
part_info$ID <- PARTICIPANT_ID
part_info <- part_info[part_info$ID %in% unique(titreDat$individual),]

##########################################
## Extract ages
ages <- part_info[,c("ID","DOB")]
ages <- ages[ages$ID %in% unique(titreDat$individual),]
colnames(ages) <- c("individual","DOB")
##########################################

samplingTimes <- unique(c(as.matrix(part_info[,c("V1","V2","V3","V4")])))
samplingTimes <- samplingTimes[!is.na(samplingTimes)]
samplingTimes <- samplingTimes[order(samplingTimes)]

tmp <- part_info[,c("V1","V2","V3","V4","ID")]
tmp <- reshape2::melt(tmp,id.vars="ID")
colnames(tmp) <- c("individual","visit","samples")
tmp <- tmp[complete.cases(tmp),]

titreDat$individual <- match(titreDat$individual, PARTICIPANT_ID)
tmp$individual <- match(tmp$individual, PARTICIPANT_ID)
ages$individual <- match(ages$individual,PARTICIPANT_ID)

finalDat <- merge(titreDat, tmp, by=c("visit","individual"))

#### NEED TO FIX THIS - REMOVING DUPLICATED DATA FOR THE SAME SAMPLE TIME
finalDat <- finalDat[!duplicated(finalDat[,c("visit","individual","virus","samples")]),]
finalDat <- finalDat[order(finalDat$indiv,finalDat$samples, finalDat$virus),]
ids <- unique(finalDat$individual)
ages$individual <- match(ages$individual, ids)
finalDat$individual <- match(finalDat$individual,ids)
finalDat <- finalDat[,c("individual","virus","samples","titre")]
finalDat$virus <- fluscape_virus_key[finalDat$virus]


#####
## Expand out all individuals and virus names
#####
all_viruses <- expand.grid("individual"=unique(finalDat$individual),"virus"=seq(firstYear,lastYear,by=1))
all_samples <- unique(finalDat[,c("individual","samples")])
all_viruses <- merge(all_viruses, all_samples)
all_viruses <- all_viruses[order(all_viruses$individual,all_viruses$samples,all_viruses$virus),]

tmp <- merge(finalDat, all_viruses, by=c("individual","virus","samples"), all=TRUE)
tmp <- tmp[order(tmp$individual,tmp$samples,tmp$individual),]
finalDat <- tmp

### NEED TO FIX THIS - NOT SURE WHICH YEAR X31 CORRESPONDS TO
finalDat <- finalDat[finalDat$virus != 1931,]

finalDat$group <- 1
finalDat[!is.na(finalDat$titre) & finalDat$titre == 0,"titre"] <- 5
finalDat[!is.na(finalDat$titre),"titre"] <- log2(finalDat[!is.na(finalDat$titre),"titre"]/5)

write.table(finalDat,"data/fluscape_data.csv",sep=",",row.names=FALSE)
write.table(ages,"data/fluscape_ages.csv",sep=",",row.names=FALSE)
