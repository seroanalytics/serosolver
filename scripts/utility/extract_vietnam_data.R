#### Extract vietnam data
options(StringsAsFactors=F)
data1=read.csv("~/Documents/Fluscape/flu-model/sero_model/datasets/HaNamCohort.csv", as.is=T)

# List test strains
nstrains=length(data1)-2 # remove subject and sample year
strain_names=names(data1)[3:(nstrains+2)]
test.index=c(1:nstrains)

# Convert to log titres and set missing data = NA
data1[data1=="*"]=NA
data1[,strain_names]=apply(data1[,strain_names],2,function(x){log2(as.numeric(x)/10)+1}) 

# Convert names into strain years
strain_years=as.numeric(sapply(strain_names,function(x){
  a1=max(which(strsplit(x, "")[[1]]=="."))
  lstr=nchar(x)
  yr1=substr(x, a1+1, lstr)
  
  if(nchar(yr1)>4){yr1=substr(yr1, 1, 4)}
  year=yr1
  if(nchar(yr1)==2 & as.numeric(yr1)>15){year=paste("19",yr1,sep="")}
  if(nchar(yr1)==2 & as.numeric(yr1)<15){year=paste("20",yr1,sep="")}
  year
}
))
#strain_years_unique=sort(unique(strain_years)) # years of samples tested against
# Gather participants and infection years
names(strain_years) <- strain_names
n_part=max(data1$Subject.number) # number of participants

test_years=unique(data1$Sample.year) # year of testing
test.n=length(test_years) # number of test years

inf_years=seq(min(strain_years),max(c(test_years,strain_years))) #annual infection model
inf.n=length(inf_years) # number of possible infecting strains

# Set up list of test data for quick access
dataMelt <- melt(data1, id.vars=c("Subject.number","Sample.year"))
dataMelt$variable <- as.character(dataMelt$variable)
dataMelt$variable <- strain_years[dataMelt$variable]
colnames(dataMelt) <- c("individual","samples","virus","titre")
dataMelt <- dataMelt[complete.cases(dataMelt),]
dataMelt <- dataMelt[order(dataMelt$individual, dataMelt$samples, dataMelt$virus),]
finalDat <- plyr::ddply(dataMelt,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))
finalDat <- finalDat[order(finalDat$individual, finalDat$run, finalDat$samples, finalDat$virus),]
write.table(finalDat,"~/Documents/Fluscape/serosolver/data/real/vietnam_data_primary.csv",sep=",",row.names=FALSE)
