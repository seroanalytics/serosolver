#use 2017 PCR data to inform attack rate time

#2017 cohort
PCR_add_dates<- read.csv("PHIRST_Serology_2016/PHIRST_2017_PCR_21Mar19.csv")

#add column for id number
PCR_add_dates$ind_id <- substr(PCR_add_dates$individual_visit_study_id, 1, 8)

#for some positives, they also defined which infection they were likely to have had
#some they didnt know which it was, some said both 
table(rowSums(PCR_add_dates[which(PCR_add_dates$npsinfa == 1), c("npsh3", "npsph1")]))

#which said neither, there are no subtype specific CT values
PCR_add_dates[which(PCR_add_dates$npsinfa == 1 & PCR_add_dates$npsph1 == 0
                    & PCR_add_dates$npsh3 == 0), c("npsinfa", "npsinfact", "npsh3", "npsh3ct", "npsph1", "npsph1ct")]
#which said both, there were high CT values for both h3 and h1
PCR_add_dates[which(PCR_add_dates$npsinfa == 1 & PCR_add_dates$npsph1 == 1
                    & PCR_add_dates$npsh3 == 1), c("npsinfa", "npsinfact", "npsh3", "npsh3ct", "npsph1", "npsph1ct")]

# keep the neither, the both and the just h1
#remove the npsinfa and npsh3 pos and npsph1 neg
PCR_H1N1 <- PCR_add_dates[- which(PCR_add_dates$npsinfa == 1 & PCR_add_dates$npsh3 == 1 &
                                    PCR_add_dates$npsph1 == 0), ]
#Extract month/year of sampling
#UPDATED DATES NOT AVAILABLE YET FOR THIS COLUMN
PCR_H1N1$npsdatecol <- as.Date(PCR_H1N1$npsdatecol, format = "%d-%m-%Y")

#Extract the negative and positive results
PCR_H1N1_neg <- PCR_H1N1[which(PCR_H1N1$npsinfa == 0), c("ind_id", "npsdatecol", "npsinfa")]
PCR_H1N1_pos <- PCR_H1N1[which(PCR_H1N1$npsinfa == 1), c("ind_id", "npsdatecol", "npsinfa")]

#Convert id to character 
PCR_H1N1_pos$ind_id <- as.character(PCR_H1N1_pos$ind_id)

#Find the first and last time positive, if all are within 30 days, then call them the same 
min_time_pos <- do.call("c", lapply(split(PCR_H1N1_pos$npsdatecol, PCR_H1N1_pos$ind_id), min))
max_time_pos <- do.call("c", lapply(split(PCR_H1N1_pos$npsdatecol, PCR_H1N1_pos$ind_id), max))

#Combine into data frame to calculate difference
time_pos <- data.frame(min_time_pos, max_time_pos)
time_pos$time_diff <- time_pos$max_time_pos - time_pos$min_time_pos

#any individuals had infections that lasted longer than 30 days?
any(time_pos$time_diff > 30) #Yes
which(time_pos$time_diff > 30) #Only two, I'll leave them in for the poster

#NOTE : the second cohort (this one!) has many more longer infections, it might be worth
#considering changing timescale or using Stefano's method of deciding which infections were separate

PCR_H1N1_pos_update <- data.frame(ind_id = names(min_time_pos), npsdatecol = min_time_pos,npsinfa = rep(1, length(min_time_pos)))
#combine pos and negative
PCR_H1N1_update <- rbind(PCR_H1N1_neg,PCR_H1N1_pos_update) 
#remove the NAs - dates missing
PCR_H1N1_update <- PCR_H1N1_update[!is.na(PCR_H1N1_update$npsdatecol),]

### then convert to the correct scale for sersolver
buckets <- 12 # monthly time scale
year_inf <- as.numeric(format(PCR_H1N1_update$npsdatecol,format="%Y"))*buckets
month_inf <- as.numeric(format(PCR_H1N1_update$npsdatecol,format="%m"))

#infection time on correct scale
infection_time <- year_inf + month_inf
#add to data.frame
PCR_H1N1_update$infection_time <- infection_time

#remove the repeated negatives (which occur from having more than one negative in a month)
PCR_H1N1_update <- unique(PCR_H1N1_update[ , c("ind_id", "infection_time", "npsinfa") ] )
dim(PCR_H1N1_update)

#remove any negative entires for the positive months
#(currently we could have positives and negatives in the same month)
pos <- PCR_H1N1_update[which(PCR_H1N1_update$npsinfa == 1),]

#check that for all ids, there are no multiple infection times entered
ids <- pos$ind_id
ids_check <- c()
for(i in 1:length(ids)){
  if(any(table(PCR_H1N1_update$infection_time[which(PCR_H1N1_update$ind_id == ids[i])])>1)){
    ids_check <- c(ids_check, i)
  }
} 

#remove those with negatives in positive months
which_to_remove <- sapply(1:dim(pos)[1], function (x) which(PCR_H1N1_update$ind_id == pos$ind_id[x]
                                                            & PCR_H1N1_update$infection_time == pos$infection_time[x]
                                                            & PCR_H1N1_update$npsinfa == 0))


PCR_H1N1_update <- PCR_H1N1_update[-which_to_remove, ]

#check there are no multiple entries left
ids <- unique(PCR_H1N1_update$ind_id)
ids_check <- c()
for(i in 1:length(ids)){
  if(any(table(PCR_H1N1_update$infection_time[which(PCR_H1N1_update$ind_id == ids[i])])>1)){
    ids_check <- c(ids_check, i)
  }
} 

# the repeats come from individuals having both a negative and a positive in the same month
#these should have been removed already

table(PCR_H1N1_update$ind_id[which(PCR_H1N1_update$ind_id %in% ids[ids_check]
                      & PCR_H1N1_update$npsinfa == 0)])

#remove the negatives, bias towards positives 
PCR_H1N1_update <- PCR_H1N1_update[-which(PCR_H1N1_update$ind_id %in% ids[ids_check]
                             & PCR_H1N1_update$npsinfa == 0),]

#then put into matrix to be passed as infection data
library(tidyr)
length(unique(PCR_H1N1_update$ind_id))*length(unique(PCR_H1N1_update$infection_time))

#expand so that we have NAs for those individuals with no PCR data 
PCR_H1N1_expanded <- PCR_H1N1_update %>% complete(infection_time, ind_id)
dim(PCR_H1N1_expanded)

#take just one year
PCR_2017 <-PCR_H1N1_expanded[which(PCR_H1N1_expanded$infection_time >= 2017*12
                                   & PCR_H1N1_expanded$infection_time < 2018*12),]

library(reshape)

infection_data <- cast(PCR_2017, ind_id ~ infection_time)

colSums(infection_data, na.rm = T)

no_inf <- colSums(infection_data, na.rm = T)


#save csv
write.csv(infection_data, file = "data/south_africa/infection_data_2017.csv", row.names = FALSE)
