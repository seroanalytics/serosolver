Sero <- read.csv("~/Dropbox/PHIRST_Serology_2016/PHIRST_2016_hai_res_long_20Jun19.csv")

PCR_add_dates <- read.csv("~/Dropbox/PHIRST_Serology_2016/PHIRST_2016_PCR_26Feb19_add_dates.csv")
PCR_add_dates$ind_id <- as.character(PCR_add_dates$ind_id)

Demog <- read.csv("~/Dropbox/PHIRST_Serology_2016/PHIRST Demog 5Jul2019.csv")


###############################################################
#########Format PCR data into infection data format#########
###############################################################

#for some positives, they also defined which infection they were likely to have had
#some they didnt know which it was, some said both 
table(rowSums(PCR_add_dates[which(PCR_add_dates$npsinfa == 1), c("npsh3", "npsh1")]))

#which said neither, there are no subtype specific CT values
PCR_add_dates[which(PCR_add_dates$npsinfa == 1 & PCR_add_dates$npsh1 == 0
                    & PCR_add_dates$npsh3 == 0), c("npsinfa", "npsinfact", "npsh3", "npsh3ct", "npsh1", "npsh1ct")]
#which said both, there were high CT values for both h3 and h1
PCR_add_dates[which(PCR_add_dates$npsinfa == 1 & PCR_add_dates$npsh1 == 1
                    & PCR_add_dates$npsh3 == 1), c("npsinfa", "npsinfact", "npsh3", "npsh3ct", "npsh1", "npsh1ct")]

# keep the neither, the both and the just h1
#remove the npsinfa and npsh3 pos and npsh1 neg
PCR_H1N1 <- PCR_add_dates[- which(PCR_add_dates$npsinfa == 1 & PCR_add_dates$npsh3 == 1 &
                      PCR_add_dates$npsh1 == 0), ]
#Extract month/year of sampling
PCR_H1N1$collection_date_no <- as.Date(PCR_H1N1$collection_date_no, format = "%d/%m/%Y")

#Extract the negative and positive results
PCR_H1N1_neg <- PCR_H1N1[which(PCR_H1N1$npsinfa == 0), c("ind_id", "collection_date_no", "npsinfa")]
PCR_H1N1_pos <- PCR_H1N1[which(PCR_H1N1$npsinfa == 1), c("ind_id", "collection_date_no", "npsinfa")]

#Convert id to character 
PCR_H1N1_pos$ind_id <- as.character(PCR_H1N1_pos$ind_id)

#Find the first and last time positive, if all are within 30 days, then call them the same 
min_time_pos <- do.call("c", lapply(split(PCR_H1N1_pos$collection_date_no, PCR_H1N1_pos$ind_id), min))
max_time_pos <- do.call("c", lapply(split(PCR_H1N1_pos$collection_date_no, PCR_H1N1_pos$ind_id), max))

#Combine into data frame to calculate difference
time_pos <- data.frame(min_time_pos, max_time_pos)
time_pos$time_diff <- time_pos$max_time_pos - time_pos$min_time_pos

#any individuals had infections that lasted longer than 30 days?
any(time_pos$time_diff > 30) #No

#NOTE : the second cohort has many more longer infections, it might be worth
#considering changing timescale or using Stefano's method of deciding which infections were separate

PCR_H1N1_pos_update <- data.frame(ind_id = names(min_time_pos), collection_date_no = min_time_pos,npsinfa = rep(1, length(min_time_pos)))
#combine pos and negative
PCR_H1N1_update <- rbind(PCR_H1N1_neg,PCR_H1N1_pos_update) 
#remove the NAs - dates missing
PCR_H1N1_update <- PCR_H1N1_update[!is.na(PCR_H1N1_update$collection_date_no),]

### then convert to the correct scale for sersolver
buckets <- 12 # monthly time scale
year_inf <- as.numeric(format(PCR_H1N1_update$collection_date_no,format="%Y"))*buckets
month_inf <- as.numeric(format(PCR_H1N1_update$collection_date_no,format="%m"))

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

#REMOVE ANY DATES EARLIER THAN 2016, THESE SHOULD NOT EXIST!
PCR_H1N1_update <- PCR_H1N1_update[-which(PCR_H1N1_update$infection_time < 2016*12), ]

#then put into matrix to be passed as infection data
library(tidyr)
length(unique(PCR_H1N1_update$ind_id))*length(unique(PCR_H1N1_update$infection_time))

#expand so that we have NAs for those individuals with no PCR data 
PCR_H1N1_expanded <- PCR_H1N1_update %>% complete(infection_time, ind_id)
dim(PCR_H1N1_expanded)

#there seems to be PCR data for individuals beyond the study design 
#for the Epidemics poster, I'm just going to extract the PCR data from 2016
#which corresponds to cohort 1
PCR_2016 <-PCR_H1N1_expanded[which(PCR_H1N1_expanded$infection_time < 2017*12),]


library(reshape)

infection_data <- cast(PCR_2016, ind_id ~ infection_time)
#ensure wd is serosolver here
#save csv
#write.csv(infection_data, file = "data/south_africa/infection_data_2016.csv", row.names = FALSE)

##########################################################
################Formatting titre data#####################
##########################################################

#add ages (or DOB)

Sero$serum_collection_date_draw_1 <- as.Date(Sero$serum_collection_date_draw_1, format = "%d/%m/%Y")

### then convert to the correct scale for sersolver
buckets <- 12 # monthly time scale
year_sero <- as.numeric(format(Sero$serum_collection_date_draw_1,format="%Y"))*buckets
month_sero <- as.numeric(format(Sero$serum_collection_date_draw_1,format="%m"))

#collection time on correct scale
collection_time <- year_sero + month_sero
#add to data.frame
Sero$collection_time <- collection_time 

#convert to the log scale function
convert_to_log_scale <- function(titre){
  titre[which(titre < 10)] <- 0
  titre[which(titre >= 10)] <- floor(log2(titre[which(titre >= 10)] / 10) + 1)
  return(titre)
}

#Just take H1N1 to begin with 
#convert to log scale
Sero$flua_h1n1pdm_1_log <- convert_to_log_scale(Sero$flua_h1n1pdm_1)
Sero$flua_h1n1pdm_2_log <- convert_to_log_scale(Sero$flua_h1n1pdm_2)

# match id numbers in Sero to master data base Demog to find date of birth (DOB)
DOB <- Demog$dob[match(Sero$ind_id, Demog$ind_id)]

#convert DOB to the serosolver timescale
DOB <- as.Date(DOB, format = "%d/%m/%Y")
buckets <- 12 # monthly time scale
year_DOB <- as.numeric(format(DOB,format="%Y"))*buckets
month_DOB <- as.numeric(format(DOB,format="%m"))
DOB_time <- year_DOB + month_DOB

#combine to make titre_dat

#is run the number of samples
virus <- 24191
titre_dat_1 <- data.frame(individual = Sero$ind_id, samples = collection_time, 
                        virus = virus, titre = Sero$flua_h1n1pdm_1_log , 
                        run = 1 , group = 1, 
                        DOB = DOB_time)

titre_dat_2 <- data.frame(individual = Sero$ind_id, samples = collection_time, 
                        virus = virus, titre = Sero$flua_h1n1pdm_2_log , 
                        run = 2 , group = 1, 
                        DOB = DOB_time)

titre_dat <- rbind(titre_dat_1, titre_dat_2)

# remove these NA's, need to go back and figure out where they came from
titre_dat <- titre_dat[which(!is.na(titre_dat$samples)),]

#ensure wd is serosolver here
#write.csv(titre_dat, file = "data/south_africa/titre_dat_full.csv", row.names = FALSE)

