#script pulls data from the data folder and then creates csv's so it is easier to manipulate in the future. 
#I might leave the forloops as is, because it is for the whole data frame rather than a summary stat. I will possibly include those
#last updated: 4/24/2022

#packages
#install.packages("tidyverse") #for a remote machine or new download of R
library(tidyverse)

#########################################################################
#    Load and check data
#########################################################################
ID_sightings_db = read.csv('my_data.csv')


head(ID_sightings_db)
dim(ID_sightings_db)
ID_sightings_db <- ID_sightings_db %>% 
    mutate(SightingTime = ifelse(is.na(SightingTime), 0000, SightingTime))

#make times that are 0 <- 2359 so that they will be ordered last. When 'sliced' for selecting the daily obs then a real time observation will be selected first
ID_sightings_db <- ID_sightings_db %>% 
  mutate(SightingTime = ifelse(SightingTime == 0, 2359, SightingTime))


########################################################################
#    Sighting pairs
#########################################################################

todayIS <- Sys.Date( )
todayIS <- format(todayIS, format="%Y.%m.%d")

tail(ID_sightings_db)

sightings_prep = ID_sightings_db %>%
  dplyr::filter(!is.na(SightingEGNo), SightingYear > 1977) %>% #remove sightings where ID - EGNo is NA, and anything before 1977
  dplyr::mutate(date=as.Date(paste0(SightingYear,'-',SightingMonth,'-',SightingDay))) %>% #revise date column
  dplyr::filter(!is.na(date)) %>%  #Gets rid of sighting records where SightingDay or SightingMonth is 0
  dplyr::arrange(SightingEGNo, date, SightingTime) #arrange by time

sightings_prep<- sightings_prep %>%
  group_by(SightingEGNo, date) %>%
  arrange(SightingEGNo, date, SightingTime) %>% #I think just grouping can mess up the order of the time
  slice(1) #this picks the first sighting of that day 
sightings_prep <- sightings_prep[sightings_prep$Latitude != 0,]
# sightings_prep <- sightings_prep %>%
#   filter(SightingYear >= 1980, SightingYear <= 2022 )

dim(sightings_prep)
min(sightings_prep$SightingYear)
 

sight_pairs = NULL
for (i in seq(dim(sightings_prep)[1]-1)){
  if(sightings_prep$SightingEGNo[i] == sightings_prep$SightingEGNo[i+1])
  {
    sight_date_diff = unclass(difftime(sightings_prep$date[i+1], sightings_prep$date[i], units="days"))
    if (sight_date_diff <= 30) # no more than 30 days apart 
    {
      pair = data.frame(EGNo = sightings_prep$SightingEGNo[i], 
                        GenderCode = sightings_prep$GenderCode[i], 
                        sight1_date = sightings_prep$date[i], 
                        sight2_date = sightings_prep$date[i+1], 
                        sight1_latitude = sightings_prep$Latitude[i],
                        sight1_longitude = sightings_prep$Longitude[i],                        
                        sight2_latitude = sightings_prep$Latitude[i+1],
                        sight2_longitude = sightings_prep$Longitude[i+1],
                        sight1_area = sightings_prep$AreaCode[i], 
                        sight2_area = sightings_prep$AreaCode[i+1],
                        sight1_id = sightings_prep$SightingId[i], 
                        sight2_id = sightings_prep$SightingId[i+1],
                        sight1_ageClassCode = sightings_prep$AgeClassCode[i], 
                        sight2_ageClassCode = sightings_prep$AgeClassCode[i+1])
      sight_pairs = rbind(sight_pairs, pair)
    }
  }
}

#double check tails match
tail(sightings_prep)
tail(sight_pairs)
dim(sight_pairs)
#dim(STASHED)
summary(sight_pairs)

sight_pairs = sight_pairs %>%
  dplyr::mutate(interval_days = unclass(difftime(sight2_date, sight1_date, units="days"))) %>%
  dplyr::mutate(year = format(sight1_date, "%Y")) %>%
  dplyr::mutate(time_period = case_when(year >=1980 & year < 1990 ~ "1980to1989", #previously was 1978to1991 but making it more cohesive 
                                 year >=1990 & year <2000 ~ "1990to1999",
                                 year >=2000 & year <2010 ~ "2000to2009",
                                 year >= 2010 & year <2015 ~ "2010to2014",
                                 year >= 2015 & year <2020~ "2015to2019",
                                 year >= 2020 ~ "2020to2022")) %>% #regime shift year 2015, will need to update timeperiods as data is added 
  dplyr::mutate(movement = paste0(sight1_area, "_to_", sight2_area)) %>%
  dplyr::filter(sight1_latitude !=0, sight2_latitude !=0, sight1_longitude!=0, sight2_longitude!=0) #removing 0 lat/long



## UPDATETHE TIMEPERIODS!!!! 
# 
# mvmnt <- mvmnt %>%
#   mutate(time_period = case_when(year < 1990 ~ "1978to1989", #previously was 1978to1991 but making it more cohesive 
#                                  year >=1990 & year <2000 ~ "1990to1999",
#                                  year >=2000 & year <2010 ~ "2000to2009",
#                                  year >= 2010 & year <2015 ~ "2010to2014",
#                                  year >= 2015 & year <2020~ "2015to2019",
#                                  year >= 2020 ~ "2020to2022"),
#          middle_time_period = case_when(middle_year < 1990 ~ "1978to1989", #previously was 1978to1991 but making it more cohesive 
#                                         middle_year >=1990 & middle_year <2000 ~ "1990to1999",
#                                         middle_year >=2000 & middle_year <2010 ~ "2000to2009",
#                                         middle_year >= 2010 & middle_year <2015 ~ "2010to2014",
#                                         middle_year >= 2015 & middle_year <2020~ "2015to2019",
#                                         middle_year >= 2020 ~ "2020to2022"))
# 
# 
# glimpse(mvmnt)


dim(sight_pairs[sight_pairs$interval_days == 0,]) #this should be 0 

#checks and summaries 
dim(sight_pairs %>% filter(interval_days > 0)) #this should match the original dimensions above
dim(sight_pairs)
summary(sight_pairs)
sight_pairs %>% dplyr::group_by(time_period) %>% dplyr::summarize(n=n())

#write csv files 
#write.csv(sight_pairs, file=(paste0('data/processedData/', todayIS, '_sight_pairs.csv')), row.names = FALSE)
#head(sight_pairs)
#sight_pairs <- read.csv("data/processedData/20240615_sight_pairs.csv", header = TRUE)
#head(sight_pairs)


# consecutive sightings
#YEAR IS A NON-NUMERIC SO CONVERT BACK TO THAT
sight_pairs$year<- as.numeric(sight_pairs$year)

consecutive <- sight_pairs %>%
  dplyr::group_by(EGNo) %>%
  dplyr::filter(interval_days != 0) %>%
  dplyr::mutate(consecutive = 1) %>%
  dplyr::mutate(decade = floor(year/10)*10)

dim(consecutive)

consecutive <- consecutive[ order(consecutive$EGNo, consecutive$sight1_date),]

for(i in 1:(length(consecutive$EGNo)-1)){#change to match seq of dim 
  if((consecutive$sight1_id[i+1] == consecutive$sight2_id[i])){ #if sight1_id of the next one = sight2_id of the last one
    consecutive$consecutive[i+1] = consecutive$consecutive[i] +1 
  }
}


#double check end 
#should look to make sure the sighting ids zigzag down the column and the number of consecutive days makes sense
summary(consecutive)
glimpse(consecutive)
tail(consecutive[,c("EGNo", "sight1_date", "sight2_date", 
                    "sight1_id", "sight2_id", "interval_days","consecutive")], n=10)

#write.csv(consecutive, file=(paste0('data/processedData/', todayIS, '_consecutive.csv')), row.names = FALSE)

#Getting total consecutive days and grouping 

#reverse by date to prep for the forloop
revConsecutive <- consecutive[order(consecutive$EGNo, consecutive$sight1_date, decreasing = TRUE), ]
#newMu_pts2 <- newMu_pts[dim(newMu_pts)[1]:1,] # this is a smarter way to do reverse order

#add columns I want
revConsecutive$totalChain = 0
revConsecutive$trackID = 0 
head(revConsecutive[,c("EGNo", "sight1_date", "sight2_date", 
                       "sight1_id", "sight2_id", "interval_days","consecutive", "totalChain", "trackID")])
tail(revConsecutive[,c("EGNo", "sight1_date", "sight2_date", 
                       "sight1_id", "sight2_id", "interval_days","consecutive", "totalChain",  "trackID")], n=10)

#doing a count to create a chain, but if there is 6 sightings it counts 5, so down in the dist section below I did a simple fix
i<-1
while(i < (dim(revConsecutive)[1])){
  if((revConsecutive$sight1_id[i] == revConsecutive$sight2_id[i+1])){ #if sight1_id of the next one = sight2_id of the last one
    tmp <- revConsecutive$consecutive[i]
    revConsecutive$totalChain[i:(i+tmp-1)] = tmp
    i <- (i+tmp)
  }else{
    revConsecutive$totalChain[i] = revConsecutive$consecutive[i]
    i <- i+1
  }
}

#This is the ultimate check!! 
table(revConsecutive$totalChain) # all total chains should be able to be divided by the length of the chain and be a whole number. 

table(revConsecutive$totalChain)/as.numeric(names(table(revConsecutive$totalChain)))
#this also gives the number of each chain I have 


#now adding trackID
totalConsecutive <- revConsecutive[order(revConsecutive$EGNo, revConsecutive$sight1_date), ]
#double check 
head(totalConsecutive[,c("EGNo", "sight1_date", "sight2_date", 
                       "sight1_id", "sight2_id", "interval_days","consecutive", "totalChain", "trackID")])
tail(totalConsecutive[,c("EGNo", "sight1_date", "sight2_date", 
                       "sight1_id", "sight2_id", "interval_days","consecutive", "totalChain",  "trackID")], n=10)

#trackID through while loop it is now consecutive 2022-07-06 :)  

i <- 1
while(i < (dim(totalConsecutive)[1])){
  tmp <- totalConsecutive$totalChain[i] 
  if(i > 1){
    prev_track_number <- totalConsecutive$trackID[i-1]
    totalConsecutive$trackID[i:(i+tmp-1)] <- prev_track_number + 1
    i <- (i+tmp)
  }else{
    totalConsecutive$trackID[i:(i+tmp-1)] <- i
    i <- (i+tmp)
  }
}

#this should be true
max(totalConsecutive$trackID) == length(unique(totalConsecutive$trackID))

#but now we add 1 to the total chain, because it is in pairs, and the first sighting doesnt get counted but its an asy fix
totalConsecutive$totalChain <- totalConsecutive$totalChain + 1 #lol I realized this was messed up but easy fix because of the double sightings layout!!!!! 
#rather than doing the whole loop over 

write.csv(totalConsecutive, file=(paste0('data/processedData/', todayIS, 'totalConsecutive.csv')), row.names = FALSE)


###############################################################################################
mvmnt <- read.csv(file="data/processedData/20250203totalConsecutive.csv", header=TRUE) 


head(mvmnt, n=10)
tail(mvmnt)
dim(mvmnt)

sight_dist_m <- raster::pointDistance(p1=mvmnt[,c("sight1_longitude", "sight1_latitude")], p2=mvmnt[,c("sight2_longitude", "sight2_latitude")],
                                      lonlat = TRUE) #elipsoidal distance
summary(sight_dist_m)
table(sight_dist_m==0) #83 entries are 0

#doooo I wanna add these all and in what order.... 
mvmnt$dist_m <- sight_dist_m

# Median Speed 1.3km/h SEUS hain et al
mvmnt$m_per_day <- mvmnt$dist_m/mvmnt$interval_days
mvmnt$dist_km <- sight_dist_m/1000
mvmnt$km_per_day <- mvmnt$dist_km/mvmnt$interval_days
mvmnt$km_per_h <- mvmnt$km_per_day/24
mvmnt$m_per_sec <- mvmnt$m_per_day/86400


head(mvmnt[order(mvmnt$km_per_day, decreasing = TRUE),])

head(mvmnt)

speed <- summary(mvmnt[,c("dist_km", "km_per_day", "km_per_h", "m_per_day", "m_per_sec")])

#write.csv(mvmnt, file= 'data/20250203_mvmnt_with_dist.csv')
##mvmnt <- read.csv(file= 'data/20220805_mvmnt_with_dist.csv')
############
mvmnt[which(mvmnt$dist_m==0),]
length(which(mvmnt$dist_m==0)) #83 observations of whales moving NOWHERE
mvmnt$year[which(mvmnt$dist_m==0)]
#the fact that the lat long is exactly the same is concerning but like that is possible the whales stay in the same area or return 

bb_mvmnt <- mvmnt %>%
  filter(totalChain >2 )
dim(bb_mvmnt)
dim(mvmnt)
length(unique(bb_mvmnt$trackID)) 
#the number of tracks 3 or longer is 8672
write.csv(bb_mvmnt, file= paste0('data/', todayIS, '_mvmnt_with_n_greaterthan2.csv'), row.names = FALSE)

#I do need to determine their areas, so the next loop can run properly 
#########################################
# Go run MoveOrStay Script!!!!
# then load in df that has the areas associated with each sighting. 

mvmnt <- read.csv("data/processedData/20250203mvmnt_withAreasLats_jDates.csv", header = TRUE)

###############################################################
##  making their own track CSVs 
###############################################################
head(mvmnt)
dim(mvmnt)


summary(mvmnt$totalChain)
whaleTrackIDs <- unique(mvmnt$trackID) 
length(whaleTrackIDs) #12486


dim(mvmnt)[1] + length(whaleTrackIDs) # 59615 = number of individual sightings associated with a track for all years 


#if I do this again maybe parallel it 
mypath <- paste0('data/processedData/', todayIS, '_indvTrack_csv/')
dir.create(mypath)
for(j in 1:length(whaleTrackIDs)){
  desiredID <- c(paste0(whaleTrackIDs[j]))
  track_df <- mvmnt[mvmnt$trackID == desiredID, ]
  lastSighting <- track_df[max(track_df$consecutive),] %>%
    dplyr::select(trackID, EGNo, GenderCode, sight2_date, sight2_latitude, sight2_longitude, area_2, lat_2, sight2_id, 
                  consecutive, decade, time_period, jdate_2, season_jdate_2, season_year_2, year_2) %>%
    dplyr::mutate(consecutive=consecutive + 1) %>%
    dplyr::rename(date = sight2_date, lat = sight2_latitude, lon = sight2_longitude, 
                  areaCode =  area_2, latCode = lat_2, sightingID = sight2_id, jdate = jdate_2, 
                  season_jdate = season_jdate_2, season_year = season_year_2, year = year_2)
  daysSincePrev <- c(0, track_df$interval_days) # the first sighting recieves a 0 because there were no prev sightings under constraints 
  fromLast_area <- c("initial", track_df$moveStay_area)
  fromLast_lat <- c("initial", track_df$moveStay_lat)
  km_per_day <- c(0, track_df$km_per_day)
  track_df <- track_df %>%
    dplyr::select(trackID, EGNo, GenderCode, sight1_date, sight1_latitude, sight1_longitude, area_1, lat_1, sight1_id, 
                  consecutive, decade, time_period, jdate_1, season_jdate_1, season_year_1, year_1) %>%
    dplyr::rename(date = sight1_date, lat = sight1_latitude, lon = sight1_longitude, 
                  areaCode = area_1, latCode = lat_1, sightingID = sight1_id,  jdate = jdate_1, 
                  season_jdate = season_jdate_1, season_year = season_year_1, year = year_1) %>%
    dplyr::bind_rows(lastSighting) %>%
    dplyr::bind_cols(daysSincePrev=daysSincePrev, 
                     fromLast_area=fromLast_area,
                     fromLast_lat=fromLast_lat,
                     km_per_day=km_per_day)
  
  #mypath <- paste0('data/processedData/', todayIS, 'indvTrack_csv/')
  write.csv(track_df, file= paste0(mypath, paste('sightingTrack',whaleTrackIDs[j], sep = "_")))
  
}









