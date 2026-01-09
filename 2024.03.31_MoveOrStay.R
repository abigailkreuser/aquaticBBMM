### Assigning areas and lat zones to each sighting to determine if the whale moves or stays. 
library(tidyverse)
library(mapdata)  # map_data worldHires coastline
library(raster)
library(sf)
library(sfnetworks)
library(ggnewscale)

todayIS <- Sys.Date()
todayIS <- format(todayIS, format="%Y.%m.%d")

#areas with GSL divided into 5 SLE, neGSL, nwGSL, sGSL and CS
rw_areas <- st_read('data/rw_polygons/20241129_rw_areasGSL5_Expanded.shp')
colnames(rw_areas) <- c("area", "geometry")


sf_use_s2(TRUE)

colnames(rw_areas) <- c("area", "geometry")


#######################################################
#mvmnt <- mvmnt_areas
dim(mvmnt)
glimpse(mvmnt)

mvmnt_sf1 <- sf::st_as_sf(mvmnt, coords=c('sight1_longitude', 'sight1_latitude'), crs=4326)
head(mvmnt_sf1)

mvmnt_sf2 <- sf::st_as_sf(mvmnt, coords=c('sight2_longitude', 'sight2_latitude'), crs=4326)
head(mvmnt_sf2)
###########################
#mvmnt_sf1$polygonName <- st_contains(mvmnt_sf1, rw_areas)

# area_mat<-st_intersects(mvmnt_sf2$geometry, rw_areas$geometry)
# 
# as.integer(area_mat) #if I get an error here there is overlap and it needs to be found otherwise i dont know how to get around it 
# 
# 
# #if the list cannot be forced into an interger, do this to find overlapping polygons
# area_mat_l <- lengths(st_intersects(mvmnt_sf1$geometry, rw_areas$geometry)) >1
# length(area_mat[area_mat_l])

#####################################################
mvmnt_sf1 <- mvmnt_sf1 %>% mutate(
  int = as.integer(st_intersects(mvmnt_sf1, rw_areas)),
  area_1 = if_else(is.na(int), 'NA', rw_areas$area[int])
) 

mvmnt_sf2 <- mvmnt_sf2 %>% mutate(
  int = as.integer(st_intersects(mvmnt_sf2, rw_areas)),
  area_2= if_else(is.na(int), 'NA', rw_areas$area[int])
  )

####################################################
sf_use_s2(TRUE)

area_list = st_intersects(mvmnt_sf1, rw_areas)

int<-NULL
hello <- NULL
for(i in 1:length(area_list)){
  if(length(area_list[[i]]) == 1){
    int <- c(int, as.integer(area_list[i]))
  }else{
    hello <- c(hello, i)
    int <- c(int, as.integer(area_list[[i]][1]))
  }
  
}

#checking for doubles 
mvmnt_sf1[hello,]
area_list[hello]

mvmnt_sf1$int <- int
mvmnt_sf1 <- mvmnt_sf1 %>% mutate(
  area_1 = if_else(is.na(int), 'NA', rw_areas$area[int])
) 


area_list = st_intersects(mvmnt_sf2, rw_areas)

int<-NULL
hello <- NULL
for(i in 1:length(area_list)){
  if(length(area_list[[i]]) == 1){
    int <- c(int, as.integer(area_list[i]))
  }else{
    hello <- c(hello, i)
    int <- c(int, as.integer(area_list[[i]][1]))
  }
  
}

mvmnt_sf2$int <- int
mvmnt_sf2 <- mvmnt_sf2 %>% mutate(
  area_2 = if_else(is.na(int), 'NA', rw_areas$area[int])
) 

#checking for doubles 
mvmnt_sf2[hello,]
area_list[hello]

####################################################


## LATS NEED SPHERICAL GEOM OFF
# if you try to assign grid cells
# to a lat zone to make sure it looks right on the maps, I turn off spherical geometry 
# otherwise the line is drawn like a parabola across the ocean, rather than along
# the latitude band. 

rw_lats <- st_read('data/rw_polygons/20240420_rw_LatZones.shp') 
#colnames(rw_lats) <- c("lat", "geometry")

#polygons to rw_areas
rw_areas <- rw_lats

colnames(rw_areas) <- c("area", "geometry")


### load in lats 

sf_use_s2(FALSE)


# mvmnt_sf1 <- mvmnt_sf1 %>% mutate(
#   int = as.integer(st_intersects(mvmnt_sf1, rw_areas)),
#   lat_1 = if_else(is.na(int), 'NA', rw_areas$area[int])
# ) 
# 
# mvmnt_sf2 <- mvmnt_sf2 %>% mutate(
#   int = as.integer(st_intersects(mvmnt_sf2, rw_areas)),
#   lat_2 = if_else(is.na(int), 'NA', rw_areas$area[int])
#   )


lat_list = st_intersects(mvmnt_sf1, rw_areas)
as.integer(lat_list)
hello<- NULL
int<-NULL
for(i in 1:length(lat_list)){
  if(length(lat_list[[i]]) == 1){
    int <- c(int, as.integer(lat_list[i]))
  }else{
    hello <- c(hello, i)
    int <- c(int, as.integer(lat_list[[i]][1]))
  }

}

lat_list[hello]
mvmnt_sf1[hello[c(1,3,9)],]

dim(mvmnt_sf1)

mvmnt_sf1$int <- int
mvmnt_sf1 <- mvmnt_sf1 %>% mutate(
  lat_1 = if_else(is.na(int), 'NA', rw_areas$area[int])
) 

lat_list = st_intersects(mvmnt_sf2, rw_areas)

int<-NULL
for(i in 1:length(lat_list)){
  if(length(lat_list[[i]]) == 1){
    int <- c(int, as.integer(lat_list[i]))
  }else{
    int <- c(int, as.integer(lat_list[[i]][1]))
  }
  
}


mvmnt_sf2$int <- int
mvmnt_sf2 <- mvmnt_sf2 %>% mutate(
  lat_2 = if_else(is.na(int), 'NA', rw_areas$area[int])
) 

head(mvmnt_sf2)
tail(mvmnt_sf1)
######################################################3
mvmnt_sf1$geometry <- NULL
mvmnt_sf2$geometry <- NULL

# need top make data frame 
mvmnt_withAreas <- left_join(mvmnt_sf1, mvmnt_sf2[,c("sight1_id", "sight1_latitude", "sight1_longitude", "lat_2", "area_2")], by = "sight1_id")


dim(mvmnt_withAreas)
glimpse(mvmnt_withAreas)


### get rid of these silly columns
drop_cols <- c("sight1_area", "sight2_area")

mvmnt_withAreas_n_lats <- mvmnt_withAreas %>% 
  dplyr::select(-one_of(drop_cols))


#write.csv(as.data.frame(mvmnt_withAreas_n_lats), file = paste0("data/processedData/", todayIS, "mvmnt_withAreas_plusLats.csv"), row.names = FALSE)

## next step will be move or stay
unique(mvmnt_withAreas_n_lats$area_1)

mvmnt_withAreas_n_lats$area_1[mvmnt_withAreas_n_lats$area_1 == 'NA'] <- "other"
mvmnt_withAreas_n_lats$area_2[mvmnt_withAreas_n_lats$area_2 == 'NA'] <- "other"

mvmnt_withAreas_n_lats$lat_1[mvmnt_withAreas_n_lats$lat_1 == 'NA'] <- "other"
mvmnt_withAreas_n_lats$lat_2[mvmnt_withAreas_n_lats$lat_2 == 'NA'] <- "other"


#stay or move 

head(mvmnt_withAreas_n_lats)


mvmnt <- mvmnt_withAreas_n_lats
mvmnt <- mvmnt %>%
  mutate(moveStay_lat = case_when(lat_1 == lat_2 ~ "stay",
                               lat_1 != lat_2 ~ "move")) %>%
  mutate(moveStay_area = case_when(area_1 == area_2 ~ "stay",
                                   area_1 != area_2 ~ "move"))   

# if using the singles added df to make more sense for everything <3
mvmnt[mvmnt$sight1_date == mvmnt$sight2_date,]

mvmnt$interval_days[mvmnt$sight1_date == mvmnt$sight2_date] <- 0
mvmnt$consecutive[mvmnt$sight1_date == mvmnt$sight2_date] <- 1
mvmnt$moveStay_lat[mvmnt$sight1_date == mvmnt$sight2_date] <- "single"
mvmnt$moveStay_area[mvmnt$sight1_date == mvmnt$sight2_date] <- "single"


head(mvmnt)
 
#Changing to julian date and then to shift to nov and then label which calving season 
mvmnt <- mvmnt %>%
  mutate(jdate_1 = 0, jdate_2 = 0,  season_jdate_1=0, season_jdate_2=0, season_year_1=0,  season_year_2=0) %>%
  mutate(year_1 = as.numeric(format(as.Date(sight1_date), "%Y"))) %>%
  mutate(year_2 = as.numeric(format(as.Date(sight2_date), "%Y"))) 
  
i<-1

for(i in seq(dim(mvmnt)[1])){
  mvmnt$jdate_1[i] = julian(as.Date(mvmnt$sight1_date[i]), 
                                         origin=as.Date(paste0((mvmnt$year_1[i]-1), '-', 12 , '-',31)))
  if(mvmnt$jdate_1[i] >304){
    mvmnt$season_jdate_1[i] <- mvmnt$jdate_1[i]-304
    mvmnt$season_year_1[i] <- mvmnt$year_1[i]+1
  }else{
    mvmnt$season_jdate_1[i] <- 365-304+mvmnt$jdate_1[i] #not counting leap years o_0
    mvmnt$season_year_1[i] <- mvmnt$year[i] #r starts jdate at 0 so 1 is jan2, and that feels wrong
    
  }
  mvmnt$jdate_2[i] = julian(as.Date(mvmnt$sight2_date[i]), 
                            origin=as.Date(paste0((mvmnt$year_2[i]-1), '-', 12 , '-',31)))
  if(mvmnt$jdate_2[i] >304){
    mvmnt$season_jdate_2[i] <- mvmnt$jdate_2[i]-304
    mvmnt$season_year_2[i] <- mvmnt$year_2[i]+1
  }else{
    mvmnt$season_jdate_2[i] <- 365-304+mvmnt$jdate_2[i] #not counting leap years o_0
    mvmnt$season_year_2[i] <- mvmnt$year_2[i] #r starts jdate at 0 so 1 is jan2, and that feels wrong
    
  }
}

head(mvmnt)

any(is.na(mvmnt$area_1))
any(mvmnt$season_year_1 != mvmnt$year_2)

mvmnt_years <- mvmnt[mvmnt$season_year_1 != mvmnt$season_year_2,]

table(mvmnt_years$area_2)
table(mvmnt_years$moveStay_lat)
table(mvmnt_years$moveStay_area)
table(mvmnt_years$lat_1)

##########
#write.csv(as.data.frame(mvmnt), file = paste0("data/processedData/", todayIS, "_mvmnt_singlesAdded_withAreasLats_jDates.csv"), row.names = FALSE)
