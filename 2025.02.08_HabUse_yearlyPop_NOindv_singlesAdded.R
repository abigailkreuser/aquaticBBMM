## Calculates the summed estimated occurrence distributions to estimate 
# the amount of time or whale days spent in a grid cell. 
## Takes a longgggg  while

library(tidyverse)
library(raster)
library(sf)
library(terra)
library(doParallel)


todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 


rw_areas <- st_read('data/rw_polygons/20241129_rw_areasGSL5_Expanded.shp')
## check to make sure desired polys look as desired 
# ggplot() + 
#   geom_sf(data=rw_areas, aes(fill=name), alpha=0.5, show.legend = "polygon")  #Set alpha at 0.5 for transparency fill= polyCol ,

rw_lats <- st_read('data/rw_polygons/20240420_rw_LatZones.shp') 
#colnames(rw_lats) <- c("zone", "geometry")


#####################################################################################
#could maybe make this be ordered north to south here so I dont have to do it in post 

daysInaYear <-  format(as.Date(as.Date("2024-01-01"):as.Date("2024-12-31"), origin="1970-01-01"), "%b%d")
`%ni%` <- Negate(`%in%`)

area_sums_df <- lat_sums_df<- indv_area_sums_df <- indv_lat_sums_df <- missingTrax <- NULL

#file path to the daily rasters
#local route 
#dailyRastPath <- "data/processedData/results/20240619YearDaily_Rasters/"
#super computer route
dailyRastPath <- "data/processedData/results/20250207YearDaily_Rasters_singlesAdded/"
year_folders <- list.files(dailyRastPath)

g <-1
for(g in 1:length(year_folders)){
  
  #bring in the raster files of that year
  files_ready = list.files(paste0(dailyRastPath, year_folders[g], '/'), full.names=TRUE)# read in file names
  file_day_names = list.files(paste0(dailyRastPath, year_folders[g], '/'))
  #read in each raster 
  r_daily <- terra::rast(files_ready)
  varnames(r_daily)<- sub('\\.nc$', '', file_day_names) 
  names(r_daily)<- sub('\\.nc$', '', file_day_names)
  dayName <- varnames(r_daily)
  
  ## areas habitat utilization 
  #extract points and sum values in each polygon
  area_sums <- terra::extract(r_daily, rw_areas) %>%
    group_by(ID) %>% 
    summarize(across(where(is.numeric),\(x) sum(x, na.rm = TRUE))) 
  
  #rename ID to the names of polygons
  area_sums$ID <- rw_areas$name

  #prep to add 0s for days with no data
  noDataDays <-  daysInaYear[daysInaYear %ni% colnames(area_sums)] #made up opposite of %in% opperator
  noDataDays <- noDataDays[noDataDays !="Feb29"] #hopefully if theres a leap year its already acounted for?? I dont want to add 0s for lack of one
  
  noDataDays_area_df <- as.data.frame(matrix(data=0, nrow = length(rw_areas$name), ncol= length(noDataDays)))
  colnames(noDataDays_area_df) <- noDataDays
  
 
  #add noDataDays to the larger area_sums
  area_sums<- cbind(area_sums, noDataDays_area_df)
  
  
  # transpose to be days as rows
  area_sums <-area_sums %>%
    t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column() 
  
  #rename and remove the top row to make it the header
  colnames(area_sums) <- c("day", area_sums[1,-1])
  area_sums <- area_sums[-1,]
  area_sums <- area_sums %>%
    mutate_at(2:dim(area_sums)[2], as.numeric) 
  
  #add julian date columns, year, and daily total
  area_sums <- area_sums %>%
    mutate(totalinPolys = rowSums(area_sums[,2:dim(area_sums)[2]]),
           j_date = yday(as.Date(area_sums$day, "%b%d")),
           year = paste0(year_folders[g]))
  
  #calculate habitat use outside of polygons
  dailyTotals <- as.data.frame(terra::global(r_daily, "sum", na.rm=TRUE)) %>%
    tibble::rownames_to_column()
 
  colnames(dailyTotals) <- c("day", "nWhales") #daily total could also be n whales modeled!!! 
  
  #add column to the df 
  area_sums <- left_join(area_sums, dailyTotals, by="day")
  
  #change NAs to 0s
  area_sums$nWhales[is.na(area_sums$nWhales)] <- 0
  
  #calculate the amount of time spent outside those polygons
  area_sums$other <- area_sums$nWhales - area_sums$totalinPolys
  
  ##latitudinal zones 
  #extract points and sum values in each polygon
  lat_sums <- terra::extract(r_daily, rw_lats) %>% 
    group_by(ID) %>% 
    summarize(across(where(is.numeric),\(x) sum(x, na.rm = TRUE))) 
  
  #rename ID to the names of polygons
  lat_sums$ID <- rw_lats$name
  
  #prep lat noDataDays
  noDataDays_lat_df <- as.data.frame(matrix(data=0, nrow = length(rw_lats$name), ncol= length(noDataDays)))
  colnames(noDataDays_lat_df) <- noDataDays
  
  
  #add noDataDays to the larger lat_sums
  lat_sums<- cbind(lat_sums, noDataDays_lat_df )
  
  # transpose to be days as rows
  lat_sums <-lat_sums %>%
    t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column() 
  
  #rename and remove the top row to make it the header
  colnames(lat_sums) <- c("day", lat_sums[1,-1])
  lat_sums <- lat_sums[-1,]
  lat_sums <- lat_sums %>%
    mutate_at(2:dim(lat_sums)[2], as.numeric) 
  
  #add julian date columns, year, and daily total
  lat_sums <- lat_sums %>%
    mutate(totalinPolys = rowSums(lat_sums[,2:dim(lat_sums)[2]]),
           j_date = yday(as.Date(lat_sums$day, "%b%d")),
           year = paste0(year_folders[g]))
  
  #add nWhales column to the df 
  lat_sums <- left_join(lat_sums, dailyTotals, by="day")
  
  #change NAs to 0s
  lat_sums$nWhales[is.na(lat_sums$nWhales)] <- 0
  
  #calculate the amount of time spent outside those polygons
  lat_sums$other <- lat_sums$nWhales - lat_sums$totalinPolys

  area_sums_df<- rbind(area_sums_df, area_sums)
  lat_sums_df <- rbind(lat_sums_df,lat_sums)

}


## add time period to df for ease 
area_sums_df <-   area_sums_df %>%
  dplyr::mutate(time_period = case_when(year < 1980 ~ "pre_1980",
                                        year >=1980 & year < 1990 ~ "1980to1989", #previously was 1978to1991 but making it more cohesive 
                                        year >=1990 & year <2000 ~ "1990to1999",
                                        year >=2000 & year <2010 ~ "2000to2009",
                                        year >= 2010 & year <2015 ~ "2010to2014",
                                        year >= 2015 & year <2020 ~ "2015to2019",
                                        year >= 2020 ~ "2020to2022",))
lat_sums_df <-   lat_sums_df %>%
  dplyr::mutate(time_period = case_when(year < 1980 ~ "pre_1980",
                                        year >=1980 & year < 1990 ~ "1980to1989", #previously was 1978to1991 but making it more cohesive 
                                        year >=1990 & year <2000 ~ "1990to1999",
                                        year >=2000 & year <2010 ~ "2000to2009",
                                        year >= 2010 & year <2015 ~ "2010to2014",
                                        year >= 2015 & year <2020 ~ "2015to2019",
                                        year >= 2020 ~ "2020to2022",))




##Need to save this as a csv

yearlyPath<- "data/processedData/results/habitatUtilizations/yearly/"
write.csv(area_sums_df, paste0(yearlyPath, todayIS, "areaSums_all_habitatUtilization_singlesAdded_1.csv"))
write.csv(lat_sums_df, paste0(yearlyPath, todayIS, "latSums_all_habitatUtilization_singlesAdded_1.csv"))
#########################################################################################
# area_sums_df <-  read.csv("data/processedData/results/habitatUtilizations/yearly/areaSums_all_habitatUtilization.csv")
# lat_sums_df <-read.csv("data/processedData/results/habitatUtilizations/yearly/latSums_all_habitatUtilization.csv")


