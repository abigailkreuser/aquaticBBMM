## how to make rasters for each day of the year that I have 

##I think it should be similar to the way I create the transisiton maps and rasters. 
#Like that df_list of the sight pairs, but do it for days of the year. 
#Years and days of the year will be the plan. 
# 
# https://www.youtube.com/watch?v=GMnuNjnXnS8&ab_channel=MilosMakesMaps
# #this is a really good tutorial and has so many good functions 

# uses terra and speking of which i need to update my code for that distance function will save my life. 

####################
library(tidyverse)
library(sf)
library(doParallel)
library(terra)

#get the date for file saving 
todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 

#where the results will be saved
#super computer route
mypath2daily <- 'data/processedData/results/20250207bbmmRasterDataDaily_res4km_bmvar464.1sqkm_locerr2.5k_singlesAdded/'

file_names<- list.files(mypath2daily)

filez<- as.data.frame(file_names)
file_dates <- as.Date(sub( "(\\d{8})_.*", "\\1", file_names), "%Y%m%d" )

length(unique(file_dates))

#Year Then Daily Raster
filez <- filez %>%
  mutate(date = as.Date(sub( "(\\d{8})_.*", "\\1", file_names), "%Y%m%d" )) %>%
  dplyr::mutate(year = format(date, "%Y")) %>%
  dplyr::mutate(month_day = format(date, "%b%d"))%>%
  dplyr::mutate(decade = case_when(year < 1990 ~ "1978-1989", #previously was 1978-1991 but making it more cohesive 
                               year >=1990 & year <2000 ~ "1990-1999",
                               year >=2000 & year <2010 ~ "2000-2009",
                               year >= 2010 & year <2015 ~ "2010-2014",
                               year >= 2015 & year <2020~ "2015-2019",
                               year >= 2020 ~ "2020-2022"))


#get all the days in a year i have a leap year for the decadal 
daysInaYear <-  format(as.Date(as.Date("2024-01-01"):as.Date("2024-12-31"), origin="1970-01-01"), "%b%d")

#create a list of which files for each day of the year 
dfDailyFiles<- as.data.frame(NULL)
dfDailylist<- as.list(NULL)
dfYearDaily_list<- as.list(NULL)
i <- 1
# can comment out this part once a list is made, 
#and then use it to direct multiple scripts across the years 
for(i in 1:length(unique(filez$year))){
  dfYear <- filez %>%
    filter(year == unique(filez$year)[i])
  for(g in 1:length(daysInaYear)){
    dfDailyFiles <- dfYear  %>%
      filter(month_day == daysInaYear[g])
    dfDailylist[[g]] <- dfDailyFiles
    names(dfDailylist)[g] <- daysInaYear[g]
  }
  dfYearDaily_list[[i]]<- dfDailylist
  names(dfYearDaily_list)[i] <- unique(filez$year)[i]
}
length(dfYearDaily_list)
#this will be a list of the Years within that a list of each day of the year 366days for including leap years
str(dfYearDaily_list)
saveRDS(dfYearDaily_list, paste0("data/processedData/dfYearDaily_list", todayIS, "_singlesAdded.Rdata"))


#If list is already made, just read in here
#dfYearDaily_list <- readRDS("data/processedData/dfYearDaily_list20250207_singlesAdded.Rdata")
#dfYearDaily_list <- readRDS("../IndvMvmntModel/data/processedData/dfYearDaily_list20250207_singlesAdded.Rdata")

##choose the correct exent based on how the figures were made 
extentMax <- rast("data/processedData/20231002_baseRaster_res0.04AquaModis.nc")


#create and set the path where the saving is going to happen for the rasters and for the figures
mypath2daily
mypath2NEWrasters <- paste0('data/processedData/results/', todayIS, 'YearDaily_Rasters_singlesAdded/')
dir.create(mypath2NEWrasters)


selectyears <- length(dfYearDaily_list)
selectyears


## HAVE TO MAKE ALL THE RASTERS FIRST before I can have the scale for visuals

for(j in 1:length(dfYearDaily_list)){ 

  YearTitle <- names(dfYearDaily_list)[j]
  mypath2NEWrasters <- paste0('data/processedData/results/', todayIS, 'YearDaily_Rasters_singlesAdded/', YearTitle,'/') 
  dir.create(mypath2NEWrasters)

  
  dfDailylist <- dfYearDaily_list[[j]]
  for(g in 1:length(dfDailylist) ){
    if(length(dfDailylist[[g]]$file_names) ==0){ #skip days with 0 
      next}
    
    #this is getting pulled from the dfDaily not the files
    dayTitle <- names(dfDailylist)[g]
 
    #bring in the raster files of that day
    files_ready = paste0(mypath2daily, dfDailylist[[g]]$file_names) # read in file names
    #length(files_ready)

    
    r_list <- sapply(files_ready, terra::rast)
    
 
    for(i in 1:length(r_list)){
      r_list[[i]] <- terra::resample(r_list[[i]], extentMax)
    }

    bbmm_stack <- terra::rast(r_list)
    sumDays <- sum(bbmm_stack, na.rm=TRUE)
 
    
    #sumDays <- terra::app(bbmm_stack, fun=base::sum, na.rm=TRUE)

    terra::varnames(sumDays) <- dayTitle
  
  # ### 2. REMOVE 0S saving with 0s so min and maxes can be found
  # ### ------------
  # sumDays <-  terra::ifel(sumDays == 0, NA, sumDays)

#sumDays
  #bb_rast
  #save raster for later computing
  new_file_name <- paste0(mypath2NEWrasters, dayTitle, ".nc")
    terra::writeCDF(sumDays, filename=new_file_name, varname=dayTitle, overwrite=TRUE)

}
}
