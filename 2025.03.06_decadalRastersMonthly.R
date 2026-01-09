## Decadal monthly rasters

## 2025 update, I am using just months, no daily!! 
## remember terra doesn't parallel or some function in it doesnt
####################
library(tidyverse)
library(sf)
library('ncdf4')
library(terra)

#######################################################
#                 Animate rasters with R
#                 Milos Popovic
#                 2023/06/29
########################################################

todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 

#file path to the daily rasters
mypath2years <- "data/processedData/results/20251211YearDaily_Rasters_singlesAdded/"
year_folders<- list.files(mypath2years,full.names = TRUE )

#for this MS I need 80s and 2020-2022 as full years but for the general scope I want the other years/time periods included
time_period_list <- list(c(1978:1989), c(1980:1989), c(1990:1999),c(2000:2009), c(2010:2014), c(2015:2019), c(2010:2019), c(2020:2022), c(2020:2024))
names(time_period_list) <- c("1978to1989", "1980to1989", "1990to1999","2000-2009", "2010-2014", "2015-2019", "2010-2019", "2020-2022", "2020-2024")
length(time_period_list)

daysInaYear <-  format(as.Date(as.Date("2024-01-01"):as.Date("2024-12-31"), origin="1970-01-01"), "%b%d")
`%ni%` <- Negate(`%in%`)
monthsInaYear <-  unique(format(as.Date(as.Date("2024-01-01"):as.Date("2024-12-31"), origin="1970-01-01"), "%b"))


timePeriod_file_list <- as.list(NULL)
DailyFileNames <- as.list(NULL)
tpYearz_daySelected <- NULL

for(i in 1:length(time_period_list)){
  tp_name <- names(time_period_list)[i]
  yearz <- time_period_list[[i]]

  for(j in 1:length(daysInaYear)){
    daySelected <- daysInaYear[j]
    tpYearz_daySelected <- NULL
    for(k in 1:length(yearz)){
      #list the full file name of the day I am looking for.
      #need to use this function for the files that don't actually exist
      filesInYear <- list.files(paste0(mypath2years, yearz[k],'/'), full.names = TRUE)

      DayFileName <- grep(daySelected, filesInYear, value=TRUE)

      tpYearz_daySelected <- c(tpYearz_daySelected, DayFileName)
    }
    DailyFileNames[[j]]<- tpYearz_daySelected
    names(DailyFileNames)[j] <- daySelected

}
    timePeriod_file_list[[i]] <- DailyFileNames
   names(timePeriod_file_list)[i] <- tp_name

}




#save this list
saveRDS(timePeriod_file_list, paste0("data/processedData/timePeriod_file_list", todayIS,"_noSingles", ".Rdata"))
#timePeriod_file_list <- readRDS("data/processedData/timePeriod_file_list20250306_singlesAdded.Rdata")

#check list length
length(timePeriod_file_list) #should be 6 
length(time_period_list)

extentMax <- terra::rast("data/processedData/20231002_baseRaster_res0.04AquaModis.nc")

#create and set the path where the saving is going to happen for the rasters and for the figures
mypath2NEWrasters <- paste0('data/processedData/results/', todayIS, '_decadalRasters_noSingles/')
dir.create(mypath2NEWrasters)


monthlyNEWrasters <- paste0(mypath2NEWrasters, 'monthly/')
dir.create(monthlyNEWrasters)


j<- 1

for(j in 1:length(timePeriod_file_list)){ 
  tp_Title <- names(timePeriod_file_list)[j]
  
  monthly_tp_path<- paste0(monthlyNEWrasters, 
                            tp_Title,'/')
  dir.create(monthly_tp_path)
  
  dfDailylist <- timePeriod_file_list[[j]]
  dfDailyUnlist<- unlist(dfDailylist, use.names = FALSE)

  for(g in 1:length(monthsInaYear)){
    monthlyList <- dfDailyUnlist[grep(monthsInaYear[g], dfDailyUnlist)]
    if(length(monthlyList)==0){
      next}
    
    monthTitle <- monthsInaYear[g]

    #bring in the raster files of that transition area
    files_ready = monthlyList # read in file names
    #length(files_ready)
    
    #all raster files listed
    r_list <- sapply(files_ready, terra::rast)
    
    #resample to max extent
    for(i in 1:length(r_list)){
      crs(r_list[[i]]) <- "+init=epsg:4326"
      r_list[[i]] <- terra::resample(r_list[[i]], extentMax)
    }
    
    #stack as spatrasters 
    bbmm_stack <- terra::rast(r_list)
    #sum
    sumDays <- terra::app(bbmm_stack, fun=base::sum, na.rm=TRUE)
    #label
    terra::varnames(sumDays) <- "days"
    
    # ### 2. REMOVE 0S saving with 0s so min and maxes can be found
    # ### ------------
    # sumDays <-  terra::ifel(sumDays == 0, NA, sumDays)

    #save raster for later computing
    new_file_name <- paste0(monthly_tp_path, monthTitle, ".nc")
    terra::writeCDF(sumDays, new_file_name, varname=monthTitle, overwrite=TRUE)

  }
}
