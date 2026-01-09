## The aquatic BBMM


####################
library(tidyverse)
library(mapdata)  # map_data worldHires coastline
library(maps) #idk how the indexing of the land map works, but loading the library mapdata, and lexical scoping maps, I get the map I want
library(raster)
library(sf)
library(sfnetworks)
library(doParallel)
library(ncdf4)
library(terra)
#####################


##############################
##### date labeling prep #####
##############################
todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 


##############
## LAND HO ###
###############
#load in the land, crop, remove country boundaries
mapBase <-st_as_sf(maps::map("worldHires", fill=TRUE, plot=F))
mapBase<- st_cast(mapBase, "MULTIPOLYGON") %>% st_cast("POLYGON")
land_all <- sf::st_transform(mapBase, crs=4326)

#Cropping to north Atlantic study area (open to changes on bounds)
sf_use_s2(FALSE) 

## add a 0 width buffer to allow the crop function. Otherwise it doesnt work because of some self intersection issue.
land <- st_buffer(land_all,dist=0)
lat_bounds= c(20,75)
lon_bounds= c(-110,30)
box <- c(xmin= lon_bounds[1], xmax= lon_bounds[2], ymin= lat_bounds[1], ymax= lat_bounds[2])
land_crop <- st_crop(land, st_bbox(box))
#plot(land_crop) #love the auto colors 
land <- st_union(land_crop) # one land mass without the country boundaries 
land <- st_union(land) #now a cropped simple polygon of points 

sf::sf_use_s2(TRUE) #turn spherical geometry back on
###############################################################################################

######################################################
#####            read in network                  ####
######################################################
#read in network
sightNet_clean <- st_read('data/processedData/20240617_sightNet_clean.nc')
##shapefile for old iMac at home, idk why but was not able to read .nc files
#sightNet_clean <- st_read('data/processedData/sightNet_clean.shp')

head(sightNet_clean)
net <- sfnetworks::as_sfnetwork(sightNet_clean, directed=FALSE)
head(net)
sightNet_clean <- net %>%
  activate("edges") %>%
  mutate(weight = edge_length()) #weighting the edges as much as their length. now theres a column for them, and shortest path will be taken.  

##############################################################################################################################
###             load functions                    #####
#######################################################
make_lineAroundLand<- function(p1, p2, net = sightNet_clean){
  ## Associate the points that have the intersection with the network by finding the nearest node
  fromPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(sightNet_clean,"nodes")),
    sf::st_coordinates(p1),k=1)$nn.idx[,1]
  
  toPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(sightNet_clean,"nodes")),
    sf::st_coordinates(p2),k=1)$nn.idx[,1]
  
  if(fromPoint == toPoint){
    stop(c(unique(singleTrack$trackID), 'nearest network neighbor for from and to point is identical'))
  }else{
    #build the shortest route 
    selected_edges <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(sightNet_clean,"edges"),
                                                    from = fromPoint,
                                                    to = toPoint,
                                                    weights = NULL,
                                                    output = "epath")$epath)
    #make a tibble from it 
    edge_geom <- lapply(selected_edges, function(e) {
      g <- sightNet_clean %>% as_tibble(spatial = FALSE)
      g[e,"geometry"]
    })
    
    #make the tibble a df, then an sf df
    df_edges <- do.call(rbind, edge_geom)
    df_edges_sf <- sf::st_as_sf(df_edges, crs = 4326)
    
    #make that df multiple linestrings 
    edge_l <- df_edges_sf %>% 
      dplyr::summarise(do_union = FALSE) %>% 
      sf::st_cast('LINESTRING')
    
    #make it a singular line string
    edge_l <- edge_l %>% 
      sf::st_sf() %>% 
      sf::st_combine() %>% 
      sf::st_sf() %>% 
      sf::st_line_merge() %>% 
      sf::st_cast("LINESTRING") 
    
    return(edge_l)
  }
}

# return points around land

make_newMuPts <- function(p1, p2, net = sightNet_clean, n_muPts = stepz){
  ## Associate the points that have the intersection with the network by finding the nearest node
  fromPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(sightNet_clean,"nodes")),
    sf::st_coordinates(p1),k=1)$nn.idx[,1]
  
  toPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(sightNet_clean,"nodes")),
    sf::st_coordinates(p2),k=1)$nn.idx[,1]
  
  if(fromPoint == toPoint){
    stop(c('j = ', unique(singleTrack$trackID), ' i=', i, ' nearest network neighbor point is identical'))
  }else{
    #build the shortest route 
    selected_edges <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(sightNet_clean,"edges"),
                                                    from = fromPoint,
                                                    to = toPoint,
                                                    weights = NULL,
                                                    output = "epath")$epath)
    #make a tibble from it 
    edge_geom <- lapply(selected_edges, function(e) {
      g <- sightNet_clean %>% as_tibble(spatial = FALSE)
      g[e,"geometry"]
    })
    
    #make the tibble a df, then an sf df
    df_edges <- do.call(rbind, edge_geom)
    df_edges_sf <- sf::st_as_sf(df_edges, crs = 4326)
    
    #make that df multiple linestrings 
    edge_l <- df_edges_sf %>% 
      dplyr::summarise(do_union = FALSE) %>% 
      sf::st_cast('LINESTRING')
    
    #make it a singular line string
    edge_l <- edge_l %>% 
      sf::st_sf() %>% 
      sf::st_combine() %>% 
      sf::st_sf() %>% 
      sf::st_line_merge() %>% 
      sf::st_cast("LINESTRING") 
    
    # sample the line.... the other parte was done before this, 
    # so i should double check if all the code before is needed
    l_pts <- sf::st_line_sample(st_transform(edge_l, 3857), n = n_muPts)
    
    
    
    new_pts <- sf::st_transform(l_pts, 4326) %>%
      sf::st_cast('POINT') %>%
      sf::st_as_sf() 
    
    
    return(new_pts)
  }
}

## sigma_sq, can use for both approaches 
get_sigma_sq <- function(intervalDays=time_lag[i], alpha, BMvar = BMvar, loc_err1=location_err[i], loc_err2=location_err[i+1]){
  sigma_squared <- intervalDays*alpha*(1-alpha)*BMvar + (1-alpha)^2*loc_err1^2 + alpha^2*loc_err2^2
  return(sigma_squared)
}


##get the segments that cross land

get_lineSegs <- function(p1, p2){
  l_temp <- sf::st_union(p1, p2) %>% #make the two points a linestring!! 
    sf::st_cast('LINESTRING') 
  noLand  <- l_temp %>%
    sf::st_intersects(land, sparse=FALSE)==FALSE
  return(noLand)
}

################################################################################
#### base raster and extent read in 
landNAs <- terra::rast("data/baseRaster/20240619_Sat4kmresRasterLandNAs_forGridDist.nc")
newExtent <- terra::ext(-100, -40, 20, 70)
watergrid_landNAs <- terra::crop(landNAs, newExtent)

# prep labeling
todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 


######################
## Variance plug in ##
######################

#set the motion variance to be using for all tracks. 
BMvarALL <- 464099334 ## THIS ONE, this is the motion var with  2.5 km location err

## location error hard coded in 
locERR <-  2500

# BMvarALL <- 464301909 this one is from the super comupter it still took 17hrs to calculate, and its slightly smaller, so i want to go with the larger one, it also removed 0s and NAs
BMvarALL_sqkm <- round(BMvarALL/(1000^2), digits=1)
BMvarALL_sqrt <- round(sqrt(BMvarALL/(1000^2)), digits=1)
print("BMvarALL_sqkm")
BMvarALL_sqkm
print("BMvarALL_sqrt")
BMvarALL_sqrt
######################################################
### create folders and set paths and variables   #####
#######################################################
#where the results will be saved
mypath2daily <- paste0("data/processedData/results/", todayIS, "bbmmRasterDataDaily_res4km_bmvar", BMvarALL_sqkm, "sqkm_locerr2.5k/")
dir.create(mypath2daily)

mypath2fulltracks <- paste0("data/processedData/results/", todayIS, "bbmmRasterDataFullTracks_res4km_bmvar", BMvarALL_sqkm, "sqkm_locerr2.5k/")
dir.create(mypath2fulltracks)


iteration_j <-  iteration_i <- err_msg <- trackID_j <- NULL


#######################
### whale track ids ###
#######################

mvmnt <- read.csv(file = "data/20250203_mvmnt_with_dist.csv")
#mvmnt <- read.csv(file="data/tracks_thatCrossLandtoRerun.csv")

desiredTrack <- unique(mvmnt$trackID) #prepping all individual tracks to be indexed in the loop by their track ID
# wanted<- c(3460, 2360, 3620, 2503, 4715, 3560, 5060)
#   #read.csv("data/DONTHAVES_fullday_02.06.csv")
# desiredTrack <- unique(wanted$x)

## when running all tracks break up into sections 
#desiredTrack <- desiredTrack[615:length(desiredTrack)]
# length(desiredTrack)
# breaksDesired <- round(seq(from = 1, to = length(desiredTrack), length.out=3))
# running_set <- 1
# desiredTrackStart <- 4 #breaksDesired[running_set] +1 #adds one to the break point, so the track isn't computed twice 
# desiredTrackEnd <- breaksDesired[running_set+1]
# print(c("running set #", running_set, ":", desiredTrackStart, "to",
#         desiredTrackEnd ))


####################################################
#run it 
#can't parallel terra  
for(j in 1:length(desiredTrack)){  #in an ideal world 
  trackpath <- "data/processedData/2025.02.04_indvTrack_csv/"
  singleTrack <- read.csv(paste0(trackpath, paste('sightingTrack',desiredTrack[j], sep = "_")))
  #then convert to sf, for intersections with and routing around land
  singleTrack_sf <- sf::st_as_sf(singleTrack, coords = c('lon','lat'), crs=4326)
  
  n_locs <- dim(singleTrack)[1] # number of rows or known locations 
  time_lag <- singleTrack$daysSincePrev[-1] #Takes out the 0 at the top, but will have (n_locs-1) values
  T_Total <- sum(time_lag) #total time in the track
  singleTrack$location_err <-locERR  #adding this location error, because 2.5km is tha mean distance whales sighted multiple times in one day moved (measured the full track)
  # The bmvar above has been calculated with this motion var 
  #tagged whales also have 2500m as their largest error, but I need to get the notes section from phill
  location_err <- singleTrack$location_err 
  lon <- singleTrack$lon
  lat <- singleTrack$lat
  max_lag = 30 #days no lags greater than 30 days, this should be taken care of in the pull and shape script 
  BMvar <- BMvarALL #assigned BMvar, possibly going to do seasonal in the future or behavior 
  
 
  wiggleroom <- 10
  #these vaules would be specific to the group of tracks or specific tracks I am using.
  # the model goes grid cell by gird cell to calculate the probability the whale was in that cell. Can be very slow if they stay in BOF the whole track 
  # but probability is for the whole North Atlantic 
  x_min <- floor(min(lon)) - wiggleroom 
  x_max <- ceiling(max(lon)) + wiggleroom
  y_min <- floor(min(lat)) - wiggleroom
  y_max <-  ceiling(max(lat)) + wiggleroom 
  
#crop to the newExtent with wiggleroom, for the watergrid
  newExtent <- terra::ext( x_min, x_max, y_min, y_max)
  watergrid_landNAs <- terra::crop(landNAs, newExtent)
  
  probability <- NULL # this is a vector not the 
  int <- NULL
  
  #this one is used to create the daily rasters
  steps_theta <-NULL
  
  #using the sf data frame to determine which segments cross land. 
  sightPoints <- singleTrack_sf$geometry
  P1 <- sightPoints[-length(sightPoints)]
  P2 <- sightPoints[-1]
  stay_in_water <- get_lineSegs(P1,P2) #make sure spherical geom is back on post land crop
  
  #testing 8 steps for the isotope whales if they have 1 day
  if(any(singleTrack_sf$km_per_day >= 100) | any(singleTrack_sf$daysSincePrev == 1)){ # run 8 stepz for whales that are traveling over 100 km in a day
    # when there is 4 steps for the whales that are traveling too fast, however, trying to split that up from 
    # just the length and not make it an alegbraic mess, I am just doing 8 steps for the whole track. 
    
    for(i in 1:(n_locs-1)){
      if(time_lag[i] <= max_lag){
        theta <- NULL
        tm <- 0
        #previously was using 25, which the only thing that im worried about is whales that travel fast, my have clumpy expected distributions, and it not smooth
        #stepz <- 25 #this is increments of the whole time lag aka interval days so that would be 25 steps for 1 day or 25 steps for 25 days (1 a day)
        stepz <- (4+(time_lag[i]-1)*8+4) # this is written ugly, but like "phonetically", and its really just time_lag*4
        #more importantly this will be used to create daily rasters. 
        # I want there to be 2 steps, i.e., expected locations the day of the known location and then 4 per day in between
        # time_lag -1 gives the days inbetween, then add 2 for the known location of the final day. 
        
        if(stay_in_water[i] == FALSE){ #if the shortest path of P1 to P2 crosses land it will do this loop 
          newMu_pts <- try(make_newMuPts(p1 = P1[i], p2 = P2[i],
                                         net = sightNet_clean, n_muPts =(stepz-2)), silent = TRUE )#make list of expected points that go around land on the shortest possible path
          
          if(class(newMu_pts)[1]== "try-error"){
            iteration_j <-  c(iteration_j, j) #which j use to index track ID
            trackID_j <- c(trackID_j, desiredTrack[j]) #should probably just use this 
            iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
            err_msg <- c(err_msg, newMu_pts) #what is the error??
            i <- i+1 #because this is a while loop it will get stuck
            next #skip this, i iteration, other points in the track will still hold value
          }
          
          firstNewPt_to_P1_dist <- geosphere::distm(st_coordinates(P1[i]), st_coordinates(newMu_pts$x[1]))
          lastNewPt_to_P1_dist <- geosphere::distm(st_coordinates(P1[i]), st_coordinates(newMu_pts$x[length(newMu_pts$x)]))
          
          # i dont have control over the order of points on the line string I think they are always ordered directionally?
          if(firstNewPt_to_P1_dist > lastNewPt_to_P1_dist){ # so this corrects the direction if needed 
            newMu_pts$x <- rev(newMu_pts$x)
          }
          
          newMu_pts <- rbind(st_as_sf(P1[i]), newMu_pts, st_as_sf(P2[i])) #add the starting and ending point
          alpha_land <- seq(from = 0, to= time_lag[i], length.out=stepz)/time_lag[i] # calculate alpha values for each point
          newPts_withSteps <- cbind(alpha_land, newMu_pts$x) #bind locations and alpha values together to prep for calculating possible variance. sigma squared
          sigma_sq <- get_sigma_sq(intervalDays=time_lag[i], alpha=alpha_land, BMvar = BMvar, 
                                   loc_err1=location_err[i], loc_err2=location_err[i+1]) #calclate the variance based on location error, time/alpha values, and the MOTION VARIANCE bmVar
          
          k<-1
          for(k in 1:length(newMu_pts$x)){ #for each expected location 
           tmpGrid <- watergrid_landNAs
            mu_lonlat <- st_coordinates(newMu_pts$x[k]) #grab the coords from the expected point
            mu_cell <- xyFromCell(watergrid_landNAs, cellFromXY(watergrid_landNAs, mu_lonlat)) #determine the lat lon of target cell
            tmpGrid[cellFromXY(watergrid_landNAs, st_coordinates(newMu_pts$x[k]))] <- 2 #assign target cell a different value in prep for gridDist
            ZTZ <- terra::gridDist(tmpGrid, target=2) #calc gridDist to target cell
            dist_mu_to_cell <- geosphere::distm(mu_cell, mu_lonlat) #calculate distance of mu (expected) location to the cell being used to measure
            #dist_mu_to_cell <- terra::distance(mu_cell, mu_lonlat, lonlat =TRUE) These are the same, but rolling with my first love geosphere
            ZTZ[ZTZ == 0] <- dist_mu_to_cell #replace the 0 value with the distance difference
            # since i dont know the directionality of mu to the cell, its not truly correct 
            # the rest of the distances, but at a 1km raster resolution it is the best possible
            ZTZ <- ZTZ^2 #square it, this is the top half of the fraction in the probability distribution fxn
            #use geospere to calculate the distance, because sf_distance uses a different method
            theta <- (1/sqrt(2*pi*sigma_sq[k]))*exp(-ZTZ/(2*sigma_sq[k]))  #calculate the estimated possible habitat use for each grid cell
            
            if(alpha_land[k] == 0 | alpha_land[k] == 1){#if it is at the beginning or the end only value it at half, 
              theta <- theta * 0.5 #because it will be calculated again for the next link
              #It is rectified for the first and last positions for the daily rasters
            }
            #20240103 update I just added the sqrt ehich is CLEARLY in horne et al, but was not in the source code I used
            int <- c(int, theta)
          }
        }else{ #if P1 to P2 does not cross land it will do this loop 
          g<-0
          for(g in seq(tm, time_lag[i], length.out=stepz)){
            alpha <- g/time_lag[i] # fraction between points 1&2 based on timestep place of total time, first iteration would be 0, this alpha is different than in the variance code
            
            #geodesic distance of P1 to P2, most accurate and geosphere can project the point I need
            dist_1to2 <- geosphere::distm(x=c(lon[i], lat[i]), y=c(lon[i+1], lat[i+1]))  
            
            #expected location ~mu~ depending on time point, as a fraction of the whole time gap 
            mu_lonlat <- geosphere::destPoint(p=c(lon[i], lat[i]), b=geosphere::bearing(p1=c(lon[i], lat[i]), p2=c(lon[i+1], lat[i+1])), d=alpha*dist_1to2)
            
            # ~sigma squared~ alpha*(1-alpha) is a normal distribution
            sigma_sq <- time_lag[i]*alpha*(1-alpha)*BMvar +  
              ((1-alpha)^2)*(location_err[i]^2) + #this is adds to the variance distance, because the start and end points may not be where they really are
              (alpha^2)*(location_err[i+1]^2) 
            mu_cell <- xyFromCell(watergrid_landNAs, cellFromXY(watergrid_landNAs, mu_lonlat)) #determine the lat lon of target cell
            tmpGrid <- watergrid_landNAs #reset the tmpGrid to the original water grid
            tmpGrid[cellFromXY(watergrid_landNAs, mu_lonlat)] <- 2 #set target cell 
            ZTZ <- terra::gridDist(tmpGrid, target=2) #calc gridDist to target cell
            dist_mu_to_cell <- geosphere::distm(mu_cell, mu_lonlat) #calculate distance of mu (expected) location to the cell being used to measure
            ZTZ[ZTZ == 0] <- dist_mu_to_cell #replace the 0 value with the distance difference
            # since i dont know the directionality of mu to the cell, its not truly correct 
            # the rest of the distances, but at a 1km raster resolution it is the best possible
            ZTZ <- ZTZ^2
            theta <- (1/sqrt(2*pi*sigma_sq))*exp(-ZTZ/(2*sigma_sq)) #this is the probability distribution fxn, and is a vector that is the same length as rows in the grid
            if(alpha == 0 | alpha == 1){
              theta <- theta * 0.5 #if it is at the beginning or the end only value it at half, 
              #because it will be calculated again for the next link
              #It is rectified for the first and last positions for the daily rasters
            }
            #20240103 update I just added the sqrt ehich is CLEARLY in horne et al, but was not in the source code I used
            int <- c(int, theta)
            
          }
        }
      }
    }
    i<-1
    
    #make the list of rasters a stacked raster, that can be divided by the number of steps 
    steps_theta <- try(terra::rast(int), silent = TRUE) #adding the try error clause, 
    #because if there is only a two point track and it doesnt work for getting around land int is null 
    
    if(class(steps_theta)[1]== "try-error"){
      iteration_j <-  c(iteration_j, j) #which j use to index track ID
      trackID_j <- c(trackID_j, desiredTrack[j]) #should probably just use this 
      iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
      err_msg <- c(err_msg, steps_theta) #what is the error??
      j <- j+1 #because this is a while loop it will get stuck
      next #skip this, j iteration, other points in the track will still hold value
    }
    #length(steps_theta)
    #length(int)
    
    ##create the first daily raster using the first two steps 
    sight1Date <- format(lubridate::date(singleTrack$date[1]), format="%Y%m%d") 
    #as.Date(sight1Date, format="%Y%m%d") for when i want to make movies
    
    whaleDays <- sum(steps_theta[[1:4]], na.rm = TRUE) #first four columns for the first four steps
    
    #this makes all the values add up to 1, because 1 day has occurred 
    #and it will be less than 1 depending on how much the distributions would be on land 
    bb_daily  <- whaleDays/sum(values(whaleDays), na.rm = TRUE) #*0.5 #but I think this weights the steps too high to be left at , and should be halved
    #plot(bb_daily)
    #points(singleTrack$lon, singleTrack$lat)
    
    #save file as 
    terra::writeCDF(bb_daily, filename=paste0(mypath2daily, paste(sight1Date, desiredTrack[j], "bb.nc", sep = "_")), varname="whaleDays", overwrite=TRUE)
    
    if(T_Total != 1){
      h<-1
      for(h in 1:((dim(steps_theta)[3]/8)-1)){ ## The third column is the raster layers, which should be divisible by 8
        # since it was 8 steps a day. Then minus 1 to take off the first two steps and last two steps of the track
        bb_date <- format(base::as.Date(sight1Date, format="%Y%m%d") + h, format="%Y%m%d") 
        
        whaleDays <- sum(steps_theta[[(5+(8*(h-1))):(5+(8*h))]], na.rm = TRUE) #8 consecutive columns for each day that passed, starting with the 3rd column
        bb_daily <- whaleDays/sum(values(whaleDays), na.rm = TRUE)
        
        terra::writeCDF(bb_daily, filename=paste0(mypath2daily, paste(bb_date, desiredTrack[j], "bb.nc", sep = "_")), varname="whaleDays", overwrite=TRUE)
      }
    }#else{
    
    ##creating daily raster for the last day of the track
    bb_date <- format(lubridate::date(singleTrack$date[length(singleTrack$date)]), format="%Y%m%d") 
    
    #again 3 is the diminsion of the number of layers in the spatraster so the last 4 correspond to the day
    whaleDays <- sum(steps_theta[[(dim(steps_theta)[3]-3):dim(steps_theta)[3]]], na.rm = TRUE) #this is to be the last day with the last four columns and the last four steps
    bb_daily  <- whaleDays/sum(values(whaleDays), na.rm = TRUE)#*0.5 
    
    terra::writeCDF(bb_daily, filename=paste0(mypath2daily, paste(bb_date, desiredTrack[j], "bb.nc", sep = "_")), varname="whaleDays", overwrite=TRUE)
    
    #Scaling probabilities so they sum to 1.0
    #sum all the layers in steps_theta, then divide by the total of all the values
    probability <- sum(steps_theta, na.rm=TRUE)/sum(values(steps_theta), na.rm=TRUE)
    
    #multiplied by the total number of days passed
    bb_rast <- probability*T_Total
    
    
    terra::writeCDF(bb_rast, filename=paste0(mypath2fulltracks, paste("bb_rast_", desiredTrack[j], ".nc", sep = "")), varname="whaleDays", overwrite=TRUE)
    
    
  }else{
    
    for(i in 1:(n_locs-1)){
      if(time_lag[i] <= max_lag){
        theta <- NULL
        tm <- 0
        #previously was using 25, which the only thing that im worried about is whales that travel fast, my have clumpy expected distributions, and it not smooth
        #stepz <- 25 #this is increments of the whole time lag aka interval days so that would be 25 steps for 1 day or 25 steps for 25 days (1 a day)
        stepz <- (2+(time_lag[i]-1)*4+2) # this is written ugly, but like "phonetically", and its really just time_lag*4
        #more importantly this will be used to create daily rasters. 
        # I want there to be 2 steps, i.e., expected locations the day of the known location and then 4 per day in between
        # time_lag -1 gives the days inbetween, then add 2 for the known location of the final day. 
        
        if(stay_in_water[i] == FALSE){ #if the shortest path of P1 to P2 crosses land it will do this loop 
          newMu_pts <- try(make_newMuPts(p1 = P1[i], p2 = P2[i],
                                         net = sightNet_clean, n_muPts =(stepz-2)), silent = TRUE )#make list of expected points that go around land on the shortest possible path
          
          if(class(newMu_pts)[1]== "try-error"){
            iteration_j <-  c(iteration_j, j) #which j use to index track ID
            trackID_j <- c(trackID_j, desiredTrack[j]) #should probably just use this 
            iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
            err_msg <- c(err_msg, newMu_pts) #what is the error??
            i <- i+1 #because this is a while loop it will get stuck
            next #skip this, i iteration, other points in the track will still hold value
          }
          
          firstNewPt_to_P1_dist <- geosphere::distm(st_coordinates(P1[i]), st_coordinates(newMu_pts$x[1]))
          lastNewPt_to_P1_dist <- geosphere::distm(st_coordinates(P1[i]), st_coordinates(newMu_pts$x[length(newMu_pts$x)]))
          
          # i dont have control over the order of points on the line string I think they are always ordered directionally?
          if(firstNewPt_to_P1_dist > lastNewPt_to_P1_dist){ # so this corrects the direction if needed 
            newMu_pts$x <- rev(newMu_pts$x)
          }
          
          newMu_pts <- rbind(st_as_sf(P1[i]), newMu_pts, st_as_sf(P2[i])) #add the starting and ending point
          alpha_land <- seq(from = 0, to= time_lag[i], length.out=stepz)/time_lag[i] # calculate alpha values for each point
          newPts_withSteps <- cbind(alpha_land, newMu_pts$x) #bind locations and alpha values together to prep for calculating possible variance. sigma squared
          sigma_sq <- get_sigma_sq(intervalDays=time_lag[i], alpha=alpha_land, BMvar = BMvar, 
                                   loc_err1=location_err[i], loc_err2=location_err[i+1]) #calclate the variance based on location error, time/alpha values, and the MOTION VARIANCE bmVar
          
          k<-1
          for(k in 1:length(newMu_pts$x)){ #for each expected location 
            #ZTZ <- (geosphere::distm(x=cbind(area_grid$x, area_grid$y), y=st_coordinates(newMu_pts$x[k])))^2 #calculate the distance of each grid cell from the expected locations
            mu_lonlat <- st_coordinates(newMu_pts$x[k])
            mu_cell <- xyFromCell(watergrid_landNAs, cellFromXY(watergrid_landNAs, mu_lonlat)) #determine the lat lon of target cell
            tmpGrid <- watergrid_landNAs
            tmpGrid[cellFromXY(watergrid_landNAs, st_coordinates(newMu_pts$x[k]))] <- 2
            ZTZ <- terra::gridDist(tmpGrid, target=2) #calc gridDist to target cell
            dist_mu_to_cell <- geosphere::distm(mu_cell, mu_lonlat) #calculate distance of mu (expected) location to the cell being used to measure
            ZTZ[ZTZ == 0] <- dist_mu_to_cell #replace the 0 value with the distance difference
            # since i dont know the directionality of mu to the cell, its not truly correct 
            # the rest of the distances, but at a 1km raster resolution it is the best possible
            ZTZ <- ZTZ^2 #square it, this is the top half of the fraction in the probability distribution fxn
            #use geospere to calculate the distance, because sf_distance uses a different method, and after squaring the difference between the two methods is too big
            theta <- (1/sqrt(2*pi*sigma_sq[k]))*exp(-ZTZ/(2*sigma_sq[k]))  #calculate the estimated possible habitat use for each grid cell
            
            if(alpha_land[k] == 0 | alpha_land[k] == 1){#if it is at the beginning or the end only vaule it at half, 
              theta <- theta * 0.5 #because it will be calculated again for the next link
              #It is rectified for the first and last positions for the daily rasters
            }
            #20240103 update I just added the sqrt ehich is CLEARLY in horne et al, but was not in the source code I used
            int <- c(int, theta)

          }
        }else{ #if P1 to P2 does not cross land it will do this loop 
          g<-0
          for(g in seq(tm, time_lag[i], length.out=stepz)){
            alpha <- g/time_lag[i] # fraction between points 1&2 based on timestep place of total time, first iteration would be 0, this alpha is different than in the variance code
            
            #geodesic distance of P1 to P2, most accurate and geosphere can project the point I need
            dist_1to2 <- geosphere::distm(x=c(lon[i], lat[i]), y=c(lon[i+1], lat[i+1]))  
            
            #expected location ~mu~ depending on time point, as a fraction of the whole time gap 
            mu_lonlat <- geosphere::destPoint(p=c(lon[i], lat[i]), b=geosphere::bearing(p1=c(lon[i], lat[i]), p2=c(lon[i+1], lat[i+1])), d=alpha*dist_1to2)
            
            # ~sigma squared~ alpha*(1-alpha) is a normal distribution
            sigma_sq <- time_lag[i]*alpha*(1-alpha)*BMvar +  
              ((1-alpha)^2)*(location_err[i]^2) + #this is adds to the variance distance, because the start and end points may not be where they really are
              (alpha^2)*(location_err[i+1]^2) 
            
            tmpGrid <- watergrid_landNAs
            tmpGrid[cellFromXY(watergrid_landNAs, mu_lonlat)] <- 2 #set target cell 
            mu_cell <- xyFromCell(watergrid_landNAs, cellFromXY(watergrid_landNAs, mu_lonlat)) #determine the lat lon of target cell
            ZTZ <- terra::gridDist(tmpGrid, target=2) #calc gridDist to target cell
            dist_mu_to_cell <- geosphere::distm(mu_cell, mu_lonlat) #calculate distance of mu (expected) location to the cell being used to measure
            ZTZ[ZTZ == 0] <- dist_mu_to_cell #replace the 0 value with the distance difference
            # since i dont know the directionality of mu to the cell, its not truly correct 
            # the rest of the distances, but at a 1km raster resolution it is the best possible
            ZTZ <- ZTZ^2
            #ZTZ <- (geosphere::distm(x=cbind(area_grid$x, area_grid$y), y=mu_lonlat))^2
            theta <- (1/sqrt(2*pi*sigma_sq))*exp(-ZTZ/(2*sigma_sq)) #this is the probability distribution fxn, and is a vector that is the same length as rows in the grid
            if(alpha == 0 | alpha == 1){
              theta <- theta * 0.5 #if it is at the beginning or the end only vaule it at half, 
              #because it will be calculated again for the next link
              #It is rectified for the first and last positions for the daily rasters
            }
            #20240103 update I just added the sqrt ehich is CLEARLY in horne et al, but was not in the source code I used
            int <- c(int, theta)
          }
        }
      }
    }
    i<-1
    
    #make the list of rasters a stacked raster, that can be divided by the number of steps 
    steps_theta <- try(terra::rast(int), silent = TRUE) #adding the try error clause, 
    #because if there is only a two point track and it doesnt work for getting around land int is null 
    
    if(class(steps_theta)[1]== "try-error"){
      iteration_j <-  c(iteration_j, j) #which j use to index track ID
      trackID_j <- c(trackID_j, desiredTrack[j]) #should probably just use this 
      iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
      err_msg <- c(err_msg, steps_theta) #what is the error??
      j <- j+1 #because this is a while loop it will get stuck
      next #skip this, j iteration, other points in the track will still hold value
    }
    #length(steps_theta)
    #length(int)
    
    ##create the first daily raster using the first two steps 
    sight1Date <- format(lubridate::date(singleTrack$date[1]), format="%Y%m%d") 
    #as.Date(sight1Date, format="%Y%m%d") for when i want to make movies
    
    whaleDays <- sum(steps_theta[[1:2]], na.rm = TRUE) #first two columns for the first two steps
    
    #this makes all the values add up to 1, because 1 day has occured 
    #and it will be less than 1 depending on how much the distributions would be on land 
    #that shant be an issue tho with the gridDist functionality and landNAs
    bb_daily  <- whaleDays/sum(values(whaleDays), na.rm = TRUE)
    
    #save file as 
    terra::writeCDF(bb_daily, filename=paste0(mypath2daily, paste(sight1Date, desiredTrack[j], "bb.nc", sep = "_")), varname="whaleDays", overwrite=TRUE)
    
    if(T_Total != 1){
      h<-1
      for(h in 1:((dim(steps_theta)[3]/4)-1)){ ## The third column is the raster layers, which should be divisible by 4
        # since it was 4 steps a day. Then minus 1 to take off the first two steps and last two steps of the track
        bb_date <- format(base::as.Date(sight1Date, format="%Y%m%d") + h, format="%Y%m%d") 
        
        whaleDays <- sum(steps_theta[[(3+(4*(h-1))):(3+(4*h))]], na.rm = TRUE) #4 consecutive columns for each day that passed, starting with the 3rd column
        bb_daily  <- whaleDays/sum(values(whaleDays), na.rm = TRUE)
        
        terra::writeCDF(bb_daily, filename=paste0(mypath2daily, paste(bb_date, desiredTrack[j], "bb.nc", sep = "_")), varname="whaleDays", overwrite=TRUE)
      }
    }#else{
    
    ##creating daily raster for the last day of the track
    bb_date <- format(lubridate::date(singleTrack$date[length(singleTrack$date)]), format="%Y%m%d") 
    
    #again 3 is the diminsion of the number of layers in the spatraster so the last 2 correspond to the day
    whaleDays <- sum(steps_theta[[(dim(steps_theta)[3]-1):dim(steps_theta)[3]]], na.rm = TRUE) #this is to be the last day with the last two columns and the last two steps
    bb_daily  <- whaleDays/sum(values(whaleDays), na.rm = TRUE)
    
    terra::writeCDF(bb_daily, filename=paste0(mypath2daily, paste(bb_date, desiredTrack[j], "bb.nc", sep = "_")), varname="whaleDays", overwrite=TRUE)
    
    #Scaling probabilities so they sum to 1.0
    #sum all the layers in steps_theta, then divide by the total of all the values
    probability <- sum(steps_theta, na.rm=TRUE)/sum(values(steps_theta), na.rm=TRUE)
    
    #multiplied by the total number of days passed
    bb_rast <- probability*T_Total
    
    
    terra::writeCDF(bb_rast, filename=paste0(mypath2fulltracks, paste("bb_rast_", desiredTrack[j], ".nc", sep = "")), varname="whaleDays", overwrite=TRUE)
  }
  
}  

err_df <- cbind(iteration_j, trackID_j, iteration_i, err_msg)
write.csv(err_df, paste0("data/processedData/", todayIS, "_", 
                         desiredTrackStart, "to", 
                         desiredTrackEnd,"_",
                         "errorDataFrame_bbmm.csv"))


#7 mins to run the one around land in 0.05 res
# 11 seconds to run around land 0.5 res
