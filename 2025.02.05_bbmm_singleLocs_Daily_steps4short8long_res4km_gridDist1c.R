# single location BBMM output 
#referred to as sightings not associated with a track in manuscript:

#Sightings not associated with a track
#Sightings where individuals were not seen within 30 days prior or after 
#were not associated with a track. However, these sightings were incorporated 
#into the model with a single day occurrence distribution (n = 9,422) to ensure 
#that all available right whale identification data were included in occurrence distributions 
#for ecological interpretation. The Brownian bridge was applied to a false two-location track created 
#with the sighting as an initial and terminal location. The maximum time interval (30 days) was used, 
#but only the first daily occurrence distribution was kept and associated with the date of the sighting. 
#This created a probability distribution of undirected movement for one day with the maximum uncertainty in 
#time allowed by the model (Figure S2). 


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


#these may just need to 




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

## using terra in certain functions is not parallelable, 
## so instead I often break up jobs into chunks to run quicker. 

# prep labeling
todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 


#load in base raster
## spatrasters of the water points. 
landNAs <- terra::rast("data/baseRaster/20240619_Sat4kmresRasterLandNAs_forGridDist.nc")
newExtent <- terra::ext(-100, -40, 20, 70)
watergrid_landNAs <- terra::crop(landNAs, newExtent)

######################
## Variance plug in ##
######################

#set the motion variance to be using for all tracks. 
BMvarALL <- 464099334 ## THIS ONE, this is the motion var with  2.5 km location err

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
mypath2daily <- paste0("data/processedData/results/", todayIS, "bbmmRasterDataDaily_res4km_bmvar", BMvarALL_sqkm, "sqkm_locerr2.5k_singleLocs/")
dir.create(mypath2daily)

mypath2fulltracks <- paste0("data/processedData/results/", todayIS, "bbmmRasterDataFullTracks_res4km_bmvar", BMvarALL_sqkm, "sqkm_locerr2.5k_singleLocs/")
dir.create(mypath2fulltracks)


iteration_j <-  iteration_i <- err_msg <- trackID_j <- NULL


#######################
### whale track ids ###
#######################
dat <- read.csv("data/processedData/2025.02.05_singleLocs_NOTinModel_fakesightpairs_with_jdate.csv")
desiredTrack <- unique(dat$trackID) #prepping all individual tracks to be indexed in the loop by their track ID
length(desiredTrack)

## when running all tracks break up into sections 
# desiredTrack <- desiredTrack[615:length(desiredTrack)]
# length(desiredTrack)
 breaksDesired <- round(seq(from = 1, to = length(desiredTrack), length.out=5))
running_set <- 4
desiredTrackStart <- breaksDesired[running_set] +1
desiredTrackEnd <- breaksDesired[running_set+1]
# print(c("running set #", running_set, ":", desiredTrackStart, "to",
#       desiredTrackEnd ))

####################################################
#run it
j <-1
for(j in 1){
  trackpath <- 'data/processedData/2025.02.05_singleLoc_indvTrack_csv/'
  singleTrack <- read.csv(paste0(trackpath, paste('sightingTrack',desiredTrack[j], sep = "_")))
  #then convert to sf, for intersections with and routing around land
  singleTrack_sf <- sf::st_as_sf(singleTrack, coords = c('lon','lat'), crs=4326)
  
  n_locs <- dim(singleTrack)[1] # number of rows or known locations 
  time_lag <- 30 #singleTrack$daysSincePrev[-1] #Takes out the 0 at the top, but will have (n_locs-1) values
  T_Total <- sum(time_lag) #total time in the track
  singleTrack$location_err <- 2500 #adding this location error, because 2.5km is the mean distance whales sighted multiple times in one day moved (measured the full track)
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
  # but probability is for the whole north Atlantic 
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
  
  for(i in 1:(n_locs-1)){
    if(time_lag[i] <= max_lag){
      theta <- NULL
      tm <- 0
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
        for(g in seq(tm, time_lag[i], length.out=stepz)[1:2]){
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
  plot(bb_daily)
  
}




test <- min(values(tenDay_int-thirtyDay_int, na.rm=TRUE) > 0, na.rm=TRUE)

err_df <- cbind(iteration_j, trackID_j, iteration_i, err_msg)
write.csv(err_df, paste0("data/processedData/", todayIS, "_", 
                         desiredTrackStart, "to", 
                         desiredTrackEnd,"_",
                         "errorDataFrame_bbmm.csv"))

