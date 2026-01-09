## Brownian Bridge Motion Variance calculations to be ran on a supercomputer/HPC
#the slowness is due to the network integration, otherwise it could likely be ran 
#a personal computer

# The Brownian bridge code was adapted from source code in the R packages 
# move (v4.1.6; Kranstauber et al. 2021) and BBMM (v3.0; Nielsen et al. 2013).
# The network was inspired by inspired by the pathroutr package (v0.1.1-beta, London 2020)

#########################################
####################
library(tidyverse)
library(mapdata)  # map_data worldHires coastline
library(raster)
library(sf)
library(sfnetworks)
library('ncdf4')

##############################
##### date labeling prep #####
##############################
todayIS <- format(Sys.Date(), format="%Y%m%d") #this is also in the functions page / loading data because it is critical to me 


##############
## LAND HO ###
###############
#load in the land if the earth is round but I need it to be flat 
mapBase <-st_as_sf(maps::map("worldHires", fill=TRUE, plot=F))
mapBase<- st_cast(mapBase, "MULTIPOLYGON") %>% st_cast("POLYGON")
land_all <- sf::st_transform(mapBase, crs=4326)

#Cropping to north atlantic study area (open to changes on bounds)
sf_use_s2(FALSE) 

## Trying to make the land bigger, 
## but the crop function doesnt work because of some self intersection issue
## now I am looking to add a 0 width buffer 
land <- st_buffer(land_all,dist=0)
lat_bounds= c(20,75) #beautiful map but need to change the lat/lon bounds what should those be 25 min
lon_bounds= c(-110,30)
box <- c(xmin= lon_bounds[1], xmax= lon_bounds[2], ymin= lat_bounds[1], ymax= lat_bounds[2])
land_crop <- st_crop(land, st_bbox(box))
#plot(land_crop) #love the auto colors 
land <- st_union(land_crop) # one land mass without the country boundaries 
land <- st_union(land) #now a cropped simple polygon of points 

sf::sf_use_s2(TRUE)

#st_write(land, "data/rw_polygons/20240324_landCroppedWholeRWRange.shp")

######################################################
#####            read in network                  ####
######################################################
#read in network
sightNet_clean <- st_read('data/processedData/20240617_sightNet_clean.nc')
head(sightNet_clean)
net <- sfnetworks::as_sfnetwork(sightNet_clean, directed=FALSE)
head(net)
sightNet_clean <- net %>%
  activate("edges") %>%
  mutate(weight = edge_length()) #%>% #weighting the edges as much as their length now theres a column for them 
#activate("nodes") %>%
#mutate(bc = centrality_betweenness(weights = weight, directed = FALSE))  ## this is the code from roxel demo, but I will only be doing the edges

make_lineAroundLand<- function(p1, p2, net = sightNet_clean){
  ## Associate the points that have the intersection with the network by finding the nearest node
  fromPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(net,"nodes")),
    sf::st_coordinates(p1),k=1)$nn.idx[,1]
  
  toPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(net,"nodes")),
    sf::st_coordinates(p2),k=1)$nn.idx[,1]
  
  if(fromPoint == toPoint){
    stop(paste0('trackID ',c(unique(singleTrack$trackID), ' nearest network neighbor point is identical')))
  }else{
    #build the shortest route 
    selected_edges <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(net,"edges"),
                                                    from = fromPoint,
                                                    to = toPoint,
                                                    weights = NULL,
                                                    output = "epath")$epath)
    #make a tibble from it 
    edge_geom <- lapply(selected_edges, function(e) {
      g <- net %>% as_tibble(spatial = FALSE)
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


# get nodes for the motion variance 
getNodesAroundLand <- function(p1, p2, net = sightNet_clean){
  ## Associate the points that have the intersection with the network by finding the nearest node
  fromPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(net,"nodes")),
    sf::st_coordinates(p1),k=1)$nn.idx[,1]
  
  toPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(net,"nodes")),
    sf::st_coordinates(p2),k=1)$nn.idx[,1]
  
  if(fromPoint == toPoint){
    stop(paste0('trackID ',c(unique(singleTrack$trackID), ' nearest network neighbor point is identical')))
  }else{
    #build the shortest route 
    selected_nodes <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(net,"nodes"),
                                                    from = fromPoint,
                                                    to = toPoint,
                                                    weights = NULL,
                                                    output = "both")$vpath)
    justNodes <- net %>%
      activate("nodes") %>%
      st_as_sf()
    
    #make a tibble from it 
    node_geom <- lapply(selected_nodes, function(v) {
      g <-   justNodes %>% as_tibble(spatial = FALSE)
      g[v,"geometry"]
    })
    
    
    #make the tibble a df, then an sf df
    df_nodes <- do.call(rbind, node_geom)
    df_nodes_sf <- sf::st_as_sf(df_nodes, crs = 4326)
    
    ## okay amazing now we have the nodes
    
    return(df_nodes_sf)
  }
}


# return points around land

make_newMuPts <- function(p1, p2, net = net, n_muPts = stepz){
  ## Associate the points that have the intersection with the network by finding the nearest node
  fromPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(net,"nodes")),
    sf::st_coordinates(p1),k=1)$nn.idx[,1]
  
  toPoint <-   nabor::knn(
    sf::st_coordinates(sfnetworks::activate(net,"nodes")),
    sf::st_coordinates(p2),k=1)$nn.idx[,1]
  
  if(fromPoint == toPoint){
    stop(c(unique(singleTrack$trackID), 'nearest network neighbor point is identical'))
  }else{
    #build the shortest route 
    selected_edges <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(net,"edges"),
                                                    from = fromPoint,
                                                    to = toPoint,
                                                    weights = NULL,
                                                    output = "epath")$epath)
    #make a tibble from it 
    edge_geom <- lapply(selected_edges, function(e) {
      g <- net %>% as_tibble(spatial = FALSE)
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


#######
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

#######################
## movement variance ##
#####################################################################################
#read in data with at least 3 locations 
mvmnt_n3 <- read.csv(file="data/20250203_mvmnt_with_n_greaterthan2.csv", header=TRUE)
desiredTracks <- unique(mvmnt_n3$trackID)
length(desiredTracks)
str(desiredTracks)
trackpath <- "data/processedData/2025.02.04_indvTrack_csv/"
indvTrack_BMvar<- data.frame()
#Creating NULL vectors to store data

#select network (loaded in the fxn doc)
net <- sightNet_clean

#liklihood function from Horne et al 2007
likelihood <- function(var){    
  v <- t_intervals*alpha*(1-alpha)*var + ((1-alpha)^2)*(loc_err1^2) + 
    (alpha^2)*(loc_err2^2)
  l <- (1/(2*pi*v))*exp(-dist_diff_sq/(2*v))
  -sum(log(l), na.rm=TRUE)
}


#############################
### FOR THE WHOLE SPECIES ###
#############################
#should parallel compute now the network is in place it is soooo slow



t_intervals <- alpha <- dist_diff_sq <- ob <- loc_err1 <- loc_err2 <- trackID <-NULL
iteration_j <-  iteration_i <- err_msg <- trackID_j <- NULL






#this should also be run in parallel to make the most of the super computer 

for(j in c(1:length(desiredTracks))){ 
  singleTrack <- read.csv(paste0(trackpath, paste('sightingTrack', desiredTracks[j], sep = "_")))
  #singleTrack <- read.csv(paste0(trackpath, paste('sightingTrack', 10897, sep = "_")))
  
  singleTrack_sf <- sf::st_as_sf(singleTrack, coords = c('lon','lat'), crs=4326)
  
  n_locs = dim(singleTrack)[1] # number of rows 
  time_lag = singleTrack$daysSincePrev # the first sighting has a 0, and the days are counted to the sighting previous
  singleTrack$location_err <- 2500 ## chose this value to incorporate the movement an individual makes in a day that it is seen. More justification in manuscript and supplement. 
  location_err = singleTrack$location_err 
  lon = singleTrack$lon
  lat = singleTrack$lat
  max_lag = 30 #days no lags greater than 30 days, this is already filtered out, but remains because I dont want to break the code
  
  i <- 2
  while(i < n_locs){
    trackID <- c(trackID, desiredTracks[j])
    ob <- c(ob, i) #ob will be a vector of i in a list c(2,4,6,8,....) 
    t <- time_lag[i]+time_lag[i+1] #lil t is the time total on interval from points  2&3+1&2 = so time passed from 1-3
    t_intervals <- c(t_intervals, t) #will be a vector of c(t1+t2, t2+t3, t3+t4, t4+t5....)
    a <- time_lag[i]/ t#a is the time lag or interval of the 2&3 point/total time from 1-3 
    alpha <- c(alpha, a) #that ratio become alpha? so for my one dayers that'll be 0.5, 0.5, 0.5
    loc_err1 <- c(loc_err1, location_err[i-1])
    loc_err2 <- c(loc_err2, location_err[i+1]) 
   
     ##OKAY HERE IS WHERE I NEED TO ADD THE NETWORK IN 
    sightPoints <- singleTrack_sf$geometry
    P1 <- sightPoints[i-1]
    P2 <- sightPoints[i]
    P3 <- sightPoints[i+1]
    ## FIRST: DOES THE 1 TO 3 DIST CROSS LAND Y/N
    stay_in_water <- get_lineSegs(P1,P3) #make sure spherical geom is back on 
    
    
    if(stay_in_water == FALSE){ #determine if the line drawn from P1 to P3 crosses land or not
      
      nodes_sf <- try(getNodesAroundLand(p1=P1, p2=P3, net = net), silent = TRUE)
      
      
      if(class(nodes_sf)[1]== "try-error"){
        iteration_j <-  c(iteration_j, j) #which j use to index track ID
        trackID_j <- c(trackID_j, desiredTracks[j]) #should probably just use this 
        iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
        err_msg <- c(err_msg, nodes_sf) #what is the error??
        i <- i+1 #because this is a while loop it will get stuck
        tmp_dist_diff_sq <- NA #nas are removed in the optimization and this is iportant for keeping all the vecotrs the same length
        dist_diff_sq <- c(dist_diff_sq, tmp_dist_diff_sq)
        next #skip this, i iteration, other points in the track will still hold value
      }
      ###finding the points that surround the mu point by comparing cumulative distance to expected distance
      P1 <- st_sf(P1, crs=4326)
      st_geometry(P1)<- "geometry"
      P3 <- st_sf(P3, crs=4326)
      st_geometry(P3)<- "geometry"
      
      nodes_sf <- rbind(P1, nodes_sf, P3)
      

      
      
      #use cumsum of distance to locate between what points to geosphere dest point 
      cumDistance <- cumsum(c(terra::distance(terra::vect(nodes_sf), sequential=TRUE))) #the first one can be ignored then 2 will be the distance from 1 to 2 
      totalDistance <-cumDistance[length(cumDistance)] 
      expectedDistance <- a*totalDistance
      p01_Distance <- max(cumDistance[cumDistance<=expectedDistance])
      
      #finding the node BEFORE u (the expected location)
      #the distance closest to the expected distance that is less than the expected distance
      p01 <- nodes_sf$geometry[cumDistance == max(cumDistance[cumDistance<=expectedDistance])][1]
      
      
      #finding the node AFTER u (the expected location)
      #the distance closest to the expected distance that is more than the expected distance
      p02 <- nodes_sf$geometry[cumDistance == min(cumDistance[cumDistance>=expectedDistance])][1]
      
      u <- geosphere::destPoint(p=st_coordinates(p01), b=geosphere::bearing(p1=st_coordinates(p01), p2=st_coordinates(p02)), d=(expectedDistance-p01_Distance))
      u <- sf::st_as_sf(as.data.frame(u), coords = c('lon','lat'), crs=4326)
      

    }else{ #both expected points `u` are sf 
      dist_1to3 <- geosphere::distm(x=c(lon[i-1], lat[i-1]), y=c(lon[i+1], lat[i+1])) #geodesic distance of S1 to S3, most accurate and geosphere can project the point I need
      u <- geosphere::destPoint(p=c(lon[i-1], lat[i-1]), b=geosphere::bearing(p1=c(lon[i-1], lat[i-1]), p2=c(lon[i+1], lat[i+1])), d=a*dist_1to3)
      u <- sf::st_as_sf(as.data.frame(u), coords = c('lon','lat'), crs=4326)
    }
    
    ## THEN: does the expected location to the actual cross land or stay in the water? 
    stay_in_water <- get_lineSegs(u, P2) 
    
    if(stay_in_water == FALSE){
      
      #u needs to be the toPoint and P2 needs to be the from point 
      nodes_sf <- try(getNodesAroundLand(p1=P2, p2=u, net = net), silent = TRUE)
    
      if(class(nodes_sf)[1]== "try-error"){
        iteration_j <-  c(iteration_j, j) #which j use to index track ID
        trackID_j <- c(trackID_j, desiredTracks[j]) #should probably just use this 
        iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
        err_msg <- c(err_msg, nodes_sf) #what is the error??
        i <- i+1 #because this is a while loop it will get stuck
        tmp_dist_diff_sq <-NA
        dist_diff_sq <- c(dist_diff_sq, tmp_dist_diff_sq)
        next #skip this, i iteration, other points in the track will still hold value
      }
      
      #check to account for any distance difference if the P2 point is not apart of the network. 
      #which should probably be seen for the code up top as well. 
      P2toNode <- sum(terra::distance(terra::vect(P2), terra::vect(nodes_sf$geometry[1])))
      
      #find the first node that doesnt cross land, when drawing a line to u 
      stay_in_water_fromNewNodes <- get_lineSegs(u, nodes_sf$geometry) 
      
      if(any(stay_in_water_fromNewNodes) != TRUE){
        iteration_j <-  c(iteration_j, j) #which j use to index track ID
        trackID_j <- c(trackID_j, desiredTracks[j]) #should probably just use this 
        iteration_i <-  c(iteration_i, i) #which i use to index specific location in chain
        err_msg <- c(err_msg, "mu is on land") #what is the error??
        i <- i+1 #because this is a while loop it will get stuck
        tmp_dist_diff_sq <-NA
        dist_diff_sq <- c(dist_diff_sq, tmp_dist_diff_sq)
        next #skip this, i iteration, other points in the track will still hold value
      }
      #the first node from P2 to u that could calculate a shortest distance
      firstClearNode <- min(which(stay_in_water_fromNewNodes  == TRUE))
      #cumulative distance to that node
      cumDistance <- cumsum(terra::distance(terra::vect(nodes_sf$geometry[1:firstClearNode]), sequential=TRUE)) #the first one can be ignored then 2 will be the distance from 1 to 2 
      #total distance from P2 to the first node, all the nodes to the first node that doesn't cross land to u, 
      #and then the final node to u, the expected location
      totalDistance <- P2toNode + cumDistance[length(cumDistance)] + terra::distance(terra::vect(nodes_sf$geometry[firstClearNode]),terra::vect(u)) 
      tmp_dist_diff_sq <- totalDistance^2 #take the diagonal without the first column, so the distance to each subsequent point, summed for the total distance of the line
      dist_diff_sq <- c(dist_diff_sq, tmp_dist_diff_sq)
      
    }else{
      #if the shortest path between u and P2 does not cross land 
      #find the geodesic distance between the actual location P2 and expected u 
      tmp_dist_diff_sq <-geosphere::distm(st_coordinates(P2), st_coordinates(u))^2
      dist_diff_sq <- c(dist_diff_sq, tmp_dist_diff_sq) #adding squared to the distance vector
    }

    i <- i + 1
  }
}
j

# RECODED Likelihood function for Brownian Motion variance estimation
# this is from horne et al 2007 pg 2357
BMvar <- optimize(likelihood, lower=1, upper=100000000000000000000)$minimum #this is in meters so I made the upper limits way big
indvTrack_BMvar <- rbind(indvTrack_BMvar, c(desiredTracks[j], BMvar))
print("this is the movement variance")
print(BMvar)
print("above")

dim(err_df) #132 locations or pairs where the nearest neighbor is identical 
err_df <- cbind(iteration_j, trackID_j, iteration_i, err_msg)

# trackID<- trackID[-1]
# ob <- recon$ob[-1]
# t_intervals <- t_intervals[-1]
# alpha <- alpha[-1]
# loc_err1 <- loc_err1[-1]
# loc_err2 <- loc_err2[-1]

length(trackID)
length(ob)
length(dist_diff_sq)
length(t_intervals)
length(loc_err1)
var_df<- as.data.frame(cbind(trackID, ob, dist_diff_sq, t_intervals, alpha, loc_err1, loc_err2))

nas_wut <- var_df %>%
  filter(is.na(dist_diff_sq))

zeros_wut <- var_df %>%
  filter(dist_diff_sq == 0)

var_df  <- var_df %>%
  filter(dist_diff_sq != 0, !is.na(dist_diff_sq))

summary(zeros_wut)
dim(zeros_wut) #6 rows or instances

write.csv(err_df, paste0("data/processedData/", todayIS,"_errorDataFrame_movementVariance.csv"))
write.csv(var_df, paste0("data/processedData/", todayIS,"_bmVar_distancesTOBECALCULATED.csv"))
write.csv(nas_wut, paste0("data/processedData/", todayIS,"_bmVar_distances_nas_wut.csv"))
write.csv(zeros_wut, paste0("data/processedData/", todayIS,"_bmVar_distances_zeros_wut.csv"))

