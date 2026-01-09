##20230712 Clean Version of building a network for sightings to navigate around land 

library(tidyverse)
library(mapdata)  # map_data worldHires coastline
library(raster)
library(sf)
library(ggspatial)
library(sfnetworks)
library(igraph)
library(tidygraph)

#########################################################################
#    Load and check data
#########################################################################

#need to bring in all of bobs sightings 
#bobs sightings include everything!! from everyone?? 
sightings_bob <- read.csv("data/fromNARWC/AbbySightingData_Apr2024_fromBobKenney.CSV")
sightings_phil <- read.csv('data/fromNARWC/Pendleton Kreuser update Right Whale Identification Catalog export- 2024-04-17.csv')
head(sightings_phil)
head(sightings_bob)
table(sightings_bob$TYPE)
unique(sightings_bob$YEAR)
#theres one NA so im just going to remove it
sightings_bob <- sightings_bob[!is.na(sightings_bob$LONG_DD),] 
dim(sightings_bob) #64705 sightings but why is that less than the identified indv data base?

sightings_bob_sf <- sf::st_as_sf(sightings_bob,coords = c('LONG_DD','LAT_DD'), crs=4326)
head(sightings_bob_sf)


#sightings_phil = read.csv('data/rw_sightings_2024-04-17.csv') #named phil because the sourced raw data comes from Philip H. my dawg.
dim(sightings_phil)
sightings_phil <- sightings_phil %>%
  filter(Latitude != 0) %>%
  filter(Longitude != 0)


sightings_phil_sf <- sf::st_as_sf(sightings_phil,coords = c('Longitude','Latitude'), crs=4326)

#This is adding the supplemental point in FL so the whales can navigate 
# around florida for the ones that went into the Gulf of Mexico, and one point in the st Lawernce esturary 
supPoints <- as.data.frame(rbind(c(-80,24.4), c( -67.251198, 49.161496)))
supPoints <- sf::st_as_sf(supPoints,coords = c('V1','V2'), crs=4326)


#join bob's sightins, phils sightings, and the supplemental points 
all_geom_sf <- c(sightings_bob_sf$geometry, sightings_phil_sf$geometry, supPoints$geometry)
# now only keep the points that are uniquw 
all_geom_df <- as.data.frame(
  do.call(rbind, unique(all_geom_sf))
)
length(all_geom_sf) #150033 points 
dim(all_geom_df) # 86956 points 

head(all_geom_df)
colnames(all_geom_df) <- c("lon", "lat")
#save as df 
write.csv(all_geom_df,"data/processedData/networkData/2024.10.28_uniqueSightingLocations_forNetwork.csv", row.names = FALSE)


#######################################
### old code below for removing dead whales, but I don't think i have those deets on bobs 
#####################################
#had to lexical scope because read.csv has been adding "i..." to first column name and then botching code
head(sightings_phil) 
glimpse(sightings_phil)

# #files_yr1978to1989_ID_Q1 <- filter(files %in% c(paste0("bb_rast_", yr1978to1989_ID_Q1$trackID, ".nc")))
# deadWhales <- c('DEAD ON BEACH','FLTG DEAD', 'FRST DEAD', 'STRAND')
# #deadWhales <- c('DEAD', 'dead', 'Dead')
# # Used just sight notes to see if dead being mentioned was leaving out important info, 
# # and I only decided to add the behavior strand, which applies to two sightings 

sightsAlive <- NULL

sightsAlive <-   sightings_phil[grep(paste(deadWhales,collapse="|"), sightings_phil$Behaviors, invert=TRUE),]

locationAccurate <- c('Class Z',  'Class 0','approximated', 'approximate', 'approx') # no class Z, but about 130 class 0 whcich is the radius is over 1500 m accuracy
test <-   sightsAlive[grep(paste(locationAccurate,collapse="|"), sightsAlive$SightingNote, invert=FALSE),]
plot(test)

dim(sightsAlive)
dim(test)

glimpse(sightsAlive)
sightsAlive <- sightsAlive %>%
  filter(Latitude != 0) %>%
  filter(SightingYear > 1979) #sighting effort became heavy in 1980, so using that as a cut off.  


#Im sure there are other extraneous sighting things I need to filter out but im happy with alive whales only for now 
#looking at the pull and shape R script 
dim(sightsAlive)
dim(sightings_phil) #~290 dead or stranded whale sightings 

plot(sightsAlive$Longitude, sightsAlive$Latitude)

sightsAlive<- sf::st_as_sf(sightsAlive,coords = c('Longitude','Latitude'), crs=4326)
#####################################################################################

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

sf_use_s2(TRUE) 
###########################################################
###   Making a network of all of the known sightings   ####
###########################################################
#dist_SightsAlive <- read.csv("data/processedData/allLIVEsightingsDistToCoast_m_km.csv")
#back to SF geom
dist_SightsAlive <- sf::st_as_sf(dist_SightsAlive,coords = c('X','Y'), crs=4326)
sightsAlive <- dist_SightsAlive
#length(dist_SightsAlive$SightingId)

dim(dist_SightsAlive)
glimpse(dist_SightsAlive)
##########################################################3
all_geom_df <- read.csv("data/processedData/networkData/2024.10.28_uniqueSightingLocations_forNetwork.csv")

#remake into an sf dataframe 
all_geom_sf<- sf::st_as_sf(all_geom_df, coords= c('lon', 'lat'), crs=4326)
head(all_geom_sf)
length(all_geom_sf$geometry)

sightsAlive <- all_geom_sf

#######################################################
##### build the network 
##########################################################
## For cleaning these up... Ideally I would like ones in the same exact location to be gone, 
## but I would like to weight them and consolidate based on their hotspot-ness, so do I include all nodes??
## So theoretically I could figure that out later, but I'm gonna just do all the sightings now

tri_sightMap <- filtered2 %>%
  activate("nodes") %>%
  st_as_sf() %>%
  st_union()%>%
  st_triangulate( bOnlyEdges = TRUE) %>% #triangulating long lat data doesnt happen correctly cool
  st_cast('LINESTRING') %>%
  st_sf()

plot(tri_sightMap)

tri_sightMap <- sightsAlive %>%
  st_union()%>%
  st_triangulate( bOnlyEdges = TRUE) %>% #triangulating long lat data doesnt happen correctly cool
  st_cast('LINESTRING') %>%
  st_sf()

plot(tri_sightMap)

sightNet <- sfnetworks::as_sfnetwork(tri_sightMap,
                                     directed = FALSE,
                                     length_as_weight = TRUE)
#started at like 4
filtered = sightNet %>%
  activate("edges") %>%
  filter(!edge_intersects(land)) %>%
  activate("nodes") %>% 
  filter(!tidygraph::node_is_isolated())
 

sightNet_clean2 <- sightNet_clean
sightNet_clean <- sfnetworks::activate(filtered,"edges") %>% sf::st_as_sf()

st_write(sightNet_clean, 'data/processedData/networkData/2024.10.28_sightNet_clean.nc', format="CDF", overwrite = TRUE) 
#old mac at home can't read .nc??? unsure, but it is not supported on the dirvers for the sf packages I have here 
st_write(sightNet_clean, 'data/processedData/networkData/2024.10.28_sightNet_clean.shp')

###############################################
#### plot and check the network 
#############################################
p2 <- ggplot() +
  ggspatial::annotation_spatial(land, fill = "cornsilk3", size = 0) +
  ggspatial::layer_spatial(sightNet_clean, size = 0.5) +
  #geom_sf(data = singleTrack_sf, color='deeppink3', alpha = 0.25, size = 0.5)+
  #geom_line(data= singleTrack, aes(x=lon, y=lat))
  theme_void()

ggsave(p2, file= 'figures/networkNaked.png', height = 4, width = 7)
ggsave(p2, file= 'figures/202306_sightingsAroundLand/landNetworks4.png', height = 4, width = 7)


## OMG ITS GORGEOUS AND BETTER THAN PATHROUTR 
# but its missing around the florida keys

######################################################
#####            read in network                  ####
######################################################
#read in network
#sightNet_clean <- st_read('data/processedData/sightNet_clean.nc')
sightNet_clean <- st_read('data/processedData/networkData/2024.10.28_sightNet_clean.nc')

head(sightNet_clean)
net <- sfnetworks::as_sfnetwork(sightNet_clean, directed=FALSE)
head(net)
sightNet_clean <- net %>%
  activate("edges") %>%
  mutate(weight = edge_length()) #%>% #weighting the edges as much as their length now theres a column for them 
#activate("nodes") %>%
#mutate(bc = centrality_betweenness(weights = weight, directed = FALSE))  ## this is the code from roxel demo, but I will only be doing the edges

lat_crop1= c(20,70) #beautiful map but need to change the lat/lon bounds
lon_crop1= c(-100,-19)

##map of all sightings and the network of all alive sightings 
p5 <- ggplot() +
  geom_sf(data = land, fill="cornsilk3") +
  ggspatial::layer_spatial(sightNet_clean, linewidth = 0.1) +
  geom_sf(data = sightsAlive, color='deeppink3', alpha = 0.25, size = 0.1)+ #this one the data is not here in this script... 
  #ggspatial::layer_spatial(l, color = "darkgrey", size = 0.5) +
  #ggspatial::layer_spatial(dat_sf, color = "deepskyblue3", size = 0.5) +
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(crs = st_crs(4326), 
           xlim = lon_crop1, ylim = lat_crop1,
           expand=FALSE)+
  theme_bw()
#ggsave(p5, file= 'figures/202306_sightingsAroundLand/rwWorldCropSightingsPlusNet.png', height = 4, width = 7)
ggsave(p5, file = paste0( "figures/", todayIS, "linerwWorldNet_w_sightsCrop2.png"), height = 4, width = 7)


##Croping station
##NE with GSL
lat_crop= c(38,52) #beautiful map but need to change the lat/lon bounds
lon_crop= c(-75,-55)

cropmap2 <- ggplot() +
  geom_sf(data = land, fill="cornsilk3") +
  ggspatial::layer_spatial(sightNet_clean, linewidth = 0.1) +
  geom_sf(data = sightsAlive, color='deeppink3', alpha = 0.25, size = 0.5)+ #this one the data is not here in this script... 
  ggspatial::layer_spatial(l, color = "black", size = 0.5) +
  ggspatial::layer_spatial(dat_sf, color =  "black", size = 0.5) +
  ggspatial::layer_spatial(edge_l, color =  "deepskyblue3", linewidth = 0.5)+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(xlim = lon_crop, ylim = lat_crop, expand=FALSE, crs=4326) + # Crop map edges
  #theme(text=element_text(size=20))+
  theme_bw()

ggsave(cropmap2,  file = paste0( "figures/", todayIS, "track4860_solution2_lineSightingsNetCrop_wnet.png"), height = 4, width = 7)




##########################################################################
####  load in test data set in need of navigating around Nova Scotia  ####
##########################################################################
#read in mvmnt data
mvmnt <- read.csv(file= 'data/20240615_mvmnt_with_dist.csv') #has the speeds and distances 

#filtering to get the long recent chain to GSL in the pathroutrWithWhales
places <- unique(mvmnt$movement)
movementGSL<- mvmnt[grepl("GSL", mvmnt$movement),]

trackIDNovaScosh <- unique(movementGSL$trackID)

length(trackIDNovaScosh)

mvmntGSL <- filter(mvmnt, trackID %in% trackIDNovaScosh)
#okay lets look at the new tracks 
recentGSLmvmnt <- filter(mvmntGSL, time_period == "yr2015to2020")

glimpse(recentGSLmvmnt)
max(recentGSLmvmnt$totalChain)
table(recentGSLmvmnt$totalChain)

recentGSLmvmnt[recentGSLmvmnt$totalChain == 16,]

BOFtoGSL <- recentGSLmvmnt[recentGSLmvmnt$movement == "GSC_to_GSL",]

dat <- recentGSLmvmnt[recentGSLmvmnt$trackID == 4283,] #4860

dat <- dat[order(dat$sight1_date),]
dat <- mvmnt[mvmnt$trackID == 4860,]
dat
dat_sf <- st_as_sf(dat, coords = c("sight1_longitude", "sight1_latitude"), 
                   crs = 4326, agr = "constant")

head(dat_sf)
tail(dat_sf)

st_line_sample(l)
l <- dat_sf %>% 
  dplyr::summarise(do_union = FALSE) %>% 
  sf::st_cast('LINESTRING')

##########################################################################
###            creating the points that go around                      ###
##########################################################################


head(dat_sf)
dat_sf[9,]

## Check if track crossed land 
test <- st_intersects(l, land, sparse=FALSE)#so only the first line crosses the land, and it is the only one that needs to be interpolated 
l


## Associate the points that have the intersection with the network by finding the nearest node
fromPoint <-   nabor::knn(
  sf::st_coordinates(sfnetworks::activate(sightNet_clean,"nodes")),
  st_coordinates(dat_sf$geometry[9]),k=1)$nn.idx[,1]

toPoint <-   nabor::knn(
  sf::st_coordinates(sfnetworks::activate(sightNet_clean,"nodes")),
  st_coordinates(dat_sf$geometry[10]),k=1)$nn.idx[,1]


## edge df are from below with the igraphs use not the sf one 

plz <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(sightNet_clean,"edges"),
                                     from = fromPoint,
                                     to = toPoint,
                                     weights = NULL,
                                     output = "epath")$epath)

edge_geom <- lapply(plz, function(e) {
  g <- sightNet_clean %>% as_tibble(spatial = FALSE)
  g[e,"geometry"]
})


##########################
##This is for nodes but I dont think thats needed - v is for verticies 
plz <- unlist(igraph::shortest_paths(graph = sfnetworks::activate(sightNet_clean,"nodes"),
                                     from = fromPoint,
                                     to = toPoint,
                                     weights = NULL,
                                     output = "vpath")$vpath)

node_geom <- lapply(plz, function(e) {
  g <- sightNet_clean %>% as_tibble(spatial = FALSE)
  g[e,"geometry"]
})
plot(unlist(node_geom))
#####################
node_geom
edge_geom

df_edges <- do.call(rbind, edge_geom)
df_edges_sf <- st_as_sf(df_edges, crs = 4326 )

edge_l <- df_edges_sf %>% 
  dplyr::summarise(do_union = FALSE) %>% 
  sf::st_cast('LINESTRING')

edge_l <- edge_l %>% 
  st_sf() %>% 
  st_combine() %>% 
  st_sf() %>% 
  st_line_merge() %>% 
  st_cast("LINESTRING") 


plot(node_geom)


df_nodes <- do.call(rbind, node_geom)
df_nodes_sf <- st_as_sf(df_nodes, crs = 4326 )

node_l <- df_nodes_sf %>% 
  dplyr::summarise(do_union = FALSE) %>% 
  sf::st_cast('POINT')

plot(edge_l)




#might not be needed at all, except to create vairables for how many pts 
st_geod_length(edge_l)
segLength <- st_geod_length(edge_l)/(22*24) #22 days and 24 hours in a day need to code in how many slots a day I have in the other code
my_segments <- st_segmentize(edge_l, units::set_units(10, km))
track_pts <- st_cast(my_segments, "POINT")

st_cast(pt_ls, "POINT")
pt_ls <- st_line_sample(st_transform(edge_l, 3857), n = 22*24)
plot(st_transform(pt_ls, 4326))

#convert to this other projection first then go back 



# Checking the distance between the nearest features on these bad boys
WINNER <- st_transform(pt_ls, 4326) %>%
  st_cast('POINT') %>%
  st_as_sf() 
naborpts <- WINNER[st_nearest_feature(WINNER),]
for(i in 1:nrow(WINNER)){
  naborpts$x[i] <- WINNER$x[i+1]
}
naborpts
distNabor <- data.frame(dist = rep(0, nrow(naborpts)))
for(i in 1:sum(nrow(WINNER)-1)){
  distNabor$dist[i] <- st_distance(WINNER$x[i], naborpts$x[i])
}
distNabor

distNabor <- cbind(distNabor, WINNER$x)

new_pts_summary <- summary(distNabor[-528,])

ggplot() + 
  ggspatial::annotation_spatial(land, fill = "cornsilk3", size = 0) +
  #ggspatial::layer_spatial(distNabor[-528,], size = 0.5, aes(color = dist)) +
  ggspatial::layer_spatial(l, color = "red", size = 0.25) +
  ggspatial::layer_spatial(dat_sf, color = "red", size = 0.5) +
  #ggspatial::layer_spatial(segs, color = "red", size = 0.5) +
  theme_bw()

###################
## I used igraph ##
###################
#then use the paths fxn to connect them through the network. 
### sf Networks fxn paths below was nice for illustrating, and it recognizes both nodes and edges
## igraph is set as edge or verticies graph 
paths = st_network_paths(sightNet_clean, from = fromPoint, to = toPoint, weights = "weight")
paths

nodesRoute <- paths %>%
  slice(1) %>%
  pull(node_paths) %>%
  unlist()

paths %>%
  slice(1) %>%
  pull(edge_paths) %>%
  unlist()

plot_path = function(node_path) {
  net %>%
    activate("nodes") %>%
    slice(node_path) %>%
    plot(cex = 1.5, lwd = 1.5, add = TRUE)
}

colors = sf.colors(3, categorical = TRUE)

plot(net, col = "grey")

paths %>%
  pull(node_paths) %>%
  walk(plot_path)

lines(l)
net %>%
  activate("nodes") %>%
  st_as_sf() %>%
  slice(c(fromPoint,toPoint)) %>%
  plot(col = colors, pch = 8, cex = 2, lwd = 2, add = TRUE)

class(paths)
paths


