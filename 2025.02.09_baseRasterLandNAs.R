## terra test drive 
library(terra)
install.packages("sfnetworks")
library(tidyverse)
library(mapdata)  # map_data worldHires coastline
library(raster)
library(sf)
library(ggspatial)
library(sfnetworks)
library(igraph)
library(tidygraph)

land <- st_read("data/rw_polygons/20240324_landCroppedWholeRWRange.shp")
#nAmerica_container <- read.csv("data/rw_polygons/newShapes/20240326_containerForFineMeshNet_nAmericaOnly.csv")


#make a 'container' polygon to remove the land points later

container_polygon <- nAmerica_container  %>%
  st_as_sf(coords = c(x="lon", y="lat"), crs = 4326) %>% 
  group_by(name) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 



# I dont use this as my base raster anymore, I use a chl raster from ocean color and delete all the data 
# so I have 4km raster cells, jason has 5000 m cells which he claims 

## pick one of the base raster resolutions below!!! 
    # my motion variance is a little less than 500kmsq so having 500 km x 500 km grid cells is annoying to me 
    # so I am doing both. 

#jason roberts res 
baseWaterRaster <- raster::raster(xmn= -110,
                                  xmx= 30, #this is going all the way to 30 degrees east, I'm going to leave for now knowning there is an iceland animal in there 
                                  ymn= 20,
                                  ymx= 75,
                                  res= 0.04889548, #0.0005  in km 0.05044509 0.04559400 0.03935756 0.03192526 for 25, 35, 45, 55 latitude respectively 
                                  #values=0,
                                  crs= 4326)


baseRaster_df<- as.data.frame(baseWaterRaster, xy=TRUE) #i think the raster package needs to be loaded to run this 
baseRaster_df <- baseRaster_df[,1:2]
head(baseRaster_df)

#sat res

bigBase <- terra::rast('data/baseRaster/AQUA_MODIS.20240101.L3m.DAY.CHL.chlor_a.4km.NRT.nc')
newExtent <- terra::ext(-110, 30, 20, 80)
xmn= 15,
xmx= 35, #this is going all the way to 30 degrees east, I'm going to leave for now knowning there is an iceland animal in there 
ymn= 20,
ymx= 80
baseWaterRaster <- terra::crop(bigBase, newExtent)
plot(baseWaterRaster)

baseRaster_df<- as.data.frame(baseWaterRaster, xy=TRUE) #i think the raster package needs to be loaded to run this 
baseRaster_df <- baseRaster_df[,1:2]
head(baseRaster_df)

## 1 km at the equator res 
baseWaterRaster <- raster::raster(xmn= -110,
                                  xmx= 30, #this is going all the way to 30 degrees east, I'm going to leave for now knowning there is an iceland animal in there 
                                  ymn= 20,
                                  ymx= 75,
                                  res= 0.008983112, #0.0005  in km 0.05044509 0.04559400 0.03935756 0.03192526 for 25, 35, 45, 55 latitude respectively 
                                  #values=0,
                                  crs= 4326)



# baseRaster_df <- read.csv(file="data/processedData/20231002_baseRaster_res0.04AquaModis.csv")
# head(baseRaster_df)
# baseWaterRaster<- terra::rast("data/processedData/20231002_baseRaster_res0.04AquaModis.csv")
# head(baseRaster_df)

basePoints <- sf::st_as_sf(baseRaster_df, coords=c('x','y'), crs=4326)
head(basePoints)
container_polygon <- land

container <- st_intersects(basePoints, container_polygon, sparse = FALSE)
noLand_basePoints <- basePoints[container,]
head(noLand_basePoints)
summary(container)




#this is going to be the function that saves me. I need to use spat rasters now in terra 
terra::gridDist()

#I just read the lore on why things are retiring and its because the guy that maintaing rgdal, rGEOS, etc is retiring 
# and i would like to say good for him.

# Remember, the distance calculated by gridDist is based on the grid resolution. 
# If you want more precise distances, you can use a finer grid resolution, but keep 
# in mind that it will increase computation time and memory usage.

## gotta make a new base raster where all the points are NA if they overlap with land

## so 

#water_landNA <- 


plot(land)
land_t <- vect(land)
water_t <- terra::rast(baseWaterRaster)
values(water_t) <- 0 
landNAs<-  terra::mask(water_t, land_t, inverse=TRUE)
plot(landNAs)

terra::writeCDF(landNAs, "data/baseRaster/20240814_1km_resRasterLandNAs_forGridDist.nc", varname="water", overwrite=TRUE)

watergrid_landNAs<- terra::rast("data/baseRaster/20240619_Sat4kmresRasterLandNAs_forGridDist.nc")
plot(watergrid_landNAs)


##################

# max grid distance between any two random points by district.
library(terra)
#terra 1.6.47

