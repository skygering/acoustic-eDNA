library(rgeos) 
library(sf) # Represents simple features as native R objects
library(rnaturalearth) # package that facilitates interaction with Natural Earth map data
library(rnaturalearthdata) # world vector map from Natural Earth
library(raster)
library(rgdal)
library(mapdata) #https://hansenjohnson.org/post/bathymetric-maps-in-r/
library (ggspatial)
library(marmap) 


df <- read_csv("/Volumes/GeringSSD/GU201905_Output/edna_abc.csv")