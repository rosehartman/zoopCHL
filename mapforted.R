#Map for Ted

library(deltamapr)
library(sf)
library(tidyverse)


#if you are cool with the EMP station file already in deltamapr, 
#all you have to do is:
EMPstas = filter(P_Stations, Source == "EMP")

#then pick your region file

#There are a bunch in deltamapr, if you want one of those,
#just make sure it's in the same CRS as your stations file
EMPregions = R_EDSM_Regions_1617P1 %>%
  st_transform(crs = st_crs(EMPstas)) #transform the coordinate reference system

#then merge the two using a spatial join
EMPstas2 = st_join(EMPstas, EMPregions)


#plot it

ggplot()+
  geom_sf(data = WW_Delta)+
  geom_sf(data = EMPstas2, aes(color = Region))


#if you are dealing with a .csv file of latitudes and longitudes,
#you first need to turn it into a spatial object

#thsese are the benthic stations from EDI
Benthic = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1036.4&entityid=4e6948186ad756dc2b6de4de41b601f3")

#we need to tell it what the coordinate columns are, and the coordinate reference system (WGS 86, wich is 4326)
Benthicsf = st_as_sf(Benthic, coords = c("Longitude", "Latitude"), crs = 4326) 

#If you have a shapefile somewhere else, you can import it using 'st_read'

#If not, use one of hte deltamapr ones
Benthic2 = st_transform(Benthicsf, crs = st_crs(EMPregions)) %>%
  st_join(EMPregions)


#plot it

ggplot()+
  geom_sf(data = WW_Delta)+
  geom_sf(data = Benthic2, aes(color = Region))
