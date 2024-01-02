#Check out chlorophyll and zooplankton correlations
library(tidyverse)
library(zooper)
library(lubridate)
library(Hmsc)
library(pscl)
library(lme4)
library(lmerTest)
library(visreg)
library(MuMIn)
library(wql)
library(MASS)
#######################################################################################
#what if we group the zoops by region and chlorophyll by region?
#Then run the model on the monthly means

require(discretewq) # water quality data https://github.com/sbashevkin/discretewq
require(deltamapr) # SubRegions dataset https://github.com/InteragencyEcologicalProgram/deltamapr
require(sf)
require(hms)
require(rlang)
require(readr)
require(dtplyr) # To speed things up
require(ggplot2) # Plotting
require(geofacet) # plotting


Regions<-read_csv("Rosies_regions.csv")

## Load Delta Shapefile from Brian
Delta<-deltamapr::R_EDSM_Subregions_Mahardja_FLOAT%>%
  filter(SubRegion%in%unique(Regions$SubRegion))%>%  #Filter to regions of interest
  dplyr::select(SubRegion)

Regs = unique(Regions[,c(1,5)])
Delta = merge(Delta, Regs)

Data<-wq(End_year=2022,
         Sources = c("EMP", "USGS_SFBS", "USGS_CAWSC", "NCRO"))

WQzoops<-function(variable, narm = "FALSE"){
  vardata<-Data%>% # Select all long-term surveys (excluding EDSM and the USBR Sacramento ship channel study)
    filter(!is.na(.data[[variable]]) & !is.na(Latitude) & !is.na(Datetime) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
    mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"), #Convert to a timezone without daylight savings time
           Date = with_tz(Date, tz="America/Phoenix"), # Calculate difference from noon for each data point for later filtering
           Station=paste(Source, Station),
           Time=as_hms(Datetime), # Create variable for time-of-day, not date. 
           Noon_diff=abs(hms(hours=12)-Time))%>% # Calculate difference from noon for each data point for later filtering
    lazy_dt()%>% # Use dtplyr to speed up operations
    group_by(Station, Source, Date)%>%
    filter(Noon_diff==min(Noon_diff))%>% # Select only 1 data point per station and date, choose data closest to noon
    filter(Time==min(Time))%>% # When points are equidistant from noon, select earlier point
    ungroup()%>%
    as_tibble()%>% # End dtplyr operation
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
    st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
    st_join(Delta, join=st_intersects)%>% # Add subregions
    filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
    st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
    lazy_dt()%>% # Use dtplyr to speed up operations
    
    
    group_by(SubRegion)%>%
    mutate(var_sd=sd(.data[[variable]]), 
           var_ext=(.data[[variable]]-mean(.data[[variable]],
                                           na.rm = narm))/var_sd)%>% # Calculate variable SD for each subregion, then the number of SD's each data point exceeds the mean
    ungroup()%>%
    as_tibble()%>% # End dtplyr operation
    filter(var_ext<10)%>% # Filter out any data points that are more than 10 SDs away from the mean of each subregion
    mutate(Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
           Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
                            Month%in%6:8 ~ "Summer",
                            Month%in%9:11 ~ "Fall",
                            Month%in%c(12, 1, 2) ~ "Winter",
                            TRUE ~ NA_character_))%>%
    lazy_dt()%>% # Use dtplyr to speed up operations
    group_by(Month, Season, Region, Year)%>%
    summarise(var_month_mean=mean(.data[[variable]], na.rm = narm), N=n())%>%
    ungroup()%>%
  as_tibble() %>%
  complete(nesting(Month, Season), Region, Year, fill=list(N=0)) # Fill in NAs for variable (and 0 for N) for any missing month, subregion, year combinations to make sure all months are represented in each season
    
  
  variable2<-sym(variable)
  
  out<-vardata%>%
    filter(Year%in%1975:2021)%>%
    group_by(Season, Year, Region, Month)%>%
    summarise({{variable}}:=mean(var_month_mean, na.rm = narm), "N_{{variable2}}":=sum(N), .groups="drop")
  
    
    
   cat(paste("\nFinished", variable, "\n"))
  
  return(out)
}


vars<-c("Temperature", "Secchi", "Chlorophyll", "Conductivity")
yearseasons<-expand_grid(Month= 1:12, 
                         Year=1970:2022, Region = unique(Regions$Region))

WQ_listNARM <-map(vars, WQzoops, narm = TRUE)

WQ_dataNARM<-map(WQ_listNARM, ~full_join(yearseasons, .x, by=c("Month", "Year", "Region"))%>%
                   dplyr::select(contains(vars)))%>%
  bind_cols(yearseasons, .)%>%
  arrange(Region, Year, Month) %>%
  mutate(Chlag = lag(Chlorophyll),
         logChl = log(Chlorophyll),
         logCon = log(Conductivity),
         logChlag = lag(logChl))




###############################


#get all the EMP and 20mm data from zooper
#If I'm using GAMs, I can do winter data after all.
Zoops = Zoopsynther(Data_type = "Community",
                     Sources = c("EMP", "20mm", "DOP", "FMWT", "STN"),
                     Size_class = "Meso",
                   #  Months = c(3:10),
                     Date_range = c("2000-01-01", "2022-12-30"),
                     Redownload_data = F)

test = group_by(Zoops, Year, Taxlifestage, Taxname) %>%
  summarise(N = n(), Mean = mean(CPUE, na.rm = T))


Zoops2 = dplyr::filter(Zoops, Undersampled == FALSE, !is.na(Latitude), Taxname != "Oithona similis", Lifestage != "larva") %>%
                      mutate(Taxname2 = case_when(
  Taxname == "Acanthocyclops vernalis" ~ "Acanthocyclops",
  Taxname == "Sinocalanus_UnID" ~ "Sinocalanus",
  Taxname == "Sinocalanus doerrii" ~ "Sinocalanus",
  Taxname == "Acartiella sinensis" ~ "Acartiella_UnID",
  Taxname == "Eurytemora" ~ "Eurytemora affinis",
  TRUE ~ Taxname
  ), Month = month(Date)) %>%
  group_by(SampleID, Taxname2, Month, Year, Date, Source, Station, Longitude, Latitude, SalSurf, Secchi) %>%
  summarize(CPUE = sum(CPUE))

# ZoopsSum = Zoops2 %>%
#   group_by(Station, Source, Date)%>%
#   st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
#   st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
#   st_join(Delta, join=st_intersects)%>% # Add subregions
#   filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
#   st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
#   mutate(Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
#          Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
#                           Month%in%6:8 ~ "Summer",
#                           Month%in%9:11 ~ "Fall",
#                           Month%in%c(12, 1, 2) ~ "Winter",
#                           TRUE ~ NA_character_))%>%
#   group_by(Month, Season, Region, Year, Taxname2)%>%
#   summarise(var_month_mean=mean(CPUE, na.rm = TRUE), N=n())%>%
#   ungroup()%>%
#   as_tibble() %>%
#   complete(nesting(Month, Season), Region, Year, fill=list(N=0)) # Fill in NAs for variable (and 0 for N) for any missing month, subregion, year combinations to make sure all months are represented in each season


ZoopsSum = Zoops2 %>%
  group_by(Station, Source, Date, SampleID, SalSurf, Secchi)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
  mutate(Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
         Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
                          Month%in%6:8 ~ "Summer",
                          Month%in%9:11 ~ "Fall",
                          Month%in%c(12, 1, 2) ~ "Winter",
                          TRUE ~ NA_character_))
#Let's filter months and salinities first, then join to water quality later.


save(ZoopsSum, WQ_dataNARM, file = "ZoopsSum.RData")

###############################
#Ok, Now the same witht the microzoops


#get all the EMP data from zooper
Zoopsm = Zoopsynther(Data_type = "Community",
                    Sources = c("EMP"),
                    Size_class = "Micro",
                   # Months = c(3:10),
                    Date_range = c("2000-01-01", "2022-12-30"),
                    Redownload_data = F)


Zoops2m = dplyr::filter(Zoopsm, !Undersampled, Order != "Calanoida", !is.na(Latitude)) %>%
  mutate(Taxname2 = case_when(
    Taxname == "Limnoithona sinensis" ~ "Limnoithona_UnID",
    Taxname == "Limnoithona tetraspina" ~ "Limnoithona_UnID",
    Taxname == "Sinocalanus_UnID" ~ "Sinocalanus",
    Taxname == "Synchaeta bicornis" ~ "Synchaeta_UnID",
    TRUE ~ Taxname
  ), Month = month(Date)) %>%
  group_by(SampleID, Taxname2, Month, Year, Date, Source, Station, Longitude, Latitude, SalSurf, Secchi) %>%
  summarize(CPUE = sum(CPUE)) %>%
  mutate(Taxname2 = str_replace(Taxname2, "_UnID", ""))

ZoopsSumm = Zoops2m %>%
  group_by(Station, Source, Date)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
  mutate(Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
         Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
                          Month%in%6:8 ~ "Summer",
                          Month%in%9:11 ~ "Fall",
                          Month%in%c(12, 1, 2) ~ "Winter",
                          TRUE ~ NA_character_))



save(ZoopsSumm, file = "ZoopsSumm.RData")

###############################
#I'm not particularly hopeful on mysids, but we'll give it a go
# 
# 
# #get all the EMP data from zooper, no winter data
# ZoopsMMx = Zoopsynther(Data_type = "Community",
#                      Sources = c("EMP", "FMWT"),
#                      Size_class = "Macro",
#                      #Months = c(3:10),
#                      Date_range = c("1995-01-01", "2019-12-30"),
#                      Redownload_data = F)
# 
# 
# Zoops2MM = dplyr::filter(ZoopsMMx, !Undersampled, !is.na(Latitude)) %>%
#  mutate(Taxname2 = Taxname, Month = month(Date)) %>%
#   group_by(SampleID, Taxname2, Month, Year, Date, Source, Station, Longitude, Latitude) %>%
#   summarize(CPUE = sum(CPUE)) %>%
#   mutate(Taxname2 = str_replace(Taxname2, "_UnID", ""))
# 
# ZoopsSumMM = Zoops2MM %>%
#   group_by(Station, Source, Date)%>%
#   st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
#   st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
#   st_join(Delta, join=st_intersects)%>% # Add subregions
#   filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
#   st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
#   mutate(Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
#          Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
#                           Month%in%6:8 ~ "Summer",
#                           Month%in%9:11 ~ "Fall",
#                           Month%in%c(12, 1, 2) ~ "Winter",
#                           TRUE ~ NA_character_))%>%
#   group_by(Month, Season, Region, Year, Taxname2)%>%
#   summarise(var_month_mean=mean(CPUE, na.rm = TRUE), N=n())%>%
#   ungroup()%>%
#   as_tibble() %>%
#   complete(nesting(Month, Season), Region, Year, fill=list(N=0)) # Fill in NAs for variable (and 0 for N) for any missing month, subregion, year combinations to make sure all months are represented in each season
# 
# #OK, join water quality to zooplankton!
# 
# ZoopMasterMM = left_join(ZoopsSumMM, WQ_dataNARM) %>%
#   rename(CPUE = var_month_mean) %>%
#   mutate(SalSurf = ec2pss(Conductivity, t=25),
#          logChl = log(Chlorophyll)) %>%
#   filter(!is.na(SalSurf), !is.nan(SalSurf), !is.na(Secchi))
# 
# 
# save(ZoopMasterMM, file = "ZoopMasterMysids.RData")

######################################################3
#mysids 
#get all the EMP data from zooper, 
ZoopsMMx2 = Zoopsynther(Data_type = "Community",
                       Sources = c("EMP", "FMWT", "DOP"),
                       Size_class = "Macro",
                       #Months = c(3:10),
                       Date_range = c("2000-01-01", "2022-12-30"),
                       Redownload_data = F)


Zoops2MM2 = dplyr::filter(ZoopsMMx2, !Undersampled, !is.na(Latitude)) %>%
  mutate(Taxname2 = Taxname, Month = month(Date)) %>%
  group_by(SampleID, Taxname2, Month, Year, Date, Source, Station, Longitude, Latitude, SalSurf, Secchi) %>%
  summarize(CPUE = sum(CPUE)) %>%
  mutate(Taxname2 = str_replace(Taxname2, "_UnID", ""))

ZoopsSumMM2 = Zoops2MM2 %>%
  group_by(Station, Source, Date, SalSurf, Secchi)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
  mutate(Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
         Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
                          Month%in%6:8 ~ "Summer",
                          Month%in%9:11 ~ "Fall",
                          Month%in%c(12, 1, 2) ~ "Winter",
                          TRUE ~ NA_character_))# Fill in NAs for variable (and 0 for N) for any missing month, subregion, year combinations to make sure all months are represented in each season



save(ZoopsSumMM2, file = "ZoopsSumMM2.RData")

####################################################################

#Definite big pile of NOPE on STN and FMWT data.
#Also nope on lag of CHLA

#Scale predictor variables, do some other stuff, log-transormations
ZoopMaster2 = mutate(ZoopMaster, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                     SalSurfsc = scale(log(SalSurf)), logChlsc = scale(logChl),
                     lagChlsc = scale(logChlag),
                     Region = as.factor(Region),
                     Secchisc = scale(Secchi),
                     Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year),
                     Yearf = as.factor(Year)) %>%
  filter(!is.na(Secchi), !is.na(logChl), !(is.na(SalSurf)))
#Scale predictor variables, do some other stuff, log-transormations
ZoopMaster2m = mutate(ZoopMasterm, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                      SalSurfsc = scale(log(SalSurf)), logChlsc = scale(logChl),
                      Region = as.factor(Region),
                      Secchisc = scale(Secchi),
                      Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year))


ZoopMaster2MM2 = mutate(ZoopMasterMM2, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                        SalSurfsc = scale(log(SalSurf)), logChlsc = scale(logChl),
                        Region = as.factor(Region),
                        Secchisc = scale(Secchi),
                        Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year))



#now subset species of interest
Pseudos = filter(ZoopMaster2, Taxname2 == "Pseudodiaptomus")
Eury = filter(ZoopMaster2, Taxname2 == "Eurytemora affinis")
Acan = filter(ZoopMaster2, Taxname2 == "Acanthocyclops")
Acart = filter(ZoopMaster2, Taxname2 == "Acartiella")
Bosmina = filter(ZoopMaster2, Taxname2 == "Bosmina longirostris")
Daphnia = filter(ZoopMaster2, Taxname2 == "Daphnia")
Hyp = filter(ZoopMaster2MM2, Taxname2 == "Hyperacanthomysis longirostris")
Limno = filter(ZoopMaster2m, Taxname2 == "Limnoithona")
Kero = filter(ZoopMaster2m, Taxname2 == "Keratella")
Syn = filter(ZoopMaster2m, Taxname2 == "Synchaeta")

save(ZoopMaster2, ZoopMaster2m, ZoopMaster2MM2, Pseudos, Eury,
     Acan, Bosmina, Daphnia, Hyp, Acart, Kero, Syn, Limno,file = "OrganizedZoops.RData")
