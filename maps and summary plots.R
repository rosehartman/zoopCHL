#maps of sample location

library(tidyverse)
library(sf)
library(deltamapr)
library(ggspatial)
#detach("package:MASS")

load("data/filteredzoopData.RData")
load("data/ZoopsSum.RData")
load("data/ZoopsSumm.RData")
load("data/lilzoopsmean.RData")
load("data/ZoopsSumMM2.RData")

zoopstas = select(ZoopsSum, Source, Station, Latitude, Longitude, Region) %>%
  filter(Source != "DOP", !Station %in% c("NZEZ6",  "NZEZ2", "NZEZ6SJR",  "NZEZ2SJR")) %>%
  distinct() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

ggplot()+
  geom_sf(data= WW_Delta) +
  geom_sf(data = Delta, aes(fill = Region), alpha = 0.5)+
  geom_sf(data = zoopstas, aes(shape = Source))+
  coord_sf(xlim = c(-122.2, -121.3), ylim = c(37.8, 38.6))


#expor thte stations df and edit it for overlapping stations
#write.csv(zoopstas, "data/zoopstas.csv", row.names = FALSE)

zoopstas = read_csv("data/zoopstas2.csv")%>%
  st_as_sf(coords = c("Longitude", "Latitue"), crs = 4326)


ggplot()+
  geom_sf(data= WW_Delta) +
  geom_sf(data = Delta, aes(fill = Region), alpha = 0.5)+
  geom_sf(data = zoopstas, aes(shape = Surveys), fill = "blue")+
  coord_sf(xlim = c(-122.2, -121.3), ylim = c(37.8, 38.6))+
  theme_bw()+
  scale_shape_manual(values = c(21,22,3,4,16,15,14))

ggsave("plots/zoopmap.tiff",device = "tiff", width = 7, height =7)
###################################################################################
#plot of samplign effort

mesomacromicro = bind_rows(mutate(ZoopsSum, size = "meso"), mutate(ZoopsSumm, size = "micro"),
                           mutate(ZoopsSumMM2, size = "macro"))

memami = select(mesomacromicro, SampleID, size, Source, Region, Year) %>%
  distinct() %>%
  mutate(Size = factor(size, levels = c("micro", "meso", "macro", "chl"), labels = c("Micro (50-60 micron mesh)",
                                                                              "Meso (150-160 micron mesh)",
                                                                              "Macro (500 micron mesh)",
                                                                              "Chlorophyll"))) %>%
  group_by(Region, Year, Size) %>%
  summarize(N = n())
#now chlorophyll sampling effort

chl = filter(WQ_dataNARM, Year>1999) %>%
  group_by(Year, Region) %>%
  summarize(N = sum(N_Chlorophyll, na.rm =T)) %>%
  mutate(Size = "Chlorophyll")

allsamps = bind_rows(memami, chl)

ggplot(allsamps, aes(x = Year, y = N, fill = Region)) +
  facet_wrap(~Size)+
  geom_col(color = "grey40")+
  theme_bw()+
  ylab("Number of Samples")+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "inside", legend.position.inside = c(.1, .8),
        legend.box.background = element_rect(color = "black"))

ggsave("plots/samplingeffort.tiff", device = "tiff", width =6, height =6)



ggplot(chl, aes(x = Year, y = N_chl, fill = Region))+
  geom_col(color = "grey40")+
  theme_bw()+
  ylab("Number of Chlorophyll Samples")+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "inside", legend.position.inside = c(.15, .8),
        legend.box.background = element_rect(color = "black"))


ggsave("plots/chlasamplingeffort.tiff", device = "tiff", width =6, height =5)

#######################################################################################

#ok, what about those chlorophyll stations?

Chlstas<-wq(End_year=2022, Start_year = 1999,
         Sources = c("EMP", "USGS_SFBS", "USGS_CAWSC", "NCRO")) %>%
  select(Source, Station, Latitude, Longitude)%>%
  filter(!Station %in% c("EZ6",  "EZ2", "EZ6-SJR",  "EZ2-SJR"), !is.na(Latitude)) %>%
  distinct() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(Delta)) %>%
  st_join(Delta) %>%
  filter(!is.na(Region))

Delta2 = filter(Delta, SubRegion != "Upper San Joaquin River") %>%
  group_by(Region) %>%
  summarize(geometry= st_union(geometry))

ggplot()+
  geom_sf(data= WW_Delta) +
  geom_sf(data = Delta2, alpha = 0.5)


ggplot()+
  geom_sf(data= WW_Delta) +
  geom_sf(data = Delta2, aes(fill = Region), alpha = 0.5)+
  geom_sf(data = Chlstas, aes(shape = Source), fill = "darkgreen")+
  coord_sf(xlim = c(-122.2, -121.3), ylim = c(37.8, 38.6))+
  theme_bw()+
  scale_shape_manual(values = c(21,22,9,23))
ggsave("plots/chlmap.tiff",device = "tiff", width = 7, height =7)
#OK, that needs work


####################################################
#can I plot both together? Would that be terrible?
zoopstas = mutate(zoopstas, Type = "Zooplankton")
Chlstas = mutate(Chlstas, Type = "Chlorophyll")

ggplot()+
  geom_sf(data= WW_Delta, color = "grey70", fill = "grey90") +
  geom_sf(data = Delta2, aes(fill = Region), alpha = 0.3)+
  geom_sf(data = zoopstas, aes(shape = Type), size =2)+
  geom_sf(data = Chlstas,  aes(shape = Type))+
  coord_sf(xlim = c(-122.15, -121.35), ylim = c(37.8, 38.55))+
  scale_shape_manual(values = c(12, 20), name = "Sample Type")+
  geom_sf_label(data = Delta2, aes(label = Region), 
                nudge_x = c(-.01, -0.07, 0.05,0,0), nudge_y = c(-.1, .1, -.2, -.06, .05))+
  scale_fill_brewer(palette = "Dark2", guide = NULL)+
  theme_bw()+
  annotation_scale()+
  annotation_north_arrow(which_north = "grid", location = "tr")+
  theme(legend.position = "inside", legend.position.inside = c(.15, .9))+
  ylab(NULL)+xlab(NULL)


ggsave("plots/chlzoopmap.tiff",device = "tiff", width = 7, height =7)


############################################################
#months of highest abundance
taxa = read_excel("Data/taxa.xlsx") %>%
  mutate(taxa, Taxlifestage = str_remove(Taxlifestage, " Adult")) %>%
           select(Taxlifestage, FeedingGuild) %>%
           distinct()

allzoopmonth = group_by(allzoops, Month, Taxlifestage, Taxname2, Lifestage) %>%
  summarize(CPUE = mean(CPUE), sal = median(SalSurf)) %>%
  mutate(Taxlifestage = str_remove(Taxlifestage, " NA"), Taxlifestage = str_remove(Taxlifestage, " Adult")) %>%
  left_join(taxa)

ggplot(allzoopmonth, aes(x = Month, y = Taxlifestage, fill = CPUE)) +
  geom_tile() +
  facet_wrap(~FeedingGuild)

taxalabels = c(expression(paste(italic("Acanthocyclops"), " sp.")),
               expression(paste(italic("A. sinensis"), " adult")),
               expression(paste(italic("A. sinensis"), " juv.")),
               expression(paste(italic("B. longirostris"))),
               expression(paste(italic("Daphnia"), " sp.")),
               expression(paste(italic("E. carolleeae"), " adult")),
               expression(paste(italic("E. carolleeae"), " juv.")),
               expression(paste(italic("H. longirostris"))),
               expression(paste(italic("Keratella"), " sp.")),
               expression(paste(italic("Limnoithona"), " sp. adult")),
               expression(paste(italic("Limnoithona"), " sp. juv.")),
               expression(paste(italic("P. forbesi"), " adult")),
               expression(paste(italic("Pseudodiaptomus"), " sp. juv.")),
               expression(paste(italic("Synchaeta"), " sp.")),
               expression(paste(italic("Tortanus"), " sp. adult")),
               expression(paste(italic("Tortanus"), " sp. juv")))

allzoopmonth = left_join(allzoopmonth, taxalabels)

ggplot(allzoopmonth, aes(x = Month, y = Taxlifestage, fill = sal)) +
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(name = "Median\nSalinity") +
  #scale_y_discrete(labels = taxalabels$Taxlab)+
  scale_x_continuous(breaks = c(1,3,5,7,9,11), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
facet_wrap(~FeedingGuild, nrow =3, scales = "free_y")+
  ylab(NULL)

###############################################
#how do I get the good looking labels?

library(ggh4x)
library(scales)


g <- 
  ggplot(allzoopmonth, aes(x = Month, y = Taxlifestage, fill = sal)) +
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(name = "Median\nSalinity") +
  scale_y_discrete(labels = taxalabels$Taxlab)+
  scale_x_continuous(breaks = c(1,3,5,7,9,11), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  facet_wrap(~FeedingGuild, nrow =3, scales = "free_y", strip.position = "right")+
  facetted_pos_scales(y = list(scale_y_discrete(labels = taxlabelsH),
                            scale_y_discrete(labels = taxlabelsO),
                            scale_y_discrete(labels = taxlabelsP)))+
  ylab(NULL)
g

ggsave("plots/monthlysalinity.tiff", device = "tiff", width =7, height =6)


#now let's try something a bit different

ggplot(allzoopmonth, aes(x = Month, y = CPUE, color = Taxlifestage))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Taxlifestage, scales = "free")+
  scale_color_discrete(guide = NULL)

allzoops = mutate(allzoops, Month2 = case_when(Taxlifestage == "Acanthocyclops_UnID NA"  ~ factor(Month, levels = c(12,1,2,3,4), 
                                                                                                  labels = c("0", "1", "2", "3", "4")),
                                               TRUE ~ factor(Month, levels = sort(unique(Month)))),
                  Month3 = as.numeric(Month2))


ggplot(allzoops, aes(x = Month3, y = log(CPUE), color = Taxlifestage))+ 
  geom_point(alpha = 0.1)+
  geom_smooth()+
  facet_wrap(~Taxlifestage, scales = "free")+
  scale_color_discrete(guide = NULL)


##############################################
#some quick checks
library(zooper)

mactest = Zoopsynther(Data_type = "Community", Size_class = "Macro", Years = c(2015:2020),
                      Time_consistency = TRUE)
unique(mactest$Taxlifestage)

topcatch = group_by(mactest, Taxlifestage) %>%
  summarize(CPUE = sum(CPUE))

mictest = Zoopsynther(Data_type = "Community", Size_class = "Micro")
unique(mictest$Taxlifestage)

topcatch = group_by(mactest, Taxlifestage) %>%
  summarize(CPUE = sum(CPUE))
