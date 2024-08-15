#evaluate the top-heaviness of the foodweb

library(tidyverse)



##########################################################################
#are we top-heavy?

ZoopsTH = bind_rows(ZoopsSum, ZoopsSumm, ZoopsSumMM2) %>%
  group_by(Region, Season, Taxname2) %>%
  summarize(CPUE = mean(CPUE))

ggplot(data = ZoopsTH, aes(x = Region, y = CPUE, fill = Taxname2))+
  facet_wrap(~Season)+geom_col()

#I need a crosswalk of feeding guilds
#write.csv(unique(ZoopsTH$Taxname2),"data/Taxnames.csv")
guilds = read.csv("data/Taxnames.csv")

ZoopsTH = left_join(ZoopsTH, guilds, by = c("Taxname2"= "Taxname"))

ggplot(data = ZoopsTH, aes(x = Season, y = CPUE, fill = Guild))+
  facet_wrap(~Region)+geom_col(position = "fill")


#hm. take out micros?
ZoopsTH2 = bind_rows(ZoopsSum,  ZoopsSumMM2) %>%
  group_by(Region, Season, Taxname2) %>%
  summarize(CPUE = mean(CPUE)) %>%
  left_join(guilds, by = c("Taxname2"= "Taxname"))


ggplot(data = ZoopsTH2, aes(x = Season, y = CPUE, fill = Guild))+
  facet_wrap(~Region)+geom_col(position = "fill")

#calvin wants it by salinity

ZoopsTHSal = bind_rows(ZoopsSum,  ZoopsSumMM2) %>%
  mutate(Salbin = case_when(SalSurf < 0.5 ~ "1Fresh <0.5",
                            SalSurf >= 0.5 &SalSurf < 2  ~ "2VeryLow 0.5-2",
                            SalSurf >= 2 &SalSurf < 6  ~ "3Low 2-6",
                            SalSurf >= 6 &SalSurf < 12  ~ "4Brackish 6-12",
                            SalSurf >= 12 &SalSurf < 20  ~ "5VeryBrackish 12-20",
                            SalSurf >= 20  ~ "6Salty >20")) %>%
  group_by( Taxname2, Season, Salbin) %>%
  summarize(CPUE = mean(CPUE), N = n()) %>%
  filter(!is.na(Salbin)) %>%
  left_join(guilds, by = c("Taxname2"= "Taxname"))


ggplot(data = ZoopsTHSal, aes(x = Salbin, y = CPUE, fill = Guild))+
  geom_col(position = "fill") +
  facet_wrap(~Season)+
  scale_x_discrete(labels = c("<0.5", "0.5-2", "2-6", "6-12", "12-20", ">20"))
##############################################################################################
#Ideally I"d get biomass rather than CPUE. But that's a bit complicated. I'll do it
#for micro-meso first, then seperately for macro because we have length-wieght regressions

#First export the taxa crosswlk
# tax = bind_rows(ZoopsSumLifestage, ZoopsSumm_lifestage) %>%
#   select(Taxname2, Lifestage) %>%
#   distinct() %>%
#   left_join(guilds, by = c("Taxname2"= "Taxname"))
# write.csv(tax, "data/tax.csv", row.names = FALSE)

guilds = read.csv("data/tax_meso_withbiomass.csv")


ZoopsTH2 = bind_rows(ZoopsSumLifestage,  ZoopsSumm_lifestage) %>%
  left_join(guilds) %>%
  mutate(BPUE = Biomass*CPUE) %>%
  group_by(Region, Season, Taxname2, Lifestage, Guild,  Year) %>%
  summarize(CPUE = mean(CPUE), BPUE = mean(BPUE)) %>%
  ungroup()

#Now do it for macros zooplankton (ugh)
#but we only have lengths for EMP. FMWT is not well organized and I just don't want to deal with it. 
Macinfo = Zoopsynther(Data_type = "Community",
                      Sources = c("EMP"),
                      Size_class = "Macro",
                      #Months = c(3:10),
                      Date_range = c("2000-01-01", "2022-12-30"),
                      Redownload_data = F)

macinfo = select(Macinfo, Source, SizeClass, Volume, SampleID, Date, Station, Year, Latitude, Longitude, SalSurf, Secchi) %>%
  distinct() %>%
  filter(Year >1999, !is.na(Longitude))%>%
  group_by(Station, Source, Date, SalSurf, Secchi)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  st_drop_geometry()%>% # Drop sf geometry column since it's no longer needed
  mutate(Month = month(Date), Year=if_else(Month==12, Year+1, Year), # Move Decembers to the following year
         Season=case_when(Month%in%3:5 ~ "Spring", # Create seasonal variables
                          Month%in%6:8 ~ "Summer",
                          Month%in%9:11 ~ "Fall",
                          Month%in%c(12, 1, 2) ~ "Winter",
                          TRUE ~ NA_character_))


#grab mysid lenghts from EDI
maclen = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.522.11&entityid=43455fa3be9f34fa745280a292801b7a")

maclen2 = mutate(maclen, Date = mdy(SampleDate)) %>%
  filter(Date > ymd("2000-01-01"))
lenweight = read.csv("data/macrotax.csv")

maclen2x = left_join(maclen2, lenweight) %>%
  mutate(Weight_mg = (a*(Size^b))) %>% #find weight and then convert to ug
  mutate(Weight_ug = (Weight_mg*1000),
         CBiomass = case_when(type == "Wet" ~ Weight_ug*.25*.4, type == "Dry" ~ Weight_ug*.4), #convert to carbon weight
         biomassTot = CBiomass*AdjustedFreq) #calculate total in sample

#ok, not total biomass and join to sample infor for BPUE. Fortunately they are all omnivores, so that makes life a little easier

maclen2a = group_by(maclen2x, Taxname2, Date, StationNZ, Guild) %>%
  summarize(Biomass = sum(biomassTot), Count = sum(AdjustedFreq)) %>%
  left_join(macinfo, by = c("StationNZ" = "Station", "Date")) %>%
  mutate(BPUE = Biomass/Volume, CPUE = Count/Volume, Year = year(Date)) %>%
  ungroup()
#now summarize by month and region

maclensum = maclen2a %>%
  group_by(Region, Season, Taxname2, Guild, Year) %>%
  summarize( BPUE = mean(BPUE), CPUE = mean(CPUE)) %>%
  ungroup()


ZoopsTH3 = bind_rows(ZoopsTH2, maclen2a) %>%
  filter(!is.na(Season), !is.na(Region))


ggplot(data = filter(ZoopsTH3, Year >2015), aes(x = Season, y = BPUE, fill = Guild))+
  facet_wrap(~Region)+geom_col(position = "fill")

ggplot(data = filter(ZoopsTH3, Year >2015), aes(x = Year, y = BPUE, fill = Guild))+
  facet_grid(Season~Region)+geom_col()

ggplot(data = filter(ZoopsTH3, Year >2015), aes(x = Year, y = CPUE, fill = Guild))+
  facet_wrap(~Region)+geom_col(position = "fill")

ggplot(data = filter(ZoopsTH3, Year >2015), aes(x = Season, y = CPUE, fill = Guild))+
  facet_wrap(~Region)+geom_col(position = "fill")

#calvin wants it by salinity

ZoopsTHSal = bind_rows(ZoopsSumLifestage,  ZoopsSumm_lifestage) %>%
  left_join(guilds) %>%
  mutate(BPUE = Biomass*CPUE) %>%
  mutate(Salbin = case_when(SalSurf < 0.5 ~ "1Fresh <0.5",
                            SalSurf >= 0.5 &SalSurf < 2  ~ "2VeryLow 0.5-2",
                            SalSurf >= 2 &SalSurf < 6  ~ "3Low 2-6",
                            SalSurf >= 6 &SalSurf < 12  ~ "4Brackish 6-12",
                            SalSurf >= 12 &SalSurf < 20  ~ "5VeryBrackish 12-20",
                            SalSurf >= 20  ~ "6Salty >20")) %>%
  group_by( Taxname2, Season, Salbin, Guild) %>%
  summarize(CPUE = mean(CPUE), N = n(), BPUE = mean(BPUE)) %>%
  filter(!is.na(Salbin)) 



ggplot(data = ZoopsTHSal, aes(x = Season, y = BPUE, fill = Guild))+
  facet_wrap(~Salbin)+geom_col(position = "fill")


#now let's look at change over time


ZoopsTHSalts = bind_rows(ZoopsSumLifestagex,  ZoopsSumm_lifestagex) %>%
  filter(Source == "EMP") %>%
  left_join(guilds) %>%
  mutate(BPUE = Biomass*CPUE) %>%
  mutate(Salbin = case_when(SalSurf < 0.5 ~ "1Fresh <0.5",
                            SalSurf >= 0.5 &SalSurf < 2  ~ "2VeryLow 0.5-2",
                            SalSurf >= 2 &SalSurf < 6  ~ "3Low 2-6",
                            SalSurf >= 6 &SalSurf < 12  ~ "4Brackish 6-12",
                            SalSurf >= 12 &SalSurf < 20  ~ "5VeryBrackish 12-20",
                            SalSurf >= 20  ~ "6Salty >20")) %>%
  group_by( Taxname2, Season, Salbin, Guild, Year) %>%
  summarize(CPUE = mean(CPUE), N = n(), BPUE = mean(BPUE)) %>%
  filter(!is.na(Salbin)) 


ggplot(data = ZoopsTHSalts, aes(x = Year, y = BPUE, fill = Guild))+
  facet_wrap(~Salbin)+geom_col(position = "fill")

#now a version that isn't summarized to put in a model

Zoopsalt =  bind_rows(ZoopsSumLifestagex,  ZoopsSumm_lifestagex) %>%
  filter(Source == "EMP") %>%
  left_join(guilds) %>%
  filter(!is.na(Biomass)) %>%
  mutate(BPUE = Biomass*CPUE) %>%
  mutate(Salbin = case_when(SalSurf < 0.5 ~ "1Fresh <0.5",
                            SalSurf >= 0.5 &SalSurf < 2  ~ "2VeryLow 0.5-2",
                            SalSurf >= 2 &SalSurf < 6  ~ "3Low 2-6",
                            SalSurf >= 6 &SalSurf < 12  ~ "4Brackish 6-12",
                            SalSurf >= 12 &SalSurf < 20  ~ "5VeryBrackish 12-20",
                            SalSurf >= 20  ~ "6Salty >20"))%>%
  pivot_wider(id_cols = c(Year, Season, Salbin, SalSurf), names_from = Guild, 
              values_from = BPUE, values_fill = 0, values_fn = sum) %>%
  ungroup()


zm = as.matrix(Zoopsalt[,5:7])
zm2 = zm/rowSums(zm)
  
hist(Zoopsalt$Predator)

#Huh... I could use teh same time of analysis i did for microcystis (maybe)
#but they aren't ordered categories... unless omnivores are sorta halfway between predtaors and herbivores?
#no, that doesn't make sense

M1 = brm(bind_cols(Zoopsalt$Predator, Zoopsalt$Omnivore, Zoopsalt$Herbivore) ~ SalSurf + Season + (1|Year), family  = categorical, data = ZoopsTHSalts,
         iter = 100,   backend = "cmdstanr", normalize = FALSE,
         #control = list(max_treedepth = 15),
         chains = 2, cores=4, threads = threading(2))

#Ugh, maybe I'll just do a permanova by salinity bin 
library(vegan)

zoopmat = ZoopsTHSalts %>%
  filter(!is.na(Guild), !is.na(BPUE)) %>%
  pivot_wider(id_cols = c(Year, Season, Salbin), names_from = Guild, 
              values_from = BPUE, values_fill = 0, values_fn = sum) %>%
  ungroup() %>%
  mutate(Salbin = as.ordered(Salbin)) %>%
  arrange(Year, Season, Salbin)

zm = as.matrix(zoopmat[,4:6])

zm2 = zm/rowSums(zm)

a1 = adonis2(zm ~ Salbin+ Season, data = zoopmat)

summary(a1)

mds1 = metaMDS(zm)
PlotNMDS(mds1, group = "Salbin", data = zoopmat)
PlotNMDS2(mds1, lines = "Year",  data = zoopmat, group = "Salbin")


#####################################################################http://127.0.0.1:15777/graphics/plot_zoom_png?width=925&height=545
#now by region and salinity instead of year
Zoopreg =  bind_rows(ZoopsSumLifestagex,  ZoopsSumm_lifestagex) %>%
  filter(Source == "EMP", Year >2015) %>%
  left_join(guilds) %>%
  mutate(BPUE = Biomass*CPUE, Sal = round(SalSurf), Region = as.factor(Region)) %>%
  filter(!is.na(BPUE)) %>%
  pivot_wider(id_cols = c(Year, Region, Sal), names_from = Guild, 
              values_from = BPUE, values_fill = 0, values_fn = sum) %>%
  ungroup() 
Zoopreg = Zoopreg[which(rowSums(Zoopreg[,4:6]) != 0),]


zmx = as.matrix(Zoopreg[,4:6])


zmx2 = zmx/rowSums(zmx)

mds1 = metaMDS(zmx2)
PlotNMDS(mds1, group = "Region", data = Zoopreg)
PlotNMDS2(mds1, lines = "Sal",  data = Zoopreg, group = "Region")

zoopreg2 = pivot_longer(Zoopreg,cols = c(Herbivore, Predator, Omnivore), names_to = "Guild", values_to = "BPUE")

ggplot(zoopreg2, aes(x = Sal, y = BPUE, fill = Guild)) + geom_col(position = "fill")+
  facet_wrap(~Region)

ggplot(zoopreg2, aes(x = Year, y = BPUE, fill = Guild)) + geom_col(position = "fill")+
  facet_wrap(~Region)

######################################################################

hist(zm2[,1])
hist(zm2[,2])
hist(zm2[,3])

#I think i need a MARSS model, which I've never done before and now is not the time.

#or maybe it is the time...
library(MARSS)

covariates <- rbind(
  Sal = as.numeric(zoopmat$Salbin)
) %>%
  zscore()

Q <- U <- x0 <- "zero"
B <- Z <- "identity"
d <- covariates
A <- "zero"
D <- "unconstrained"
y <- t(zm2) # to show relationship between dat & the equation
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A,
  D = D, d = d, x0 = x0
)
mar1 <- MARSS(y, model = model.list)

tidy(mar1)

plot(mar1)

####################################################################
#ok, let's just do a model with proportion of biomass as each group by salinity and year

Zoopsalt =  bind_rows(ZoopsSumLifestagex,  ZoopsSumm_lifestagex) %>%
  filter(Source == "EMP") %>%
  left_join(guilds) %>%
  filter(!is.na(Biomass)) %>%
  mutate(BPUE = Biomass*CPUE) %>%
  mutate(Salbin = case_when(SalSurf < 0.5 ~ "1Fresh <0.5",
                            SalSurf >= 0.5 &SalSurf < 2  ~ "2VeryLow 0.5-2",
                            SalSurf >= 2 &SalSurf < 6  ~ "3Low 2-6",
                            SalSurf >= 6 &SalSurf < 12  ~ "4Brackish 6-12",
                            SalSurf >= 12 &SalSurf < 20  ~ "5VeryBrackish 12-20",
                            SalSurf >= 20  ~ "6Salty >20"))%>%
  pivot_wider(id_cols = c(Year, Season, Salbin, SalSurf), names_from = Guild, 
              values_from = BPUE, values_fill = 0, values_fn = sum) %>%
  ungroup()


zm = as.matrix(Zoopsalt[,5:7])
zm2 = zm/rowSums(zm)



zoopProp = bind_cols(Zoopsalt[,-c(5:7)], zm2)

predlm = lm(log(Predator+0.01) ~ Year*SalSurf, data = zoopProp)
summary(predlm)
plot(predlm)
#Gross.

omlm = lm(log(Omnivore+0.01) ~ Year*SalSurf, data = zoopProp)
summary(omlm)
plot(omlm)
#Gross.

herblm = lm(log(Herbivore+0.01) ~ Year*SalSurf, data = zoopProp)
summary(herblm)
plot(herblm)


library(effects)

plot(allEffects(predlm))
#OK, so there is definitely an effect of salinity, but that effect doesn't change over time for predators

plot(allEffects(omlm))
#the effect of salinity does change over time for omnivores, more in recent years

plot(allEffects(herblm))
#Fewer herbivores at higher salinities and fewer over time.


#is it possible to combine all these in one analysis with a three-way interaction? Would that suck?

zoopProp2 = zoopProp %>%
  mutate(ID = paste(Year, Season, SalSurf)) %>%
  pivot_longer(cols = c(Herbivore, Predator, Omnivore), names_to = "Guild",
                         values_to = "Proportion") %>%
  filter(!is.nan(Proportion)) %>%
  mutate(Prop = case_when(Proportion ==1 ~ Proportion -0.0001,
                          Proportion ==0 ~ Proportion + 0.0001,
                          TRUE ~ Proportion))

alllm = lm(log(Proportion*100+1) ~ Guild*SalSurf*Year, data = zoopProp2)
summary(alllm)
plot(allEffects(alllm))

alllmeff = allEffects(alllm, xlevels=list(SalSurf=c(0.01, 1, 2, 3, 5, 7, 9, 11, 13, 
                                               15, 17, 19, 21, 23, 25)))

effectsdat = data.frame(alllmeff$`Guild:SalSurf:Year`$x,
                        Fit = alllmeff$`Guild:SalSurf:Year`$fit)

ggplot(effectsdat, aes(x = SalSurf, y = exp(Fit)-1)) +
  facet_grid(Guild~Year)+
  geom_point()+
  geom_line()+
  coord_cartesian(ylim = c(0,100))
#########################################################################
#maybe a simplex regression will get something better?
library(simplexreg)

sim1 = simplexreg(Prop ~ Guild,  link = "logit", 
                 
                   id = ID, data = as.data.frame(zoopProp2))


#this is frustrating, try by year and salinity bin to make life easier.


zoopsalt2 = bind_rows(ZoopsSumLifestagex,  ZoopsSumm_lifestagex) %>%
  filter(Source == "EMP") %>%
  left_join(guilds) %>%
  filter(!is.na(Biomass)) %>%
  mutate(BPUE = Biomass*CPUE) %>%
  mutate(Salbin = case_when(SalSurf < 0.5 ~ "1Fresh <0.5",
                            SalSurf >= 0.5 &SalSurf < 2  ~ "2VeryLow 0.5-2",
                            SalSurf >= 2 &SalSurf < 6  ~ "3Low 2-6",
                            SalSurf >= 6 &SalSurf < 12  ~ "4Brackish 6-12",
                            SalSurf >= 12 &SalSurf < 20  ~ "5VeryBrackish 12-20",
                            SalSurf >= 20  ~ "6Salty >20"))%>%
  pivot_wider(id_cols = c(Year, Salbin), names_from = Guild, 
              values_from = BPUE, values_fill = 0, values_fn = sum) %>%
  ungroup()


zm = as.matrix(zoopsalt2[,3:5])
zm2 = zm/rowSums(zm)



zoopProp3 = bind_cols(zoopsalt2[,-c(3:5)], zm2) %>%
  mutate(ID = paste(Year, Salbin)) %>%
  pivot_longer(cols = c(Herbivore, Predator, Omnivore), names_to = "Guild",
               values_to = "Proportion") %>%
  filter(!is.nan(Proportion), !is.na(Salbin), Year >1979) %>%
  mutate(Prop = case_when(Proportion ==1 ~ Proportion -0.0001,
                          Proportion ==0 ~ Proportion + 0.0001,
                          TRUE ~ Proportion)) 

sim2 = simplexreg(Prop ~ Year*Guild*Salbin,  link = "logit",
                  id = ID, data = as.data.frame(filter(zoopProp3, Salbin != "6Salty >20")), na.rm =T)
summary(sim2, type = "appstdPerr")
sim2
#I cant use the effects package, I guess, which isa bummer
plot(sim2, type = "residuals")
test = residuals(sim2)
dat = data.frame(sim2$model, residuals =  residuals(sim2), 
                 appstpr = sim2$appstdPerr, stdper = sim2$stdPerr,
                 predict = sim2$predict,
                 Probability = exp(sim2$predict)/(1+exp(sim2$predict)))

ggplot(dat, aes(x = Year, y = Probability, color = Guild)) +
  geom_line()+
  facet_wrap(~Salbin)

ggplot(dat, aes(x = Year, y = Probability, fill = Guild)) +
  geom_area(position = "fill")+
  facet_wrap(~Salbin)

ggplot() +
  geom_col(data =filter(zoopProp3, Salbin != "6Salty >20"),
           aes(x = Year, y = Proportion, fill = Guild), alpha = 0.5)+
  geom_area(data =dat, aes(x = Year, y = Probability, color = Guild), alpha =0, position = "fill",
            size =1)+
  facet_wrap(~Salbin)+
  scale_fill_manual(values = c("green", "slateblue", "darkred"))+
  scale_color_manual(values = c("green", "slateblue", "darkred"))+
  theme_bw()

