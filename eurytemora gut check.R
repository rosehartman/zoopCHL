#quick check for Brian - do zooplankton ever really drop out of an area?

library(zooper)
library(tidyverse)
library(lubridate)



#get all the EMP and 20mm data from zooper
#If I'm using GAMs, I can do winter data after all.
Zoops = Zoopsynther(Data_type = "Community",
                    Sources = c("EMP", "20mm", "FMWT", "STN"),
                    Size_class = "Meso",
                    #  Months = c(3:10),
                    Date_range = c("2000-01-01", "2022-12-30"),
                    Redownload_data = F)

euryall = filter(Zoops, Genus == "Eurytemora") %>%
  mutate(Year = year(Date), Month = month(Date), doy = yday(Date)) %>%
  group_by(Station, Year, SampleID, Month, Date,doy, SalSurf, Genus) %>%
  summarize(CPUE = sum(CPUE)) %>%
  mutate(Salbin = case_when(SalSurf < 1 ~ "1 Fresh",
                            SalSurf >= 1 & SalSurf <6 ~ "2-6 LSZ",
                            SalSurf >=6 & SalSurf <12 ~ "6-12 Brackish",
                            SalSurf >= 12 ~ ">12 Salty"))


hist(euryall$CPUE)

ggplot(euryall, aes(x = doy, y = log(CPUE+1), color = SalSurf)) + geom_point()

salbins = group_by(euryall, Month, Salbin) %>%
  summarize(Zero = length(CPUE[which(CPUE==0)]), NotZero = length(CPUE[which(CPUE>0)])) %>%
  mutate(Percent0 = Zero/(NotZero+Zero)*100)

ggplot(salbins, aes(x = Month, y = Percent0)) + geom_col()+ facet_wrap(~Salbin)


#check on flow-abundance relationships       
load("data/Dayflow_allw2023.RData")


springzoops = filter(ZoopsSum, month(Date) %in% c(3,4,5,6), CPUE < 20000,
                     Region %in% c("Suisun Bay", "Suisun Marsh","Confluence")) %>%
  left_join(select(Dayflow, Date, OUT)) 

ggplot(springzoops, aes(x = OUT, y = CPUE)) +geom_point()+
  facet_wrap(~Taxname2)

ggplot(springzoops, aes(x = log(OUT), y = log(CPUE+1))) +geom_smooth(method = "lm")+
  facet_wrap(~Taxname2, scales = "free_y")


#macrozoops
springzoops = filter(ZoopsSumMM2, month(Date) %in% c(3,4,5,6), CPUE < 20000,
                     Region %in% c("Suisun Bay", "Suisun Marsh","Confluence"))%>%
  left_join(select(Dayflow, Date, OUT))

ggplot(springzoops, aes(x = OUT, y = CPUE)) +geom_point()+
  facet_wrap(~Taxname2)

ggplot(springzoops, aes(x = log(OUT), y = log(CPUE+1))) +geom_smooth(method = "lm")+
  facet_wrap(~Taxname2, scales = "free_y")


#quick check of limnoithona

limnoall = Zoopsynther(Data_type = "Taxa", Taxa = "Limnoithona")

limnoall1 = filter(limnoall, CPUE !=0)
ggplot(limnoall1, aes(x = SalSurf, y = CPUE, color = Taxname))+
  geom_smooth()


ggplot(limnoall1, aes(x = Date, y = CPUE, color = Taxname))+
  geom_smooth()+
  facet_wrap(~Taxname, scales = "free")

limnoave = group_by(limnoall, Year, Taxname) %>%
  summarize(CPUE = mean(CPUE))
