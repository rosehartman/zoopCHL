#New version of the analysis for reboot of project.

#just use the 5 months of greatest abundance
#three taxa from each functional group
#linear models are better.

library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)
library(glmmTMB)
library(DHARMa)
library(performance)
library(car)
library(zooper)
library(effects)
library(MuMIn)
library(wql)
library(deltamapr)
library(sf)
library(wql)
library(corrplot)

load("data/ZoopsSum.RData")
load("data/lilzoopsmean.RData")

#Note: i tried a lag of chlorophyll a earlier and it didn't work very well
#start with Pseudodiaptomus

Regions<-read_csv("data/Rosies_regions.csv")

Delta<-deltamapr::R_EDSM_Subregions_Mahardja_FLOAT%>%
  filter(SubRegion%in%unique(Regions$SubRegion))%>%  #Filter to regions of interest
  dplyr::select(SubRegion)

Regs = unique(Regions[,c(1,5)])
Delta = merge(Delta, Regs)

#I need to re-run this sorting out the RMEL things and using BIC instead of AIC.

#Maybe filter by time of year and salinity range things are most abundant
Filters_LS = function(taxname, lifestage, data = ZoopsSum){
  data = dplyr::filter(data, Taxname2 == taxname, Lifestage == lifestage) #filter to taxon of interest
  
  #summarize by month and select the highest 5 months
  monthsdat = data %>%
    group_by(Month, Taxname2, Lifestage) %>%
    summarize(CPUE = mean(CPUE))%>%
    arrange(desc(CPUE)) %>%
    ungroup() %>%
    slice(1)
  
  months = c(monthsdat$Month, monthsdat$Month+1, monthsdat$Month+2,monthsdat$Month-1,monthsdat$Month-2)
  months = case_when(months < 1 ~ months+12,
                     months > 12 ~ months-12,
                     TRUE ~ months)
  
  
  #find the salinities with the most catch
  test = dplyr::filter(data,  CPUE !=0)
   limsSurf = quantile(test$SalSurf, probs = c(0.05, 0.95), na.rm =T)
   
   #filter data to months and salinities with highest catch, take mean by region and month, and summarize
  datafiltered = dplyr::filter(data, Month %in% months, 
                 between(SalSurf, limsSurf[1], limsSurf[2])) %>%
    group_by(Month, Region, Season, Year, Taxname2, Lifestage) %>%
    summarize(CPUE = mean(CPUE, na.rm =T)) %>%
    left_join(WQ_dataNARM, by = c("Month", "Year", "Region")) %>%
    left_join(lilzoopsmean) %>%
    mutate_all( ~ case_when(!is.nan(.x) ~ .x)) %>%
    ungroup() %>%
    mutate(SalSurf = ec2pss(Conductivity/1000, t=25),
           logChl = log(Chlorophyll),
           logCPUE = log(CPUE +1),
           rCPUE = round(CPUE))
  datafiltered$Secchisc = scale(datafiltered$Secchi)
  datafiltered$logChlsc = scale(datafiltered$logChl)
  datafiltered$SalSurfsc = scale(datafiltered$SalSurf)
  datafiltered$lilzoops = scale(log(datafiltered$lilzoopbiomass+1))
  return(datafiltered)
}


Pseudox_juv = Filters_LS("Pseudodiaptomus_UnID", "Juvenile", ZoopsSumLifestage) %>%
  filter(!is.na(SalSurf), !is.na(Secchi),  !is.na(lilzoops))
Pseudox_ad = Filters_LS("Pseudodiaptomus forbesi", "Adult", ZoopsSumLifestage) %>%
  filter(!is.na(SalSurf), !is.na(Secchi), !is.na(lilzoops))
Eury_juv = Filters_LS("Eurytemora affinis", "Juvenile", data =ZoopsSumLifestage)%>%
  filter(!is.na(Secchi))
Eury_ad = Filters_LS("Eurytemora affinis", "Adult", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi))
Acan_juv = Filters_LS("Acanthocyclops_UnID", "Juvenile", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi))
#Ok, we don't do juvenile acanthocyclops, so never mind

Acan_ad = Filters_LS("Acanthocyclops_UnID", "Adult", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi))

Limno_juv = Filters_LS("Limnoithona", "Juvenile", data = ZoopsSumm_lifestage)%>%
  filter(!is.na(Secchi))
Limno_ad = Filters_LS("Limnoithona", "Adult", data = ZoopsSumm_lifestage)%>%
  filter(!is.na(Secchi))

Tort_juv = Filters_LS( "Tortanus_UnID" , "Juvenile", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi))
Tort_ad = Filters_LS( "Tortanus_UnID" , "Adult", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi))


Acart_juv = Filters_LS( "Acartiella_UnID" , "Juvenile", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi), Year >2005)
Acart_ad = Filters_LS( "Acartiella_UnID" , "Adult", data = ZoopsSumLifestage)%>%
  filter(!is.na(Secchi), Year >2005)


#####################################################33
#quick plot of which months were included for wchich taxa

allzoops = bind_rows(Acart_juv, Eury_juv,  
                     Pseudox_juv,  Tort_juv, Limno_juv, Limno_ad, Eury_ad,
                     Tort_ad, Pseudox_ad, Acan_ad, Acart_ad)

allzoopmonth = group_by(allzoops, Month, Taxname2, Lifestage) %>%
  summarize(N= n(), CPUE = sum(CPUE, na.rm =T), logCPUE = sum(logCPUE, na.rm=T))


ggplot(allzoopmonth, aes(x = Month, y = paste(Taxname2, Lifestage), fill = logCPUE)) + geom_tile()


ggplot(allzoops, aes(x = SalSurf, y = CPUE)) + geom_point()+
  geom_smooth()+ 
  facet_wrap(~Taxname2, scales = "free")
#################################################################
#pesudodiaptomus, juveniles and adults


#zero inflated model - adults
p5ad = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + lilzoops+  (1|Month)+ (1|Year),
             ziformula = ~log(SalSurf),
             data= Pseudox_ad, family = nbinom2, na.action = "na.fail")
summary(p5ad)


#go through all the modesl
pdad = dredge(p5ad, rank = "BIC", trace = 2,extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
pBestad = get.models(pdad, 1)[[1]]
summary(pBestad)

pBestad = update(pBestad, REML = TRUE)
summary(pBestad)




#now the juvenikels
p5juv = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + lilzoops+  (1|Month)+ (1|Year),
               ziformula = ~log(SalSurf),
               data= Pseudox_juv, family = nbinom2, na.action = "na.fail")
summary(p5juv)



pdjuv = dredge(p5juv, rank = "BIC", trace = 2,extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
pBestjuv = get.models(pdjuv, 1)[[1]]
summary(pBestjuv)

pBestjuv = update(pBestjuv, REML = TRUE)
summary(pBestjuv)
plot(simulateResiduals(pBestjuv))

#Huh, well, that's odd.

hist(Pseudox_ad$logCPUE)
hist(Pseudox_juv$logCPUE)

#hmmm.
######################################################


#Eurytemora

#adults
Eury_ad2 = filter(Eury_ad, !is.na(lilzoops))
e5ad = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + lilzoops+ Secchisc + (1|Month) + (1|Year), 
             zi = ~ log(SalSurf) , data = Eury_ad2, 
             family= "nbinom2", na.action = "na.fail")

summary(e5ad)
plot(allEffects(e5ad))

summary(e5ad)
edad = dredge(e5ad,  trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
eBestad = get.models(edad, 1)[[1]]
eBestad = update(eBestad, REML = T)
summary(eBestad)
plot(simulateResiduals(eBestad))
#Huh. 


#juvenikels
Eury_juv2 = filter(Eury_juv, !is.na(lilzoops))
e5juv = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + lilzoops+ Secchisc + (1|Month) + (1|Year), 
               zi = ~ log(SalSurf) , data = Eury_juv2, 
               family= "nbinom2", na.action = "na.fail")

summary(e5juv)

plot(allEffects(e5juv))

summary(e5juv)
edjuv = dredge(e5juv,  trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
eBestjuv = get.models(edjuv, 1)[[1]]
eBestjuv = update(eBestjuv, REML = T)
summary(eBestjuv)
plot(simulateResiduals(eBestjuv))
#Huh. Little zoops is in the juveniles one too. Odd. 
plot(allEffects(eBestjuv))

ggplot(Eury_juv, aes(x = lilzoops, y = logCPUE)) + geom_point()+geom_smooth()
ggplot(Eury_juv, aes(x = Month, y = log(CPUE)))+ geom_point()

Euryall = bind_rows(Eury_juv, Eury_ad)

ggplot(Euryall, aes(x = Month, y = log(CPUE), color = Lifestage))+
  geom_point()+
  facet_grid(Year~Region)

#hm. Something is odd in the north in recent years
#######################################################################


#Acartiella
#Acartiella copepodids weren't seperated until 2006, so I'm going to subset
#to just that time period.
Acart_ad2 = filter(Acart_ad, !is.na(lilzoops))
a5ad = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + lilzoops+ logChlsc+ + (1|Year), 
             zi = ~ log(SalSurf), data = Acart_ad2, 
             family= "nbinom2", na.action = "na.fail")
summary(a5ad)
adad = dredge(a5ad,trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
aBestad = get.models(adad, 1)[[1]]
aBestad = update(aBestad, REML =T)
summary(aBestad)

###################################################################
#Tortanus

Tort_ad1 = filter(Tort_ad, !is.na(lilzoops))
t5ad = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + lilzoops + (1|Month) + logChlsc + (1|Year), 
             zi = ~ log(SalSurf), data = Tort_ad1,
             family= "nbinom2", na.action = "na.fail")

summary(t5ad)
tdad = dredge(t5ad, trace = 2, rank = "BIC",  
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
tBestad = get.models(tdad, 1)[[1]]
tBestad = update(tBestad, REML = T)
summary(tBestad)
#oh, nice, when it's just adults lilzoops are better than chlorophyll!

#Tortatuns juveniles
Tort_juv1 = filter(Tort_juv, !is.na(lilzoops))
t5juv = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + lilzoops + (1|Month) + logChlsc + (1|Year), 
               zi = ~ log(SalSurf), data = Tort_juv1,
               family= "nbinom2", na.action = "na.fail")

summary(t5juv)
tdjuv = dredge(t5juv, trace = 2, rank = "BIC",  
              extra = list("Rs" = function(x) {
                s <- performance(x)
                c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
              }))
tBestjuv = get.models(tdjuv, 1)[[1]]
tBestjuv = update(tBestjuv, REML = T)
summary(tBestjuv)
#whereas juveniles like chlorophyll

################################################

#limnoithona adulat
l5ad = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + lilzoops+ logChlsc +(1|Year), 
             zi = ~ log(SalSurf), data = Limno_ad,
             family= "nbinom2", na.action = "na.fail")

summary(l5ad)
ldad = dredge(l5ad, trace = 2, rank = "BIC",  
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
lBestad = get.models(ldad, 1)[[1]]
lBestad = update(lBestad, REML =T)
summary(lBestad)

#Juvenile limnos
l5juv = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + lilzoops+ logChlsc +(1|Year), 
               zi = ~ log(SalSurf), data = Limno_juv,
               family= "nbinom2", na.action = "na.fail")

summary(l5juv)
ldjuv = dredge(l5juv, trace = 2, rank = "BIC",  
              extra = list("Rs" = function(x) {
                s <- performance(x)
                c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
              }))
lBestjuv = get.models(ldjuv, 1)[[1]]
lBestjuv = update(lBestjuv, REML =T)
summary(lBestjuv)



################################################
#acartiella adulat
acad = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + lilzoops+ logChlsc +(1|Year), 
               zi = ~ log(SalSurf), data = Acart_ad2,
               family= "nbinom2", na.action = "na.fail")

summary(acad)
acdad = dredge(acad, trace = 2, rank = "BIC",  
              extra = list("Rs" = function(x) {
                s <- performance(x)
                c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
              }))
acBestad = get.models(acdad, 1)[[1]]
acBestad = update(acBestad, REML =T)
summary(acBestad)

#Juvenile acartiella
Acart_juv = filter(Acart_juv, ! is.na(lilzoops))
acjuv = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + lilzoops+ logChlsc +(1|Year), 
                zi = ~ log(SalSurf), data = Acart_juv,
                family= "nbinom2", na.action = "na.fail")

summary(acjuv)
acdjuv = dredge(acjuv, trace = 2, rank = "BIC",  
               extra = list("Rs" = function(x) {
                 s <- performance(x)
                 c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
               }))
acBestjuv = get.models(acdjuv, 1)[[1]]
acBestjuv = update(acBestjuv, REML =T)
summary(acBestjuv)





#################################################################
#


#########################################################################

save(ad, pd2, ed, acd, bd,dd,  dh, ds, kd, ld, aBest,dBest, pBest2, eBest, hBest2,  acBest, tBest, bBest, sBest, kbest, lBest, file = "outputs/SpeciesModels_21Jan2024.RData")
#load("outputs/SpeciesModels_20Jan2024.RData")

#################################################################################
#OK, now I need a way of displaying all these predictors and parameters and stuff in a 
#way that makes sense. I'm going to have to think about this. 


getcoefs =  function(x, species) {
  dat = as.data.frame(summary(x)$coefficients[[1]])
  dat =   mutate(dat, params = row.names(dat), Species = species, paramtype = "count")
  if(!is_null(x[[3]])) { 
    dat2 = as.data.frame(summary(x)$coefficients[[2]])
    dat2 = mutate(dat2, params = row.names(dat2), Species = species, paramtype = "zi")
    dat = bind_rows(dat, dat2)
  }
  return(dat)
}


rsqs =  function(y, species) {
  R2 = r2(y)[[2]]
  dat = data.frame(Species = species, R2a = R2)
  
  return(dat)
}


species = "Pseudodiaptomus"

agrigate = bind_rows(getcoefs(pBestad, "Pseudodiaptomus adult"),
                     getcoefs(pBestjuv, "Pseudodiaptomus juvenikle"),
                     getcoefs(eBestad, "Eurytemor adult"),
                     getcoefs(eBestjuv, "Eurytemor juvenile"),
                     getcoefs(tBestad, "Tortanus adult"),
                     getcoefs(tBestjuv, "Tortanus juvenikle"),
                     getcoefs(lBestjuv, "Limnoithona juvenile"),
                     getcoefs(lBestad, "Limnoithona adult"),
                     getcoefs(acBestad, "Acatrtilla adult"),
                     getcoefs(acBestjuv, "Acartiella Juvinile"))%>%
  mutate(pp = paste(paramtype, params))

test = agrigate  %>%
  pivot_wider(id_cols = Species, names_from = pp, values_from = Estimate) %>%
  pivot_longer(cols = 2:last_col(), names_to = "pp", values_to = "Estimate")


test2 = left_join(test, agrigate)

arg= str_split_fixed(test2$pp, " ", n = 2)
test3 = cbind(test2, arg) %>%
  rename(Paramtype = `1`, Param = `2`) %>%
  mutate(params = NULL, paramtype = NULL, Estimate2 = case_when(
    `Pr(>|z|)` < 0.05 & Estimate < 0 ~ "Decrease",
    `Pr(>|z|)` < 0.05 & Estimate >0 ~ "Increase",
    `Pr(>|z|)` > 0.05 ~ "Not Significant",
    TRUE ~ "Not Included"
  ))

ggplot(filter(test3, Param != "(Intercept)"), aes(x = Species, y = Param, fill = Estimate)) + geom_tile() +
  facet_wrap(~Paramtype) +
  scale_fill_gradient2()


test3 = mutate(test3, ParamX = factor(Param, levels = c("(Intercept)", "log(SalSurf)", "Secchisc",
                                                        "logChlsc", "lilzoops"),
                                      labels = c("Intercept", "Salinity", "Secchi Depth", 
                                                  "Chlorophyll", "Microzooplankton")),
               Paramtype = factor(Paramtype, levels = c("count", "zi"), labels =
                                    c("Count", "Zero Inflation")))


ggplot(filter(test3, Paramtype == "Count"), aes(x = ParamX, y = Species, fill = Estimate2)) + geom_tile() +
 # facet_wrap(~Paramtype) +
  scale_fill_manual(values = c("lightblue", "red", "grey", "white"), name = "Model \nEstimate") +
  ylab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_text(aes(label = round(Estimate, dig = 2)))

#Maybe it makes more sense to plot. the estimates on the response scale
test3 = mutate(test3, Estimate3 = case_when(
  `Pr(>|z|)` < 0.05 ~ exp(Estimate)
))
ggplot(filter(test3, Param != "Month", Param != "(Intercept)"), 
       aes(x = Param, y = Species, fill = Estimate3)) + geom_tile() +
  facet_wrap(~Paramtype) +
  scale_fill_gradient2()


#Get size classes and groups to sort by
library(readxl)
taxa = read_excel("Data/taxa.xlsx")

test3 = left_join(test3, taxa, by = c("Species"= "Taxon")) %>%
  mutate(Estimate2 = case_when(ParamX == "Microzooplankton" & Species %in% c("Kerottela","Syncheata") ~ "Not Tested",
                               TRUE ~ Estimate2))

ggplot(filter(test3, !is.na(Paramtype)), 
       aes(x = ParamX, y = Species, fill = Estimate2)) + geom_tile() +
   facet_grid(~Paramtype, scales = "free_x", space = "free") +
  scale_fill_manual(values = c("lightblue", "red", "grey", "white", "grey20"), name = "Model \nEstimate") +
  ylab(NULL) + xlab(NULL)+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  geom_text(aes(label = round(Estimate, dig = 2)))+
  facet_grid(FeedingGuild~Paramtype, scales = "free", space = "free")


#Think about how best to display this. May not be fair to put everyone on the same scale IDK
test3 = mutate(test3, sig = case_when(
  `Pr(>|z|)` < 0.05 ~ "p <0.05",
  `Pr(>|z|)` > 0.05 ~ "NS"
))

#Point plot with standard error
ggplot(filter(test3,Paramtype == "Count", ParamX != "Intercept"), 
       aes(x = Species, y = Estimate, color = ParamX)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Estimate- `Std. Error`, ymax = Estimate + `Std. Error`))+
  facet_wrap(~ParamX, scales = "free_y")+
  scale_alpha_manual(values = c(0.3, 1))+ 
  theme(axis.text.x = element_text(angle = 90))


#Bar plot
ggplot(filter(test3, ParamX != "Intercept"), aes(x = Species, y = Estimate, fill = ParamX, alpha = sig)) + 
  geom_bar(stat = "identity") + facet_grid(ParamX~Paramtype, scales = "free_y")+
  scale_alpha_manual(values = c(0.3, 1)) 

plot(allEffects(eBest))

#######################################################

#can I put all the effects plot together? Would that be aweful?

effectsfun = function(mod, taxon) {
  effs = allEffects(mod)
  df = data.frame(NA) 
  for(i in 1:length(effs)){
    df2 = data.frame(effs[[i]]) %>%
      mutate(param = names(effs)[i])
    names(df2) = c("level", "fit", "se", "lower", "upper", "param")
    df = bind_rows(df, df2)
  }
  df$taxa = taxon
return(df)
}

effectsall = bind_rows(effectsfun(pBest2, "Pseudodiaptomus"),
                     effectsfun(eBest, "Eurytemora"),
                     effectsfun(tBest, "Tortanus"),
                     effectsfun(bBest, "Bosmina"),
                     effectsfun(kbest, "Kerottela"),
                     effectsfun(lBest, "Limnoithona"),
                     effectsfun(hBest2, "H. Logirostris"),
                     effectsfun(aBest, "Acartiella"),
                     effectsfun(sBest, "Syncheata"),
                     effectsfun(acBest, "Acanthocyclops")) %>%
  dplyr::select(-NA.) %>%
  filter(!is.na(param)) %>%
  mutate(param = case_when(param == "I(Month^2)" ~ "Month",
                           TRUE ~ param))

#something is weird with the x-axis on salsurf
ggplot(effectsall, aes(x = level, y = fit))+ 
  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+
  geom_smooth()+
  facet_grid(taxa~param, scales = "free")

ggplot(effectsall, aes(x = level, y = fit))+ 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+
  geom_line()+
  facet_grid(taxa~param, scales = "free")+
  ylab("Partial Residuals")+
  xlab("Value of Predictor variable")

##########################################################################################

#look at correlations between predictor variables
WQandLilzoops = left_join(WQ_dataNARM, lilzoopsmean, by = c("Month", "Year", "Region")) %>%
  filter(Year >1999)
mat = dplyr::select(WQandLilzoops, Temperature, Secchi, Chlorophyll, Conductivity, lilzoopbiomass) %>%
  dplyr::filter(!is.na(Temperature), !is.na(Secchi), !is.na(Chlorophyll), !is.na(Conductivity), !is.na(lilzoopbiomass))
wqcorr = cor(as.matrix(mat))

corrplot.mixed(wqcorr, lower.col = c("blue", "red"))
