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
source("HelperFuncitons.R")


#Note: i tried a lag of chlorophyll a earlier and it didn't work very well
#start with Pseudodiaptomus

Regions<-read_csv("data/Rosies_regions.csv")

Delta<-deltamapr::R_EDSM_Subregions_Mahardja_FLOAT%>%
  filter(SubRegion%in%unique(Regions$SubRegion))%>%  #Filter to regions of interest
  dplyr::select(SubRegion)

Regs = unique(Regions[,c(1,5)])
Delta = merge(Delta, Regs)


Pseudox_juv = Filters_LS("Pseudodiaptomus_UnID", "Juvenile", ZoopsSumLifestage) 

Pseudox_ad = Filters_LS("Pseudodiaptomus forbesi", "Adult", ZoopsSumLifestage)

Eury_juv = Filters_LS("Eurytemora affinis", "Juvenile", data =ZoopsSumLifestage)

Eury_ad = Filters_LS("Eurytemora affinis", "Adult", data = ZoopsSumLifestage)
  
Acan_juv = Filters_LS("Acanthocyclops_UnID", "Juvenile", data = ZoopsSumLifestage)
  
#Ok, we don't do juvenile acanthocyclops, so never mind

Acan_ad = Filters_LS("Acanthocyclops_UnID", "Adult", data = ZoopsSumLifestage)
 

Limno_juv = Filters_LS("Limnoithona", "Juvenile", data = ZoopsSumm_lifestage)

Limno_ad = Filters_LS("Limnoithona", "Adult", data = ZoopsSumm_lifestage)


Tort_juv = Filters_LS( "Tortanus_UnID" , "Juvenile", data = ZoopsSumLifestage)

Tort_ad = Filters_LS( "Tortanus_UnID" , "Adult", data = ZoopsSumLifestage)

#acartiella juveniles weren't speticated until after 2005
Acart_juv = Filters_LS( "Acartiella_UnID" , "Juvenile", data = ZoopsSumLifestage)%>%
  filter( Year >2005)
Acart_ad = Filters_LS( "Acartiella_UnID" , "Adult", data = ZoopsSumLifestage) %>%
  filter(Year >2005)


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
p5adx = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc +  (1|Month)+ (1|Year),
             ziformula = ~log(SalSurf),
             data= Pseudox_ad, family = nbinom2, na.action = "na.fail")
summary(p5adx)


#go through all the modesl
pdadx = dredge(p5adx, rank = "BIC", trace = 2,extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
pBestadx = get.models(pdadx, 1)[[1]]
summary(pBestadx)

pBestadx = update(pBestadx, REML = TRUE)
summary(pBestadx)




#now the juvenikels
p5juvx = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc +  (1|Month)+ (1|Year),
               ziformula = ~log(SalSurf),
               data= Pseudox_juv, family = nbinom2, na.action = "na.fail")
summary(p5juvx)



pdjuvx = dredge(p5juvx, rank = "BIC", trace = 2,extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
pBestjuvx = get.models(pdjuvx, 1)[[1]]
summary(pBestjuvx)

pBestjuvx = update(pBestjuvx, REML = TRUE)
summary(pBestjuvx)
plot(simulateResiduals(pBestjuvx))


hist(Pseudox_ad$logCPUE)
hist(Pseudox_juv$logCPUE)

#hmmm.
######################################################


#Eurytemora

#adults
e5adx = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + Secchisc + (1|Month) + (1|Year), 
             zi = ~ log(SalSurf) , data = Eury_ad, 
             family= "nbinom2", na.action = "na.fail")


summary(e5adx)
plot(allEffects(e5adx))

summary(e5adx)
edadx = dredge(e5adx,  trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
eBestadx = get.models(edadx, 1)[[1]]
eBestadx = update(eBestadx, REML = T)
summary(eBestadx)
plot(simulateResiduals(eBestadx))
#Huh. 


#juvenikels

e5juvx = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + Secchisc + (1|Month) + (1|Year), 
               zi = ~ log(SalSurf) , data = Eury_juv, 
               family= "nbinom2", na.action = "na.fail")

summary(e5juvx)

plot(allEffects(e5juvx))

summary(e5juvx)
edjuvx = dredge(e5juvx,  trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
eBestjuvx = get.models(edjuvx, 1)[[1]]
eBestjuvx = update(eBestjuvx, REML = T)
summary(eBestjuvx)
plot(simulateResiduals(eBestjuvx))

plot(allEffects(eBestjuvx))

Euryall = bind_rows(Eury_juv, Eury_ad)

ggplot(Euryall, aes(x = Month, y = log(CPUE), color = Lifestage))+
  geom_point()+
  facet_grid(Year~Region)

#hm. Something is odd in the north in recent years
#######################################################################


#Acartiella
#Acartiella copepodids weren't seperated until 2006, so I'm going to subset
#to just that time period.

a5adx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + logChlsc+ + (1|Year), 
             zi = ~ log(SalSurf), data = Acart_ad, 
             family= "nbinom2", na.action = "na.fail")
summary(a5adx)
adadx = dredge(a5adx,trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
aBestadx = get.models(adadx, 1)[[1]]
aBestadx = update(aBestadx, REML =T)
summary(aBestadx)

a5juvx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + logChlsc+ + (1|Year), 
                zi = ~ log(SalSurf), data = Acart_juv, 
                family= "nbinom2", na.action = "na.fail")
summary(a5juvx)
adjuvx = dredge(a5juvx,trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
aBestjuvx = get.models(adjuvx, 1)[[1]]
aBestjuvx = update(aBestjuvx, REML =T)
summary(aBestjuvx)

###################################################################
#Tortanus
t5adx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc +  (1|Month) + logChlsc + (1|Year), 
             zi = ~ log(SalSurf), data = Tort_ad,
             family= "nbinom2", na.action = "na.fail")

summary(t5adx)
tdadx = dredge(t5adx, trace = 2, rank = "BIC",  
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
tBestadx = get.models(tdadx, 1)[[1]]
tBestadx = update(tBestadx, REML = T)
summary(tBestadx)

#Tortatuns juveniles
t5juvx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc  + (1|Month) + logChlsc + (1|Year), 
               zi = ~ log(SalSurf), data = Tort_juv,
               family= "nbinom2", na.action = "na.fail")

summary(t5juvx)
tdjuvx = dredge(t5juvx, trace = 2, rank = "BIC",  
              extra = list("Rs" = function(x) {
                s <- performance(x)
                c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
              }))
tBestjuvx = get.models(tdjuvx, 1)[[1]]
tBestjuvx = update(tBestjuvx, REML = T)
summary(tBestjuvx)
#whereas juveniles like chlorophyll

################################################

#limnoithona adulat
l5adx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + logChlsc +(1|Year), 
             zi = ~ log(SalSurf), data = Limno_ad,
             family= "nbinom2", na.action = "na.fail")

summary(l5adx)
ldadx = dredge(l5adx, trace = 2, rank = "BIC",  
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
lBestadx = get.models(ldadx, 1)[[1]]
lBestadx = update(lBestadx, REML =T)
summary(lBestadx)

#Juvenile limnos
l5juvx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + logChlsc +(1|Year), 
               zi = ~ log(SalSurf), data = Limno_juv,
               family= "nbinom2", na.action = "na.fail")

summary(l5juvx)
ldjuvx = dredge(l5juvx, trace = 2, rank = "BIC",  
              extra = list("Rs" = function(x) {
                s <- performance(x)
                c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
              }))
lBestjuvx = get.models(ldjuvx, 1)[[1]]
lBestjuvx = update(lBestjuvx, REML =T)
summary(lBestjuvx)



################################################
#acartiella adulat
acadx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + logChlsc +(1|Year), 
               zi = ~ log(SalSurf), data = Acart_ad,
               family= "nbinom2", na.action = "na.fail")

summary(acadx)
acdadx = dredge(acadx, trace = 2, rank = "BIC",  
              extra = list("Rs" = function(x) {
                s <- performance(x)
                c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
              }))
acBestadx = get.models(acdadx, 1)[[1]]
acBestadx = update(acBestadx, REML =T)
summary(acBestadx)

#Juvenile acartiella

acjuvx = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + (1|Month) + logChlsc +(1|Year), 
                zi = ~ log(SalSurf), data = Acart_juv,
                family= "nbinom2", na.action = "na.fail")

summary(acjuvx)
acdjuvx = dredge(acjuvx, trace = 2, rank = "BIC",  
               extra = list("Rs" = function(x) {
                 s <- performance(x)
                 c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
               }))
acBestjuvx = get.models(acdjuvx, 1)[[1]]
acBestjuvx = update(acBestjuvx, REML =T)
summary(acBestjuvx)





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

agrigate = bind_rows(getcoefs(pBestadx, "Pseudodiaptomus adult"),
                     getcoefs(pBestjuvx, "Pseudodiaptomus juvenile"),
                     getcoefs(eBestadx, "Eurytemora adult"),
                     getcoefs(eBestjuvx, "Eurytemora juvenile"),
                     getcoefs(tBestadx, "Tortanus adult"),
                     getcoefs(tBestjuvx, "Tortanus juvenile"),
                     getcoefs(lBestjuvx, "Limnoithona juvenile"),
                     getcoefs(lBestadx, "Limnoithona adult"),
                     getcoefs(acBestadx, "Acartiella adult"),
                     getcoefs(acBestjuvx, "Acartiella juvenile"))%>%
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
