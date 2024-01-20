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

load("data/ZoopsSum.RData")

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
Filters = function(taxname, data = ZoopsSum){
  data = filter(data, Taxname2 == taxname) #filter to taxon of interest
  
  #summarize by month and select the highest 5 months
  monthsdat = data %>%
    group_by(Month, Taxname2) %>%
    summarize(CPUE = mean(CPUE))%>%
    arrange(desc(CPUE)) %>%
    ungroup() %>%
    slice(1:5)
  
  months = monthsdat$Month
  
  #find the salinities with the most catch
  test = dplyr::filter(data,  CPUE !=0)
   limsSurf = quantile(test$SalSurf, probs = c(0.05, 0.95), na.rm =T)
   
   #filter data to months and salinities with highest catch, take mean by region and month, and summarize
  datafiltered = dplyr::filter(data, Month %in% months, 
                 between(SalSurf, limsSurf[1], limsSurf[2])) %>%
    group_by(Month, Region, Season, Year) %>%
    summarize(CPUE = mean(CPUE)) %>%
    left_join(WQ_dataNARM, by = c("Month", "Year", "Region")) %>%
    mutate(SalSurf = ec2pss(Conductivity/1000, t=25),
           logChl = log(Chlorophyll),
           logCPUE = log(CPUE +1),
           rCPUE = round(CPUE),
           Secchisc = scale(Secchi),
           logChlsc = scale(logChl),
           SalSurfsc = scale(SalSurf)) 
}

#start with p. forbesi adults. just the appropraite salinity range, most abundant months
Pseudox = Filters("Pseudodiaptomus forbesi") %>%
  filter(!is.na(Secchi), !is.na(logChl), !is.na(SalSurf))


#quick exploritory plots
hist(Pseudox$CPUE)
hist(log(Pseudox$CPUE +1))


ggplot(Pseudox, aes(x = Year, y = logCPUE)) + geom_point() 
ggplot(Pseudox, aes(x = as.factor(Month), y = logCPUE)) + geom_boxplot() + facet_wrap(~Region)
ggplot(Pseudox, aes(x = SalSurf, y = logCPUE)) + geom_point() + facet_wrap(~Region)

#initial model
m2 = lmer(logCPUE ~ log(SalSurf)+logChlsc + Month + I(Month^2)+ (1|Year), data = Pseudox)
summary(m2)
plot(m2)
plot(allEffects(m2))
plot(simulateResiduals(m2))
testZeroInflation(m2)

#zero inflated model (better)
p5 = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc +  Month + I(Month^2)+ (1|Year),
             ziformula = ~log(SalSurf),
             data= Pseudox, family = nbinom2, na.action = "na.fail")
summary(p5)
plot(simulateResiduals(p5))
testZeroInflation(p5)
plot(allEffects(m2))

#do we really need all those terms in the model?
#I was going to do something with flow...
#but I already have year as a random effect
pd = dredge(p5, rank = "BIC", trace = 2,extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
pBest = get.models(pd, 1)[[1]]
summary(pBest)

p2 = glmmTMB(rCPUE ~log(SalSurf)+ Secchisc + logChlsc+ Month + I(Month^2) + (1|Region) + (1|Year),  
             ziformula = ~log(SalSurf),
             data= Pseudox, family = "nbinom2", na.action = "na.fail")
summary(p2)

pd2 = dredge(p2, trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))

pBest2 = get.models(pd2, 1)[[1]]

pBest2 = update(pBest2, REML = TRUE)
summary(pBest2)
plot(simulateResiduals(pBest2))
testZeroInflation(pBest2)
r2(pBest2)
performance(pBest2)
icc(pBest2)
plot(allEffects(pBest2))

######################################################


#Eurytemora
Eury = Filters("Eurytemora affinis", data = ZoopsSum)%>%
  filter(!is.na(Secchi))

hist(Eury$CPUE)
hist(log(Eury$CPUE +1))
#hm. That one is more zero-inflated. May be trouble.


e5 = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + Secchisc + Month + I(Month^2) +  (1|Region) + (1|Year), 
             zi = ~ log(SalSurf) , data = Eury, 
             family= "nbinom2", na.action = "na.fail")


summary(e5)
ed = dredge(e5,  trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
eBest = get.models(ed, 1)[[1]]
eBest = update(eBest, REML = T)
summary(eBest)
plot(simulateResiduals(eBest))

plot(allEffects(eBest))
#Yeah, I think I like this one best.
#check the diagnostic plots

performance(eBest)
plot(allEffects(eBest))
#pesudos can eat everything, so no effect of chlorophyll, eurys much more narrow diet, 

########################################################################

#Acanthocyclops
Acan = Filters("Acanthocyclops_UnID", data = ZoopsSum)%>%
  filter(!is.na(Secchi))

hist(Acan$CPUE)
hist(log(Acan$CPUE +1))
#bleh. Looks promlematic


ac5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
              zi = ~log(SalSurf), 
              data = Acan, 
              family= "nbinom2", na.action = "na.fail")
summary(ac5)
acd = dredge(ac5,  trace = 2, rank = "BIC",  extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
acBest = get.models(acd, 1)[[1]]
acBest = update(acBest, REML = TRUE)
summary(acBest)
plot(allEffects(acBest))
performance(acBest)

#check the diagnostic plots

simresac5 = simulateResiduals(acBest)
plot(simresac5)



#Acartiella
#Acartiella copepodids weren't seperated until 2006, so I'm going to subset
#to just that time period.
Acart = Filters( "Acartiella_UnID" , data = ZoopsSum)%>%
  filter(!is.na(Secchi), Year >2005)

hist(Acart$CPUE)
hist(log(Acart$CPUE +1))
ggplot(Eury, aes(x = Year, y = logCPUE)) + geom_point() + facet_wrap(~Region)
ggplot(Eury, aes(x = as.factor(Month), y = logCPUE)) + geom_boxplot() + facet_wrap(~Region)
ggplot(Acart, aes(x = SalSurfsc, y = logCPUE)) + geom_point() + facet_wrap(~Region)
ggplot(Acart, aes(x = log(SalSurf), y = logCPUE)) + geom_point() + facet_wrap(~Region)

#log salinity works way better... should I try that for everythign?
#Hmm. Maybe we are having problems because we almost never catch them in the North Delta
#Acart = droplevels(filter(Acart, Region != "North"))

a5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc+(1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Acart, 
             family= "nbinom2", na.action = "na.fail")
summary(a5)
ad = dredge(a5,trace = 2, rank = "BIC", extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
aBest = get.models(ad, 1)[[1]]
aBest = update(aBest, REML =T)
summary(aBest)

simresa5 = simulateResiduals(aBest)
plot(simresa5)
#gross! WTF?
#maybe we just don' tcatch enough of them or something?
performance(aBest)

###################################################################
#Bosmina
Bos = Filters("Bosmina longirostris", data = ZoopsSum)%>%
  filter(!is.na(Secchi))

hist(Bos$CPUE)
hist(log(Bos$CPUE +1))
#hm.


b5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc +(1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Bos,
             family= "nbinom2", na.action = "na.fail")

summary(b5)
bd = dredge(b5, trace = 2,  rank = "BIC",
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
bBest = get.models(bd, 1)[[1]]
bBest = update(bBest, REML =T)
summary(bBest)

simresb5 = simulateResiduals(b5)
plot(simresb5)
performance(bBest)
plot(allEffects(bBest))



###################################################################
#Tortanus
Tort = Filters("Tortanus_UnID", data = ZoopsSum)%>%
  filter(!is.na(Secchi), !is.nan(SalSurfsc))

hist(Tort$CPUE)
hist(log(Tort$CPUE +1))
#hm.


t5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc +(1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Tort,
             family= "nbinom2", na.action = "na.fail")

summary(t5)
td = dredge(t5, trace = 2, rank = "BIC",  
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
tBest = get.models(td, 1)[[1]]
tBest = update(tBest, REML = T)
summary(tBest)

simrest5 = simulateResiduals(t5)
plot(simrest5)
performance(tBest)
plot(allEffects(tBest))

################################################
load("data/ZoopsSumm.RData")
Limno = Filters("Limnoithona", data = ZoopsSumm)%>%
  filter(!is.na(Secchi))


hist(Limno$CPUE)
hist(log(Limno$CPUE +1))
#hm.


l5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc +(1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Tort,
             family= "nbinom2", na.action = "na.fail")

summary(l5)
ld = dredge(l5, trace = 2, rank = "BIC",  
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
lBest = get.models(ld, 1)[[1]]
lBest = update(lBest, REML =T)
summary(lBest)

simrest5 = simulateResiduals(l5)
plot(simrest5)
performance(lBest)
plot(allEffects(lBest))


################################################
Asplanchna = Filters("Asplanchna", data = ZoopsSumm)%>%
  filter(!is.na(Secchi))


hist(Asplanchna$CPUE)
hist(log(Asplanchna$CPUE +1))
#This one might be rought


as5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc +(1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Asplanchna,
             family= "nbinom2", na.action = "na.fail")

summary(as5)
asd = dredge(as5, trace = 2, rank = "BIC", 
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
asBest = get.models(asd, 1)[[1]]
asBest = update(asBest, REML = T)
#might not have enough catch of these guys to model it.

summary(asBest)

simrest5 = simulateResiduals(as5)
plot(simrest5)
performance(asBest)
plot(allEffects(asBest))

#####################################################################
#Some random rotifer

Kerat = Filters("Keratella", data = ZoopsSumm) %>%
  filter(!is.na(Secchisc))
hist(Kerat$CPUE)
hist(log(Kerat$CPUE +1))
#perfectly bell-shaped except for the zeros

k5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc+ Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Kerat, 
             family= "nbinom2", na.action = "na.fail")

summary(k5)
kd = dredge(k5, trace = 2, rank = "BIC",
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
kbest = get.models(kd, 1)[[1]]
kbest = update(kbest, REML =T)
summary(kbest)
simresk5 = simulateResiduals(kbest)
plot(simresk5)



#another rotifer

Sync = Filters("Synchaeta", data = ZoopsSumm) %>%
  filter(!is.na(Secchisc))
hist(Sync$CPUE)
hist(log(Sync$CPUE +1))
#perfectly bell-shaped except for the zeros

s5 = glmmTMB(rCPUE ~log(SalSurf) +Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ log(SalSurf), data = Sync, 
             family= "nbinom2", na.action = "na.fail")

summary(s5)
ds = dredge(s5, trace = 2, rank = "BIC", 
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
sBest = get.models(ds, 1)[[1]]
sBest = update(sBest, REML =T)

summary(sBest)
simresd = simulateResiduals(sBest)
plot(simresd)
performance(sBest)


#################################################################
#Hyperacanthomysis logirostris
load("data/ZoopsSumMM2.RData")
Hyp2 = Filters("Hyperacanthomysis longirostris", data = ZoopsSumMM2) %>%
  filter(!is.na(Secchisc))

hist(Hyp2$CPUE)
hist(log(Hyp2$CPUE+1))
# a little flat and zero-inflated, bu tnot bad
h52 = glmmTMB(rCPUE ~log(SalSurf) +Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
              zi = ~ log(SalSurf), data = Hyp2, 
              family= "nbinom2", na.action = "na.fail")

summary(h52)
dh = dredge(h52, trace = 2, rank = "BIC",  extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))

hBest2 = get.models(dh, 1)[[1]]
hBest2 = update(hBest2, REML =T)

summary(hBest2)
simresd = simulateResiduals(hBest2)
plot(simresd)
check_collinearity(hBest2)
performance(hBest2)



#########################################################################

save(ad, pd2, ed, acd, bd, dh, ds, kd, ld, aBest, pBest2, eBest, hBest2,  acBest, tBest, bBest, sBest, kbest, lBest, file = "SpeciesModels_2Jan2024.RData")
load("outputs/SpeciesModels_2Jan2024.RData")

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

agrigate = bind_rows(getcoefs(pBest2, "Pseudodiaptomus"),
                     getcoefs(eBest, "Eurytemora"),
                     getcoefs(tBest, "Tortanus"),
                     getcoefs(bBest, "Bosmina"),
                     getcoefs(kbest, "Kerottela"),
                     getcoefs(lBest, "Limnoithona"),
                     getcoefs(hBest2, "H. Logirostris"),
                     getcoefs(aBest, "Acartiella"),
                     getcoefs(sBest, "Syncheata"),
                     getcoefs(acBest, "Acanthocyclops"))%>%
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
foo = data.frame(Species = unique(test2$Species), Paramtype = "Zero Inflation", 
                 ParamX = "Chlorophyll", Estimate2 = "Not Included")
test3 = bind_rows(test3, foo)

test3 = mutate(test3, ParamX = factor(Param, levels = c("(Intercept)", "log(SalSurf)", "Secchisc",
                                                        "Month", "I(Month^2)", "logChlsc"),
                                      labels = c("Intercept", "Salinity", "Secchi Depth", 
                                                 "Month", "Month^2", "Chlorophyll")),
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


ggplot(filter(test3, Paramtype == "Count", Param != "Month", Param != "I(Month^2)", Param != "(Intercept)"), 
       aes(x = ParamX, y = Species, fill = Estimate2)) + geom_tile() +
  # facet_wrap(~Paramtype) +
  scale_fill_manual(values = c("lightblue", "red", "grey", "white"), name = "Model \nEstimate") +
  ylab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_text(aes(label = round(Estimate3, dig = 2)))


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

ggplot(effectsall, aes(x = level, y = fit))+ 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+
  geom_smooth()+
  facet_grid(taxa~param, scales = "free")

ggplot(effectsall, aes(x = level, y = fit))+ 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+
  geom_line()+
  facet_grid(taxa~param, scales = "free")+
  ylab("Partial Residuals")+
  xlab("Value of Predictor variable")

