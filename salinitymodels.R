
library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)
library(Hmsc)
library(glmmTMB)
library(DHARMa)
library(vegan)
library(performance)
library(car)

load("ZoopMaster.RData")
load("ZoopMasterMicro.RData")
load("ZoopMasterMysids.RData")

#Definite big pile of NOPE on STN and FMWT data.
#Also nope on lag of CHLA

#Scale predictor variables, do some other stuff, log-transormations
ZoopMaster2 = mutate(ZoopMaster, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                     SalSurfsc = scale(SalSurf), logChlsc = scale(logChl),
                     lagChlsc = scale(logChlag),
                     Region = as.factor(Region),
                     Secchisc = scale(Secchi),
                     Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year))
#Scale predictor variables, do some other stuff, log-transormations
ZoopMaster2m = mutate(ZoopMasterm, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                     SalSurfsc = scale(SalSurf), logChlsc = scale(logChl),
                     Region = as.factor(Region),
                     Secchisc = scale(Secchi),
                     Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year))

#Scale predictor variables, do some other stuff, log-transormations
ZoopMaster2MM = mutate(ZoopMasterMM, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                      SalSurfsc = scale(SalSurf), logChlsc = scale(logChl),
                      Region = as.factor(Region),
                      Secchisc = scale(Secchi),
                      Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year))

ZoopMaster2MM2 = mutate(ZoopMasterMM2, rCPUE = round(CPUE), logCPUE = log(CPUE+1),
                       SalSurfsc = scale(SalSurf), logChlsc = scale(logChl),
                       Region = as.factor(Region),
                       Secchisc = scale(Secchi),
                       Chlorophyllsc = scale(Chlorophyll), Yearsc = scale(Year))


#Maybe filter by time of year and salinity range things are most abundant
Filters = function(taxname, data = EMPallc){
 # test = dplyr::filter(data, Taxname2 == taxname, CPUE !=0)
 # limsSurf = quantile(test$SalSurf, probs = c(0.05, 0.95))
  #limsMon = quantile(test$Month, probs = c(0.05, 0.95))
  limsYear = 1996 
  dplyr::filter(data,  Taxname2 == taxname,#between(Month, limsMon[1], limsMon[2]), 
               # between(SalSurf, limsSurf[1], limsSurf[2]),
                Year >= limsYear)
}

#Meh, lag of chlorophyll not looking so good

#now subset just pseudodiaptomus to play with
Pseudos = Filters("Pseudodiaptomus", data = ZoopMaster2) %>%
  filter(!is.na(Secchi))
hist(Pseudos$CPUE)
hist(log(Pseudos$CPUE +1))


ggplot(Pseudos, aes(x = Year, y = logCPUE)) + geom_point() + facet_wrap(~Region)
ggplot(Pseudos, aes(x = as.factor(Month), y = logCPUE)) + geom_boxplot() + facet_wrap(~Region)
ggplot(Pseudos, aes(x = SalSurf, y = logCPUE)) + geom_point() + facet_wrap(~Region)

p5 = glmmTMB(rCPUE ~log(SalSurf)+ Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year),  
             data= Pseudos, family = "nbinom2", na.action = "na.fail")

p5t = glmmTMB(rCPUE ~log(SalSurf)+ Secchisc + Month + I(Month^2) + logChlsc,  
             data= Pseudos, family = "nbinom2", na.action = "na.fail")

zinb(rCPUE ~ SalSurf + Month, Pseudos)

p5x = glmmTMB(rCPUE ~log(SalSurf)+ Secchisc + lagChlsc+ Month + I(Month^2) + (1|Region) + (1|Year),  
             data= Pseudos, family = "nbinom2", na.action = "na.fail")
summary(p5)
summary(p5x)
#OK, so the lag of chlorophyll really isn't helping any.


pd = dredge(p5, trace = 2,extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
pBest = get.models(pd, 1)[[1]]
p2 = glmmTMB(rCPUE ~log(SalSurf)+ Secchisc + lagChlsc+ Month + I(Month^2) + (1|Region) + (1|Year),  
             data= Pseudos, family = "nbinom2", na.action = "na.fail")
summary(p2)
foo = summary(pBest)
Anova(pBest, test.statistic = "F", compontent = "cond")
#all the terms are in the top model

#check the diagnostic plots
simresp5 = simulateResiduals(pBest)
plot(simresp5)
r2(pBest)
performance(pBest)
icc(pBest)

#Eurytemora
Eury = Filters("Eurytemora affinis", data = ZoopMaster2)%>%
  filter(!is.na(Secchi))

hist(Eury$CPUE)
hist(log(Eury$CPUE +1))
#hm. That one is more zero-inflated. May be trouble.


e5 = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + Secchisc + Month + I(Month^2) +  (1|Region) + (1|Year), 
             zi = ~ log(SalSurf)+ logChlsc + Secchisc + Month + I(Month^2) + (1|Region) + (1|Year), data = Eury, 
             family= "nbinom2", na.action = "na.fail")


e5t = glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + Secchisc + Month + I(Month^2), 
             zi = ~ log(SalSurf)+ logChlsc + Secchisc + Month + I(Month^2), data = Eury, 
             family= "nbinom2", na.action = "na.fail")


summary(e5t)
ed = dredge(e5,  trace = 2, extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
eBest = get.models(ed, 1)[[1]]
eurybest =  glmmTMB(rCPUE ~log(SalSurf)+ logChlsc + Month + I(Month^2) +  (1|Region) + (1|Year), 
                    zi = ~  Secchisc + Month + I(Month^2) + (1|Region) + (1|Year), 
                     data = Eury, 
                    family= "nbinom2", na.action = "na.fail")
#hmmm, Chl

summary(eurybest)
#WTF

#Yeah, I think I like this one best.
#check the diagnostic plots
simrese5 = simulateResiduals(eBest)
plot(simrese5)
performance(eBest)

#Acanthocyclops
Acan = Filters("Acanthocyclops", data = ZoopMaster2)%>%
  filter(!is.na(Secchi))

hist(Acan$CPUE)
hist(log(Acan$CPUE +1))
#bleh. Looks promlematic


ac5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~log(SalSurf) + Secchisc + Month + I(Month^2) +  (1|Region)+ logChlsc + (1|Year), 
             data = Acan, 
             family= "nbinom2", na.action = "na.fail")
summary(ac5)
acd = dredge(ac5,  trace = 2, fixed = c("cond(log(SalSurf))","zi(Month)"), extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
acBest = get.models(acd, 1)[[1]]
summary(acBest)
performance(acBest)

#check the diagnostic plots
simresac5 = simulateResiduals(acBest)
plot(simresac5)



#Acartiella
#Acartiella copepodids weren't seperated until 2006, so I'm going to subset
#to just that time period.
Acart = Filters("Acartiella", data = ZoopMaster2)%>%
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
             zi = ~ log(SalSurf)+ logChlsc+ Secchisc +Month + I(Month^2) + (1|Region) + (1|Year), data = Acart, 
             family= "nbinom2", na.action = "na.fail")
summary(a5)
ad = dredge(a5,trace = 2, fixed = c("cond(Month)", "cond(log(SalSurf))", "zi(log(SalSurf))"), extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
aBest = get.models(ad, 1)[[1]]
summary(aBest)
a5x = glmmTMB(rCPUE ~log(SalSurf) + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ log(SalSurf) + Month +  (1|Region) + (1|Year), data = Acart, 
             family= "nbinom2", na.action = "na.fail")
summary(a5x)

summary(aBest)
#Those zero inflation coefficients are GROSS
simresa5 = simulateResiduals(aBest)
plot(simresa5)
performance(aBest)

#Bosmina
Bos = Filters("Bosmina longirostris", data = ZoopMaster2)%>%
  filter(!is.na(Secchi))

hist(Bos$CPUE)
hist(log(Bos$CPUE +1))
#hm.

#linear model.

b2 = glm(logCPUE~ SalSurf +  Region + logChl, data = Bos)
summary(b2)

b3 = zeroinfl(rCPUE~ SalSurfsc + I(Month^2) + Month + Temperature + logChlsc, data = Bos, dist = "negbin")
summary(b3)


b5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc + Month + I(Month^2) + logChlsc +(1|Region) + (1|Year), 
             zi = ~ log(SalSurf)+ Secchisc +Month + I(Month^2)+ (1|Year) + (1|Region), data = Bos,
             family= "nbinom2")

summary(b5)
bd = dredge(b5, trace = 2, fixed = c("cond(log(SalSurf))", "zi(log(SalSurf))"), 
            extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
bBest = get.models(bd, 1)[[1]]

summary(bBest)

simresb5 = simulateResiduals(b5)
plot(simresb5)
performance(bBest)

#OK, what zoops do we really want to test?
#Pseudodiaptomus
#Eurytemora
#Acartiella
#Daphnia
#BOsmina
#Two rotifers
#Limnoithona
#Acanthocyclops
#make sure there is at least one predatory copepod. 

#limnoithona

Limno = Filters("Limnoithona", data = ZoopMaster2m)
Limno2 = filter(ZoopMaster2m, Taxname2 == "Limnoithona")
hist(Limno$CPUE)
hist(log10(Limno$CPUE +1))
hist(sqrt(Limno$CPUE))


hist(log(Limno2$CPUE +1))
#hm. This one is just odd looking.
#not particularly zero-inflated though


l5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc+ Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ log(SalSurf) + Secchisc+ Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), data = Limno, 
             family= "nbinom2")

summary(l5)

ld = dredge(l5, trace = 2, fixed = c("cond(Month)", "cond(log(SalSurf))", "zi(log(SalSurf))"), extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
lbest = get.models(ld, 1)[[1]]
summary(lbest)


simresl5 = simulateResiduals(lbest)
plot(simresl5)
performance(lbest)

#Some random rotifer

Kerat = Filters("Keratella", data = ZoopMaster2m) %>%
  filter(!is.na(Secchisc))
hist(Kerat$CPUE)
hist(log(Kerat$CPUE +1))
#perfectly bell-shaped except for the zeros

k5 = glmmTMB(rCPUE ~log(SalSurf) + Secchisc+ Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ log(SalSurf) + Secchisc+ Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), data = Kerat, 
             family= "nbinom2", na.action = "na.fail")

summary(k5)
kd = dredge(k5, trace = 2,fixed = c("cond(Month)", "cond(I(Month^2))", "zi(log(SalSurf))"), 
            extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
kbest = get.models(kd, 1)[[1]]
summary(kbest)
simresk5 = simulateResiduals(kbest)
plot(simresk5)



#another rotifer

Sync = Filters("Synchaeta", data = ZoopMaster2m) %>%
  filter(!is.na(Secchisc))
hist(Sync$CPUE)
hist(log(Sync$CPUE +1))
#perfectly bell-shaped except for the zeros

s5 = glmmTMB(rCPUE ~log(SalSurf) +Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ log(SalSurf) +Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), data = Sync, 
             family= "nbinom2", na.action = "na.fail")

summary(s5)
ds = dredge(s5, trace = 2, fixed = c("cond(Month)","cond(I(Month^2))"),
            extra = list("Rs" = function(x) {
              s <- performance(x)
              c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
            }))
sBest = get.models(ds, 1)[[1]]
summary(sBest)
simresd = simulateResiduals(sBest)
plot(simresd)
performance(sBest)

#Hyperacanthomysis logirostris
#this is the version without FMWT
Hyp2 = Filters("Hyperacanthomysis longirostris", data = ZoopMaster2MM2)

hist(Hyp2$CPUE)
hist(log(Hyp$CPUE+1))
# a little flat and zero-inflated, bu tnot bad
h52 = glmmTMB(rCPUE ~log(SalSurf) +Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ Month + I(Month^2)+log(SalSurf) +Secchisc + logChlsc +(1|Region) + (1|Year), data = Hyp2, 
             family= "nbinom2", na.action = "na.fail")

summary(h52)
dh = dredge(h52, trace = 2, fixed = c("cond(Month)"),  extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
hBest2 = get.models(dh, 1)[[1]]
summary(hBest2)
simresd = simulateResiduals(hBest2)
plot(simresd)
check_collinearity(hBest2)
performance(hBest2)
#visreg(hBest)
#OK, works WAY better without FMWT

#Now with FMWT
Hyp = Filters("Hyperacanthomysis longirostris", data = ZoopMaster2MM) %>%
  filter(Year >2009)

hist(Hyp$CPUE)
hist(log(Hyp$CPUE+1))
# a little flat and zero-inflated, bu tnot bad
h5 = glmmTMB(rCPUE ~log(SalSurf) +Secchisc + Month + I(Month^2) + logChlsc + (1|Region) + (1|Year), 
             zi = ~ Month + I(Month^2)+log(SalSurf) +Secchisc + logChlsc +(1|Region) + (1|Year), data = Hyp, 
             family= "nbinom2", na.action = "na.fail")

summary(h5)
dh = dredge(h5, trace = 2, fixed = c("cond(Month)"),  extra = list("Rs" = function(x) {
  s <- performance(x)
  c(Rsqc = s$R2_conditional, Rsqm = s$R2_marginal)
}))
hBest = get.models(dh, 1)[[1]]
summary(hBest)
simresd = simulateResiduals(hBest)
plot(simresd)
check_collinearity(hBest)
performance(hBest)
#visreg(hBest)

save(ad, pd, ed, acd, bd, dh, ds, kd, ld, aBest, pBest, eBest, hBest,  acBest, bBest, sBest, kbest, lbest, file = "SpeciesModels_5SEP2021.RData")
#load("SpeciesModels_5SEP2021.RData")

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
x = p5
agrigate = bind_rows(getcoefs(pBest, "Pseudodiaptomus"),
                     getcoefs(eBest, "Eurytemora"),
                     getcoefs(bBest, "Bosmina"),
                     getcoefs(kbest, "Kerottela"),
                     getcoefs(lbest, "Limnoithona"),
                     getcoefs(hBest, "H. Logirostris"),
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

test3 = mutate(test3, ParamX = factor(Param, levels = c("(Intercept)", "log(SalSurf)", "Secchisc",
                                                       "Month", "I(Month^2)", "logChlsc"),
                                     labels = c("Intercept", "Salinity", "Secchi Depth", 
                                            "Month", "Month^2", "Chlorophyll")),
               Paramtype = factor(Paramtype, levels = c("count", "zi"), labels =
                                    c("Count", "Zero Inflation")))
foo = data.frame(Species = unique(test2$Species), Paramtype = "Zero Inflation", 
                 ParamX = "Chlorophyll", Estimate2 = "Not Included")
test3 = bind_rows(test3, foo)

ggplot(test3, aes(x = ParamX, y = Species, fill = Estimate2)) + geom_tile() +
  facet_wrap(~Paramtype) +
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


#Think about how best to display this. May not be fair to put everyone on the same scale IDK
test3 = mutate(test3, sig = case_when(
  `Pr(>|z|)` < 0.05 ~ "p <0.05",
  `Pr(>|z|)` > 0.05 ~ "NS"
))

#Point plot with standard error
ggplot(filter(test3, ParamX != "Intercept"), aes(x = Species, y = Estimate, color = ParamX, alpha = sig)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Estimate- `Std. Error`, ymax = Estimate + `Std. Error`))+
  facet_grid(ParamX~Paramtype, scales = "free_y")+
  scale_alpha_manual(values = c(0.3, 1))+ 
  theme(axis.text.x = element_text(angle = 90))


#Bar plot
ggplot(filter(test3, ParamX != "Intercept"), aes(x = Species, y = Estimate, fill = ParamX, alpha = sig)) + 
  geom_bar(stat = "identity") + facet_grid(ParamX~Paramtype, scales = "free_y")+
  scale_alpha_manual(values = c(0.3, 1)) 

######################################################3
################################################################################

#stick them all together and plot them

ZoopSub = bind_rows(filter(ZoopMaster2, Taxname2 %in% c("Pseudodiaptomus", "Eurytemora affinis",
                                               "Bosmina longirostris", "Acartiella", "Acanthocyclops")),
                    filter(ZoopMaster2m, Taxname2 %in% c("Synchaeta","Keratella",
                                                        "Limnoithona")),
                    filter(ZoopMaster2MM2, Taxname2 == "Hyperacanthomysis longirostris"))

ggplot(filter(ZoopSub, Region != "North"), aes(x = Month, y = logCPUE, color = Taxname2)) + geom_point()+
  geom_smooth() + facet_wrap(~Region)

ggplot(filter(ZoopSub, Region != "North"), aes(x = Month, y = logCPUE, color = Region)) + geom_point()+
  geom_smooth() + facet_wrap(~Taxname2)

ggplot(filter(ZoopSub, Region != "North"), aes(x = logChl, y = logCPUE, color = Taxname2)) + geom_point()+
  geom_smooth() + facet_wrap(~Month)

ggplot(filter(ZoopSub, Region != "North"), aes(x = logChl, y = logCPUE)) + geom_point()+
  geom_smooth() + facet_grid(Month~Taxname2)

#####################################################################################3
#multivariates, this time with summary data
#convert to a community matrix

#######################################################\#
#Mesozoops
#Zoop = bind_rows(ZoopMasterm, ZoopMaster, ZoopMasterMM)
Zoop = dplyr::select(ZoopMaster, Month, Season, Region, Year, CPUE, Taxname2, SalSurf,
                     Temperature, Secchi, Chlorophyll) %>%
  mutate(Region = factor(Region, levels = c("Confluence", "North","SouthCentral","Suisun Bay","Suisun Marsh")) ) %>%
  filter(Region != "North", !is.na(Secchi))




MasterCom = Zoop %>%
  dplyr::filter(!is.na(Taxname2)) %>%
  pivot_wider(names_from =  Taxname2, values_from = CPUE, 
              values_fn = sum, id_cols = c(Month, Season, Region, Year, SalSurf, 
                                           Temperature, Region, Secchi, Chlorophyll))# %>%
#  filter(!is.na(`Neomysis mercedis`))

#seperate environmental and species matrix
MastenvX = dplyr::select(MasterCom, Month, Season, Region, Year, 
                        SalSurf, Temperature, Secchi, Chlorophyll) 


Mastenv = dplyr::select(MasterCom, Month, Season, Region, Year, Secchi, 
                        SalSurf, Temperature, Chlorophyll) %>%
  mutate(logChl = log(Chlorophyll +1)) %>%
  mutate(across(where(is.numeric), scale))
Mastenv[] <- lapply(Mastenv, function(x) { attributes(x) <- NULL; x })
Mastenv = mutate(Mastenv, 
                 Region = factor(Region, levels = c(1,2,3,4,5), labels = c("Confluence", "North",
                                                    "SouthCentral","Suisun Bay",
                                                    "Suisun Marsh")),
                 Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
                 logChl = log(Chlorophyll +1))

MastCOM = dplyr::select(MasterCom, "Acanthocyclops":last_col())

adonis(MastCOM ~logChl+ SalSurf + Month + Temperature, data = Mastenv)
adonis(MastCOM ~Chlorophyll+ SalSurf + Month + Temperature, data = Mastenv)

nmds = metaMDS(MastCOM, trymax = 200)
source("plotNMDS.R")
PlotNMDS2(nmds, Mastenv, textp = F, lines = "SalSurf", group = "Region")
PlotNMDS(nmds, Mastenv, textp = T, group = "Region")

#permutation design
h1 <- how(within = Within(type = "free"), plots = Plots(strata = Mastenv$Region),
          blocks = Mastenv$Year)

ef = envfit(nmds ~ SalSurf + logChl + Month  + Secchi, Mastenv, na.rm = T,
            permutations = h1)
ef
plot(ef)
#####################################################################
#Microzoops
ZoopMicro = dplyr::select(ZoopMasterm, Month, Season, Region, Year, CPUE, Taxname2, SalSurf,
                     Temperature, Secchi, Chlorophyll) %>%
  mutate(Region = factor(Region, levels = c("Confluence", "North","SouthCentral","Suisun Bay","Suisun Marsh")) ) %>%
  filter(!is.na(Secchi), Year >2000)

MasterComMicro = ZoopMicro %>%
  dplyr::filter(!is.na(Taxname2)) %>%
  pivot_wider(names_from =  Taxname2, values_from = CPUE, 
              values_fn = sum, id_cols = c(Month, Season, Region, Year, SalSurf, 
                                           Temperature, Region, Secchi, Chlorophyll)) %>%
  mutate(Oithona = Oithona + `Oithona davisae`, `Oithona davisae`= NULL)

#seperate environmental and species matrix
MastenvXMicro = dplyr::select(MasterComMicro, Month, Season, Region, Year, 
                         SalSurf, Temperature, Secchi, Chlorophyll) 


MastenvMicro = dplyr::select(MasterComMicro, Month, Season, Region, Year, Secchi, 
                        SalSurf, Temperature, Chlorophyll) %>%
  mutate(logChl = log(Chlorophyll +1)) %>%
  mutate(across(where(is.numeric), scale))
MastenvMicro[] <- lapply(MastenvMicro, function(x) { attributes(x) <- NULL; x })
MastenvMicro = mutate(MastenvMicro, 
                 Region = factor(Region, levels = c(1,2,3,4,5), labels = c("Confluence", "North",
                                                                           "SouthCentral","Suisun Bay",
                                                                           "Suisun Marsh")),
                 Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
                 logChl = log(Chlorophyll +1))

MastCOMMicro = dplyr::select(MasterComMicro, "Asplanchna":last_col())

adonis(MastCOMMicro ~logChl+ SalSurf + Month + Temperature + Secchi, data = MastenvMicro)
adonis(MastCOMMicro ~Chlorophyll+ SalSurf + Month + Temperature, data = MastenvMicro)

nmdsM = metaMDS(MastCOMMicro, k=2,trymax = 400)
nmdsM3 = metaMDS(MastCOMMicro, k=3,trymax = 400)

PlotNMDS(nmdsM3, MastenvMicro, textp = T, group = "Region")
nmdsM
#so, it works with K=3 and a smaller data set. 
#permutation design
h1Mic <- how(within = Within(type = "free"), plots = Plots(strata = MastenvMicro$Region),
             blocks = MastenvMicro$Year)
ef3 = envfit(nmdsM3 ~ Secchi+ SalSurf + logChl + Month, MastenvMicro, na.rm = T,
             permutations = h1Mic)
ef3
plot(ef3)
########################################################################

#Mysids
ZoopMac = dplyr::select(ZoopMasterMM2, Month, Season, Region, Year, CPUE, Taxname2, SalSurf,
                          Temperature, Secchi, Chlorophyll) %>%
  mutate(Region = factor(Region, levels = c("Confluence", "North","SouthCentral","Suisun Bay","Suisun Marsh"))) %>%
  filter(!is.na(Secchi), Year > 2005) 

MasterComMac = ZoopMac %>%
  dplyr::filter(!is.na(Taxname2)) %>%
  pivot_wider(names_from =  Taxname2, values_from = CPUE, 
              values_fn = sum, id_cols = c(Month, Season, Region, Year, SalSurf, 
                                           Temperature, Region, Secchi, Chlorophyll))
MasterComMac$tot = rowSums(MasterComMac[,9:16])
MasterComMac = filter(MasterComMac, tot != 0)

#seperate environmental and species matrix
MastenvXMac = dplyr::select(MasterComMac, Month, Season, Region, Year, 
                              SalSurf, Temperature, Secchi, Chlorophyll) 


MastenvMac = dplyr::select(MasterComMac, Month, Season, Region, Year, Secchi, 
                             SalSurf, Temperature, Chlorophyll) %>%
  mutate(logChl = log(Chlorophyll +1)) %>%
  mutate(across(where(is.numeric), scale))
MastenvMac[] <- lapply(MastenvMac, function(x) { attributes(x) <- NULL; x })
MastenvMac = mutate(MastenvMac, 
                      Region = factor(Region, levels = c(1,2,3,4,5), labels = c("Confluence", "North",
                                                                                "SouthCentral","Suisun Bay",
                                                                                "Suisun Marsh")),
                      Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
                      logChl = log(Chlorophyll +1))

MastCOMMac = dplyr::select(MasterComMac, `Alienacanthomysis macropsis`:last_col(), -tot)

adonis(MastCOMMac ~logChl+ log(SalSurf+1) + Month + Temperature + Secchi, data = MastenvMicro)
adonis(MastCOMMac ~Chlorophyll+ log(SalSurf) + Month + Temperature, data = MastenvMicro)

nmdsMac = metaMDS(MastCOMMac, k=2,trymax = 1000)
nmdsMac3 = metaMDS(MastCOMMac, k=3,trymax = 400)

PlotNMDS(nmdsMac3, MastenvMac, textp = T, group = "Region")


#permutation design
h1Mac <- how(within = Within(type = "free"), plots = Plots(strata = MastenvMac$Region),
          blocks = MastenvMac$Year)
ef2 = envfit(nmdsMac3 ~ SalSurf + logChl + Month +  Secchi, MastenvMac, na.rm = T,
            permutations = h1Mac)
ef2
plot(ef2)
################################################################

############################################################
#quick plot of all the cladocerans
clads = ZoopMaster2 %>%
  filter(Taxname2 %in% c("Bosmina longirostris", "Daphnia", "Diaphanosoma", "Cladocera"))
ggplot(clads, aes(x = Taxname2, y = CPUE)) + geom_col()

######################################################################
#what is the relationship between chlorophyll and month?
hist(ZoopMaster$Chlorophyll)
hist(ZoopMaster$logChl)

p5 = lmer(logChl ~log(SalSurf)+ Secchi + Month + I(Month^2) + (1|Region) + (1|Year),  
             data= ZoopMaster)
summary(p5)
visreg(p5)

library(mgcv)
WQ = filter(WQ_dataNARM, Year > 1995 , Month %in% 3:10)
g1 = gam(logChl ~ s(Conductivity, k = 3)+ s(Secchi, k=3) + s(Month, k = 5), data = WQ)
plot(g1)


