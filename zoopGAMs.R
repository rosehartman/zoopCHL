#OK, Apparently I really need to use GAMs. Let's do this. 

#

library(lubridate)
library(tidyverse)
library(brms)
library(tidybayes)
library(ggdist)
library(rstan)
library(cmdstanr)


source("HelperFuncitons.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("OrganizedZoops.RData")

P1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Pseudos, family = "hurdle_lognormal", 
         control=list(adapt_delta=0.995),
         backend = "cmdstan")
pairs(P1)
plot(P1)
summary(P1)
conditional_smooths(P1)
plot(conditional_effects(P1), points = TRUE)
pp_check(P1)
pp_check(P1, type = "ecdf_overlay")

#try to add the predictions
#Pseudos2 = add_epred_draws(P1, Pseudos)
PseudoPred = zoop_predict(P1, Pseudos)

#graph differences in chlorophyll by salinity for 
pxchl = filter(PseudoPred, Month == 6, Secchisc == unique(Secchisc)[10], 
            logChlsc %in% unique(logChlsc)[c(1,7,14,20)])
pchl<-ggplot(pxchl, aes(x=SalSurfsc, y=Pred, ymin=lowerCI, ymax=upperCI, 
                          fill= as.factor(logChlsc), group=as.factor(logChlsc)))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=as.factor(logChlsc)))+
  scale_y_continuous(expand=c(0,0), limits = c(0, NA))+
  ylab("CPUE")+
  #facet_wrap(~Secchisc, scales = "free_y")+
  theme_bw()
pchl


#Hmm. I might just have to include a plot of the month when the critter is most 
#abundant to avoid this shit. 

#Also, do I need to scale predictor variables with this new structure? 
#can I back-transform scaled predictors?

pxsec = filter(PseudoPred, Month == 6, logChlsc == unique(logChlsc)[10], 
               Secchisc%in% unique(Secchisc)[c(1,7,13,20)])

psec<-ggplot(pxsec, 
          aes(x=SalSurfsc, y=Pred, ymin=lowerCI, ymax=upperCI, 
              fill= as.factor(Secchisc), group=as.factor(Secchisc)))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=as.factor(Secchisc)))+
  scale_y_continuous(expand=c(0,0), limits = c(0, NA))+
  ylab("CPUE")+
  #facet_wrap(~Secchisc, scales = "free_y")+
  theme_bw()
psec



px2 = filter(PseudoPred, Month == 6, Secchisc == unique(Secchisc)[10], 
             SalSurfsc %in% unique(SalSurfsc)[c(1,7,14,20)])

p2<-ggplot(px2, 
             aes(x=logChlsc, y=Pred, ymin=lowerCI, ymax=upperCI, 
                 fill= as.factor(SalSurfsc), group=as.factor(SalSurfsc)))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=as.factor(SalSurfsc)))+
  scale_y_continuous(expand=c(0,0), limits = c(0, NA))+
  ylab("CPUE")+
  #facet_wrap(~Secchisc, scales = "free_y")+
  theme_bw()
p2

#####################################################################
#OK, Eurytemora

E1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Eury, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(E1)
plot(E1)
conditional_smooths(E1)

EuryPred = zoop_predict(E1, Eury)

#graph differences in chlorophyll by salinity for 
Exchl = filter(EuryPred, Month == 6, Secchisc == unique(Secchisc)[10], 
               logChlsc %in% unique(logChlsc)[c(1,7,14,20)])
echl<-ggplot(Exchl, aes(x=SalSurfsc, y=Pred, ymin=lowerCI, ymax=upperCI, 
                        fill= as.factor(logChlsc), group=as.factor(logChlsc)))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=as.factor(logChlsc)))+
  scale_y_continuous(expand=c(0,0), limits = c(0, NA))+
  ylab("CPUE")+
  #facet_wrap(~Secchisc, scales = "free_y")+
  theme_bw()
echl




#####################################################################
#OK, Acanthocyclops

A1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Acan, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(A1)
plot(A1)
conditional_smooths(A1)

AcanPred = zoop_predict(A1, Acan)

#graph differences in chlorophyll by salinity for 
Acanx = filter(AcanPred, Month == 6, Secchisc == unique(Secchisc)[10], 
               logChlsc %in% unique(logChlsc)[c(1,7,14,20)])
Achl<-ggplot(Acanx, aes(x=SalSurfsc, y=Pred, ymin=lowerCI, ymax=upperCI, 
                        fill= as.factor(logChlsc), group=as.factor(logChlsc)))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=as.factor(logChlsc)))+
  scale_y_continuous(expand=c(0,0), limits = c(0, NA))+
  ylab("CPUE")+
  #facet_wrap(~Secchisc, scales = "free_y")+
  theme_bw()
echl



#####################################################################
#Daphnia

D1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Daphnia, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(D1)
plot(D1)
conditional_smooths(D1)

DaphPred = zoop_predict(D1, Daphnia)



#####################################################################
#Bosmina

B1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Bosmina, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(B1)
plot(B1)
conditional_smooths(B1)

BosPred = zoop_predict(B1, Bosmina)

#####################################################################
#Syncheate
Syn$Yearf = as.factor(Syn$Year)
S1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Syn, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(S1)
plot(S1)
conditional_smooths(S1)

SynPred = zoop_predict(S1, Syn)

#####################################################################
#Kerotella
Kero$Yearf = as.factor(Kero$Year)
K1 = brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Kero, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(K1)
plot(K1)
conditional_smooths(K1)

KeroPred = zoop_predict(K1, Kero)


#####################################################################
#limnos
Limno$Yearf = as.factor(Limno$Year)
L1= brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
           s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
         data= Limno, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(L1)
plot(L1)
conditional_smooths(L1)

LimnoPred = zoop_predict(L1, Limno)

########################################
#hyperacanthomysis
Hyp$Yearf = as.factor(Hyp$Year)
H1= brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
          s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
        data= Hyp, family = "hurdle_lognormal", control=list(adapt_delta=0.995))
summary(H1)
plot(H1)
conditional_smooths(H1)

########################################
#acartiella
Acart$Yearf = as.factor(Acart$Year)
Ac1= brm(CPUE ~ s(SalSurfsc, k = 3)+ s(Secchisc, k = 3) + s(Month, bs = "cc", k =6) + 
          s(logChlsc, k = 3) + (1|Region) + (1|Yearf),  
        data= Acart, family = "hurdle_lognormal", 
        control=list(adapt_delta=0.995, max_treedepth = 15))
summary(Ac1)
plot(Ac1)
conditional_smooths(Ac1)










