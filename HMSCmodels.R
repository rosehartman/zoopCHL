#HMSC models



##################################################################33
#put it in a format that HMSC will like  


#I might want to try presence/absence and abundance sperately (hurdle model) but
#i'll start with just a basic model

#set up HMSC model
EMPm2 = Hmsc(Y = as.matrix(MastCOM), XData = Mastenv, 
             XFormula = ~Month + Chlorophyll + SalSurf + Year, distr = "lognormal poisson")

#go!
EMPmm2 = sampleMcmc(EMPm2, thin = 5, samples = 3000, transient = 100,
                    nChains = 2)

#now with log-transformed chlorophyll for comparison
EMPm2.1 = Hmsc(Y = as.matrix(MastCOM), XData = Mastenv, 
               XFormula = ~I(Month^2) + logChl + SalSurf + Year, distr = "lognormal poisson")

EMPmm2.1 = sampleMcmc(EMPm2.1, thin = 5, samples = 3000, transient = 100,
                      nChains = 2)


EMPm2.2 = Hmsc(Y = as.matrix(MastCOM), XData = Mastenv, 
               XFormula = ~I(Month^2) + logChl + SalSurf + Region, distr = "lognormal poisson")

EMPmm2.2 = sampleMcmc(EMPm2.2, thin = 50, samples = 5000, transient = 1000,
                      nChains = 2)


# #Check MCMC convergence diagnostics
EMPmpost = convertToCodaObject(EMPmm2)
diags = data.frame(effectiveSize(EMPmpost$Beta), 
                   gelman.diag(EMPmpost$Beta, multivariate=TRUE, transform = T)$psrf)
# 
hist(diags$effectiveSize.EMPmpost.Beta.)

hist(diags$Point.est.)
#yeah, we'll need more samples, but it's a start

postBeta = getPostEstimate(EMPmm2.2, parName = "Beta")
plotBeta(EMPmm2.2, post = postBeta, param = "Sign", 
         supportLevel = 0.95) 
plotBeta(EMPmm2.1, post = postBeta, param = "Mean")
plotBeta(EMPmm2.1, post = postBeta, param = "Support")
#Huh. I don't know if this deals with zero inflation. Also might want to scale everything.


Specieslevels = c( "Acartia", "Acartiella", "Acartiella", "Diaptomidae",
                   "Eurytemora affinis", "Pseudodiaptomus","Sinocalanus", 
                   "Tortanus", "Calanoida",
                   "Acanthocyclops",   "Cyclopoida",                            
                   "Harpacticoida","Bosmina longirostris",  "Daphnia" , "Diaphanosoma",
                   "Cladocera","Cirripedia", "Decapoda")

Specieslabels = c( "Acartia", "Acartiella", "Acartiella", "Diaptomidae",
                   "Eurytemora affinis", "Pseudodiaptomus","Sinocalanus", 
                   "Tortanus", "Other Calanoida",
                   "Acanthocyclops",   "Other Cyclopoida",                             
                   "Harpacticoida","Bosmina",  "Daphnia" , "Diaphanosoma",
                   "Other Cladocera","Cirripedia", "Decapoda")
CovLevels = EMPmm2$covNames
CovLabels = c("Intercept", "Month", "Chlorophyll", "Salinity", "Year")

source("plotBeta.R")
plotBeta2(EMPmm2, post = postBeta, supportLevel = 0.9, covlabels = CovLabels)
plotBeta3(EMPmm2, post = postBeta, supportLevel = 0.9,  covlabels = CovLabels)


#################################################################
#try it as a hurdle model

MastCOMPA = MastCOM
MastCOMPA[MastCOMPA>0] = 1
MastCOMab = MastCOM
MastCOMab[MastCOM==0] = NA

#first run presence/absence
Hm1 = Hmsc(Y = as.matrix(MastCOMPA), XData = Mastenv, 
           XFormula = ~I(Month^2) + logChl + SalSurf + Year, distr = "probit")

Hmm1 = sampleMcmc(Hm1, thin = 5, samples = 1000, transient = 100,
                  nChains = 2)

#Now abundance conditional on presence
Hm2 = Hmsc(Y = as.matrix(MastCOMab), XData = Mastenv, 
           XFormula = ~I(Month^2) + logChl + SalSurf + Year, distr = "lognormal poisson")

Hmm2 = sampleMcmc(Hm2, thin = 50, samples = 10000, transient = 1000,
                  nChains = 2)


#swap out year for regiona

#first run presence/absence
Hm1.1 = Hmsc(Y = as.matrix(MastCOMPA), XData = Mastenv, 
             XFormula = ~I(Month^2) + logChl + SalSurf + Region, distr = "probit")

Hmm1.1 = sampleMcmc(Hm1.1, thin = 50, samples = 10000, transient = 1000,
                    nChains = 2)

#Now abundance conditional on presence
Hm2.1 = Hmsc(Y = as.matrix(MastCOMab), XData = Mastenv, 
             XFormula = ~I(Month^2) + logChl + SalSurf + Region, distr = "lognormal poisson")

Hmm2.1 = sampleMcmc(Hm2.1, thin = 50, samples = 10000, transient = 1000,
                    nChains = 2)



postBetaH.1 = getPostEstimate(Hmm1.1, parName = "Beta")
plotBeta(Hmm1.1, post = postBetaH.1, param = "Sign", 
         supportLevel = 0.95) 
plotBeta(Hmm1, post = postBetaH, param = "Sign", 
         supportLevel = 0.95) 
#Chlorophyll story is similar by Year and Region


postBetaH2.1 = getPostEstimate(Hmm2.1, parName = "Beta")
plotBeta(Hmm2, post = postBetaH2, param = "Sign", 
         supportLevel = 0.95) 
plotBeta(Hmm2.1, post = postBetaH2.1, param = "Sign", 
         supportLevel = 0.95) 

#I feel like "region" is confusing and doesn't get at the basic question.
#My basic quesiton is about abiotic variables.

save(Hmm1, Hmm2, EMPmm2, EMPmm2.1, EMPmm2.2, Hmm1.1, Hmm2.1, file = "ZoopSal9AUG2021.RData")

# #Check MCMC convergence diagnostics
HM = convertToCodaObject(Hmm2.1)
diags = data.frame(effectiveSize(HM$Beta), 
                   gelman.diag(HM$Beta, multivariate=TRUE, transform = T)$psrf)
# 
hist(diags$effectiveSize.HM.Beta.)

hist(diags$Point.est.)
#That's pretty good.

################################################################3
#I think I want to concentrate on aboitic varialbes.
#But maybe use "region" as a blocking term? Year too?

Mastenv = droplevels(Mastenv)
StudyDesign = data.frame(Sample = as.factor(1:nrow(Mastenv)), Region = as.factor(Mastenv$Region))
rL = HmscRandomLevel(units=as.factor(StudyDesign$Sample)) 
rL2 = HmscRandomLevel(units=as.factor(StudyDesign$Region)) 


#first run presence/absence
Hm1.2 = Hmsc(Y = as.matrix(MastCOMPA), XData = Mastenv, 
             XFormula = ~Month + I(Month^2) + logChl + SalSurf +Secchi, distr = "probit",
             ranLevels=list("Sample"=rL, "Region" = rL2), studyDesign = StudyDesign)

Hmm1.2 = sampleMcmc(Hm1.2, thin = 10, samples = 3000, transient = 100,
                    nChains = 2)

#Now abundance conditional on presence
Hm2.2 = Hmsc(Y = as.matrix(MastCOMab), XData = Mastenv, 
             XFormula = ~Month + I(Month^2) + logChl + SalSurf + Secchi, distr = "lognormal poisson",
             ranLevels=list("Sample"=rL, "Region" = rL2), studyDesign = StudyDesign)

Hmm2.2 = sampleMcmc(Hm2.2, thin = 10, samples = 3000, transient = 100,
                    nChains = 2)


# #Check MCMC convergence diagnostics
HM = convertToCodaObject(Hmm1.2)
diags = data.frame(effectiveSize(HM$Beta), 
                   gelman.diag(HM$Beta, multivariate=TRUE, transform = T)$psrf)
# 
hist(diags$effectiveSize.HM.Beta.)
hist(diags$Point.est.)




postBetaH1.2 = getPostEstimate(Hmm1.2, parName = "Beta")
plotBeta(Hmm1.2, post = postBetaH1.2, param = "Sign", 
         supportLevel = 0.95) 


postBetaH.2 = getPostEstimate(Hmm2.2, parName = "Beta")
plotBeta(Hmm2.2, post = postBetaH.2, param = "Sign", 
         supportLevel = 0.95) 
#Chlorophyll story is similar by Year and Region


plotBeta2(Hmm1.2, post = postBetaH1.2, supportLevel = 0.9, 
          covlabels = c("Intercept","Month", "Month^2", "Chlorophyll",
                        "Salinity", "Secchi"))
plotBeta2(Hmm2.2, post = postBetaH.2, supportLevel = 0.9, 
          covlabels = c("Intercept","Month", "Month^2", "Chlorophyll",
                        "Salinity", "Secchi"))

#########################################################################
#Amount of variance explained by the JSDM

#Count model
VP = computeVariancePartitioning(Hmm2.2)
valsx = as.data.frame(t(VP$vals))
vals = mutate(valsx, Species = rownames(valsx)) %>%
  pivot_longer(cols = c(-Species), names_to = "Factor", values_to = "Variance") %>%
  mutate(Factor = factor(Factor, levels = names(valsx),
                         labels = c("Month", "Month^2", "Chlorophyll", "Salinity", "Secchi",
                                    "Random(Sample)", "Random(Region)")),
         Species = factor(Species, levels = Specieslevels, labels = Specieslabels))

#To assess the model’s explanatory power, we apply the evaluateModelFit 
#function to the posterior predictive
#distribution simulated by the function computePredictedValues.
preds = computePredictedValues(Hmm2.2, nParallel = 2)
fit = evaluateModelFit(hM = Hmm2.2, predY = preds)
fit
#
fit2 = as.data.frame(fit) %>%
  mutate(Species = Hmm2.2$spNames)


#put it in a format so we can plot it
fit2.1 = mutate(fit2, Species = factor(Species, levels = Specieslevels, 
                                       labels = Specieslabels),)
vals2 = left_join(vals, fit2.1) %>%
  mutate(RawVarience = SR2*Variance)

mycols = rainbow(nrow(VP$vals))
plotVariancePartitioning(hM=Hmm2.2, VP=VP,cols = mycols, 
                         args.leg=list(bg="white",cex=0.7),
                         main = "Proportion of explained variance, ",cex.main=0.8,
                         las = 2)
#call to ggplot
ggplot(vals, aes(x=Species, y = Variance, fill = Factor)) + 
  geom_col()+
  scale_fill_brewer(type = "qual")+
  ylab("Proportion Explained Varience ")+
  theme_bw()+ 
  theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))

mytheme =  theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))


#Raw variance plot - count model
ggplot(vals2, aes(x=Species, y = RawVarience, fill = Factor)) + geom_col()+
  scale_fill_brewer(type = "qual")+
  ylab("Proportion Raw Varience")+
  mytheme


# theme_bw()+ theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))



#Zeros model
VP0 = computeVariancePartitioning(Hmm1.2)
valsx0 = as.data.frame(t(VP0$vals))
vals0 = mutate(valsx0, Species = rownames(valsx0)) %>%
  pivot_longer(cols = c(-Species), names_to = "Factor", values_to = "Variance") %>%
  mutate(Factor = factor(Factor, levels = names(valsx0),
                         labels = c("Month", "Month^2", "Chlorophyll", "Salinity", "Secchi",
                                    "Random(Sample)", "Random(Region)")),
         Species = factor(Species, levels = Specieslevels, labels = Specieslabels))

#To assess the model’s explanatory power, we apply the evaluateModelFit 
#function to the posterior predictive
#distribution simulated by the function computePredictedValues.
preds0 = computePredictedValues(Hmm1.2, nParallel = 2)
fit0 = evaluateModelFit(hM = Hmm1.2, predY = preds0)
fit0
#
fit20 = as.data.frame(fit0) %>%
  mutate(Species = Hmm1.2$spNames)


#put it in a format so we can plot it
fit2.10 = mutate(fit20, Species = factor(Species, levels = Specieslevels, 
                                         labels = Specieslabels),)
vals20 = left_join(vals0, fit2.10) %>%
  mutate(RawVarience = TjurR2*Variance)

mycols = rainbow(nrow(VP$vals))
plotVariancePartitioning(hM=Hmm1.2, VP=VP0,cols = mycols, 
                         args.leg=list(bg="white",cex=0.7),
                         main = "Proportion of explained variance, ",cex.main=0.8,
                         las = 2)

#proportional variance plot - zeros
ggplot(vals0, aes(x=Species, y = Variance, fill = Factor)) + geom_col()+
  scale_fill_brewer(type = "qual")+
  ylab("Proportion Explained Varience ")+
  theme_bw()+ theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))

#Raw variance plot - Zeros
ggplot(vals20, aes(x=Species, y = RawVarience, fill = Factor)) + geom_col()+
  scale_fill_brewer(type = "qual")+
  ylab("Proportion Raw Varience")+
  theme_bw()+ theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))+
  coord_cartesian(ylim = c(0,1))


#Now plot just the fixed effects
#Raw variance plot - Zeros

filter(vals20, Factor != "Random(Sample)" & Factor != "Random(Region)") %>%
  ggplot( aes(x=Species, y = RawVarience, fill = Factor)) + geom_col()+
  scale_fill_brewer(type = "qual")+
  ylab("Proportion Raw Varience")+
  theme_bw()+ theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))+
  coord_cartesian(ylim = c(0,1))


filter(vals2, Factor != "Random(Sample)" & Factor != "Random(Region)") %>%
  ggplot( aes(x=Species, y = RawVarience, fill = Factor)) + geom_col()+
  scale_fill_brewer(type = "qual")+
  ylab("Proportion Raw Varience")+
  theme_bw()+ theme(axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90))+
  coord_cartesian(ylim = c(0,1))


df = mutate(df, Species = case_when(
  Species == "Pseudoiaptomus copepodid" ~ "Calanoid copepodid",
  Species == "Sinocalanus copepodid" ~ "Calanoid copepodid",
  TRUE ~ Species
))