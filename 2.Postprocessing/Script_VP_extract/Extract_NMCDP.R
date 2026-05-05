#
#
#           TEST VP full model & all taxa
#

## Function prepared
setwd("~/Code/3.Postprocessing/Script_VP_extract")
source('computeVarPartSummaries_111125.R')

setwd("~/Code/x.Extra_functions")
source('cHmsc.R')

## libraries
library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)

## extractions

# #
# #    Occurrence
# #
##.......................................................................................................
##  Rodents
##.......................................................................................................
setwd("~/Data/3.Modelrun/occ")

load('Rodtest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('Rodtest2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('Rodtest2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('Rodtest2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models_rod = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

## get performances for weighting
preds = computePredictedValues(models_rod, expected = FALSE)
MF = evaluateModelFit(hM=models_rod, predY=preds)

TabMF_rod = cbind(as.data.frame(MF), species = models_rod$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_rod, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_rod = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_rod, TabMF_rod, file = 'VPCP_rod.RDATA')



##.......................................................................................................
##  Butterflies
##.......................................................................................................
setwd("~/Data/3.Modelrun/occ")

load('bftest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bftest_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bftest_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bftest_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models_bf = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_bf, expected = FALSE)
MF = evaluateModelFit(hM=models_bf, predY=preds)

TabMF_bf = cbind(as.data.frame(MF), species = models_bf$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_bf, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = TRUE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_bf = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_bf, TabMF_bf, file = 'VPCP_bf.RDATA')



##.......................................................................................................
##  Moths
##.......................................................................................................
setwd("~/Data/3.Modelrun/occ")

load('MothTest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothTest2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothTest2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothTest2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models_moth = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_moth, expected = FALSE)
MF = evaluateModelFit(hM=models_moth, predY=preds)

TabMF_Moth = cbind(as.data.frame(MF), species = models_moth$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_moth, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_moth = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_moth, TabMF_Moth, file = 'VPCP_moth.RDATA')


##.......................................................................................................
##  Birds
##.......................................................................................................
setwd("~/Data/3.Modelrun/occ")

load('Birdtest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN1 = models
load('Birdtest_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN2 = models
load('Birdtest_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN3 = models
load('Birdtest_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN4 = models

models_bd = c.Hmsc(modelsBirdN1[[1]], modelsBirdN2[[1]], modelsBirdN3[[1]], modelsBirdN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_bd, expected = FALSE)
MF = evaluateModelFit(hM=models_bd, predY=preds)

TabMF_bd = cbind(as.data.frame(MF), species = models_bd$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_bd, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_bd = VP.test3
summary(VP.test3_bd$Vnorm)


# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_bd, TabMF_bd, file = 'VPCP_bd.RDATA')



##.......................................................................................................
##  Winter game
##.......................................................................................................
setwd("~/Data/3.Modelrun/occ")

load('WintGtest_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGtest_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGtest_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGtest_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models

models_wg = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_wg, expected = FALSE)
MF = evaluateModelFit(hM=models_wg, predY=preds)

TabMF_wg = cbind(as.data.frame(MF), species = models_wg$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_wg, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_wg = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_wg, TabMF_wg, file = 'VPCP_wg.RDATA')



# #
# #    Abundance
# #
##.......................................................................................................
##  Rodents
##.......................................................................................................
setwd("~/Data/3.Modelrun/ab")

load('RodtestAB_ch1_models_thin_1000_samples_500_chains_1.Rdata')
modelsN1 = models
load('RodtestAB_ch2_models_thin_1000_samples_500_chains_1.Rdata')
modelsN2 = models
load('RodtestAB_ch3_models_thin_1000_samples_500_chains_1.Rdata')
modelsN3 = models
load('RodtestAB_ch4_models_thin_1000_samples_500_chains_1.Rdata')
modelsN4 = models


models_rodAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

## get performances for weighting
preds = computePredictedValues(models_rodAB, expected = FALSE)
MF = evaluateModelFit(hM=models_rodAB, predY=preds)

TabMF_rod = cbind(as.data.frame(MF), species = models_rodAB$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_rodAB, exclude = 2, group = group1,
                                   groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_rod = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_rod, TabMF_rod, file = 'VPCP_ab_rod.RDATA')



##.......................................................................................................
##  Butterflies
##.......................................................................................................
setwd("~/Data/3.Modelrun/ab")

load('bfAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bfAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bfAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bfAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models_bfAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_bfAB, expected = FALSE)
MF = evaluateModelFit(hM=models_bfAB, predY=preds)

TabMF_bf = cbind(as.data.frame(MF), species = models_bfAB$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_bfAB, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_bf = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_bf, TabMF_bf, file = 'VPCP_ab_bf.RDATA')



##.......................................................................................................
##  Moths
##.......................................................................................................
setwd("~/Data/3.Modelrun/ab")


load('MothAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models_mothAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_mothAB, expected = FALSE)
MF = evaluateModelFit(hM=models_mothAB, predY=preds)

TabMF_Moth = cbind(as.data.frame(MF), species = models_mothAB$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_mothAB, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_moth = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_moth, TabMF_Moth, file = 'VPCP_ab_moth.RDATA')


##.......................................................................................................
##  Birds
##.......................................................................................................
setwd("~/Data/3.Modelrun/ab")

load('BirdAB2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN1 = models
load('BirdAB2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN2 = models
load('BirdAB2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN3 = models
load('BirdAB2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN4 = models

models_bdAB = c.Hmsc(modelsBirdN1[[1]], modelsBirdN2[[1]], modelsBirdN3[[1]], modelsBirdN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_bdAB, expected = FALSE)
MF = evaluateModelFit(hM=models_bdAB, predY=preds)

TabMF_bd = cbind(as.data.frame(MF), species = models_bdAB$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_bdAB, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_bd = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_bd, TabMF_bd, file = 'VPCP_ab_bd.RDATA')



##.......................................................................................................
##  Winter game
##.......................................................................................................
setwd("~/Data/3.Modelrun/ab")

load('WintGAB2_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGAB2_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGAB2_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGAB2_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models

models_wgAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models_wgAB, expected = FALSE)
MF = evaluateModelFit(hM=models_wgAB, predY=preds)

TabMF_wg = cbind(as.data.frame(MF), species = models_wgAB$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

##.......................................................................................................
## marginal

VP.test3 = computeVarPartSummaries(hM=models_wgAB, exclude = 2, group = group1,
                                           groupnames = groupnames1, marginal=TRUE, conditional = FALSE, Global = FALSE)
dim(VP.test3$Vnorm)
dim(VP.test3$Vdiag)
dim(VP.test3$Vmarg)
dim(VP.test3$Vpart)
dim(VP.test3$Cnorm)

VP.test3_wg = VP.test3

# Save
setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(VP.test3_wg, TabMF_wg, file = 'VPCP_ab_wg.RDATA')


