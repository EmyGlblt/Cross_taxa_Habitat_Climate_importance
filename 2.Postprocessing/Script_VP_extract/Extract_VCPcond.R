#
#
#            VPcond full model & all taxa extract
#

## Function prepared
source('computeVarPartSummaries_111125.R')
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
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('Rodtest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('Rodtest2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('Rodtest2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('Rodtest2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

## get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_rod = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1


# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)
#cond
#MB  NB  SB 
#382 101 561 

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_rod = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                             groupnames = groupnames1, all=FALSE, marginal=FALSE, 
                                             conditional = TRUE,
                                             cond = cond, wb = NULL) # cannot be both specified at the same time!

dim(VP.test1_rod$Vnorm)  # 3    4   sp 8  nsamp 1000 cond 3  

## ............................................................................................
### grouped VP - WITHIN - BTW

VP.test2_rod = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                             groupnames = groupnames1, all=FALSE, marginal=FALSE, 
                                             conditional = TRUE,
                                             cond = cond, wb = TRUE)

dim(VP.test2_rod$Vnorm)  # 3    4   sp 10  wt nsamp 1000

# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_rod, VP.test2_rod, TabMF_rod, file = 'VPcondwb_Rod.RDATA')



##.......................................................................................................
##  Butterflies
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('bftest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bftest_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bftest_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bftest_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_bf = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_bf = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                 groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                 cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW

VP.test2_bf = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                 groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                 cond = cond, wb = TRUE)


# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_bf, VP.test2_bf, TabMF_bf, file = 'VPcondwb_bf.RDATA')


##.......................................................................................................
##  Moths
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('MothTest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothTest2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothTest2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothTest2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_Moth = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_moth = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW

VP.test2_moth = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                cond = cond, wb = TRUE)


# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_moth, VP.test2_moth, TabMF_Moth, file = 'VPcondwb_moth.RDATA')



##.......................................................................................................
##  Birds
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('Birdtest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN1 = models
load('Birdtest_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN2 = models
load('Birdtest_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN3 = models
load('Birdtest_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN4 = models

models = c.Hmsc(modelsBirdN1[[1]], modelsBirdN2[[1]], modelsBirdN3[[1]], modelsBirdN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_bd = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_bd = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW
VP.test2_bd = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                cond = cond, wb = TRUE)


# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_bd, VP.test2_bd, TabMF_bd, file = 'VPcondwb_bd.RDATA')



##.......................................................................................................
##  Game triangle
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('WintGtest_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGtest_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGtest_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGtest_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_wg = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_wg = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW
VP.test2_wg = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                                groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                                cond = cond, wb = TRUE)

# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_wg, VP.test2_wg, TabMF_wg, file = 'VPcondwb_wg.RDATA')


## all occurrence
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_moth, TabMF_Moth, VP.test1_bf, TabMF_bf,
     VP.test1_rod, TabMF_rod, VP.test1_wg, TabMF_wg,
     VP.test1_bd, TabMF_bd, VP.test2_moth, VP.test2_bf,
     VP.test2_rod, VP.test2_wg,
     VP.test2_bd,file = 'VPCPcond_all_occ.RDATA')




# #
# #    Abundance
# #
##.......................................................................................................
##  Rodents
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('RodtestAB_ch1_models_thin_1000_samples_500_chains_1.Rdata')
modelsN1 = models
load('RodtestAB_ch2_models_thin_1000_samples_500_chains_1.Rdata')
modelsN2 = models
load('RodtestAB_ch3_models_thin_1000_samples_500_chains_1.Rdata')
modelsN3 = models
load('RodtestAB_ch4_models_thin_1000_samples_500_chains_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

## get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_rod = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1


# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)
#cond
#MB  NB  SB 
#382 101 561 

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_rod = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                       groupnames = groupnames1, all=FALSE, marginal=FALSE, 
                                       conditional = TRUE,
                                       cond = cond, wb = NULL)

## ............................................................................................
### grouped VP - WITHIN - BTW

VP.test2_rod = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                       groupnames = groupnames1, all=FALSE, marginal=FALSE, 
                                       conditional = TRUE,
                                       cond = cond, wb = TRUE)

# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_rod, VP.test2_rod, TabMF_rod, file = 'VPcondwb_ab_Rod.RDATA')


##.......................................................................................................
##  Butterflies
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('bfAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bfAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bfAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bfAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_bf = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_bf = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                      groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                      cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW

VP.test2_bf = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                      groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                      cond = cond, wb = TRUE)


# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_bf, VP.test2_bf, TabMF_bf, file = 'VPcondwb_ab_bf.RDATA')



##.......................................................................................................
##  Moths
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")


load('MothAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_Moth = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_moth = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                        groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                        cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW

VP.test2_moth = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                        groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                        cond = cond, wb = TRUE)


# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_moth, VP.test2_moth, TabMF_Moth, file = 'VPcondwb_ab_moth.RDATA')



##.......................................................................................................
##  Birds
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('BirdAB2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN1 = models
load('BirdAB2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN2 = models
load('BirdAB2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN3 = models
load('BirdAB2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsBirdN4 = models

models = c.Hmsc(modelsBirdN1[[1]], modelsBirdN2[[1]], modelsBirdN3[[1]], modelsBirdN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_bd = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_bd = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                      groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                      cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW
VP.test2_bd = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                      groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                      cond = cond, wb = TRUE)


# Save
##setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_bd, VP.test2_bd, TabMF_bd, file = 'VPcondwb_ab_bd.RDATA')



##.......................................................................................................
##  Winter game
##.......................................................................................................
setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('WintGAB2_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGAB2_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGAB2_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGAB2_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models

models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

# get performances for weighting
preds = computePredictedValues(models, expected = FALSE)
MF = evaluateModelFit(hM=models, predY=preds)

TabMF_wg = cbind(as.data.frame(MF), species = models$spNames)


### grouped VP

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

group = group1
groupnames = groupnames1

# prep conditions .................................

models$studyDesign
models$nr # without bg

catII = models$studyDesign[, c(1,3)]

head(catII)
tail(catII)

cond = as.character(catII$bg)
table(cond)

##.......................................................................................................
### grouped VP - CONDITION

VP.test1_wg = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                      groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                      cond = cond, wb = NULL)


## ............................................................................................
### grouped VP - WITHIN - BTW
VP.test2_wg = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                      groupnames = groupnames1, all=FALSE, marginal=FALSE,  conditional = TRUE,
                                      cond = cond, wb = TRUE)


# Save
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_wg, VP.test2_wg, TabMF_wg, file = 'VPcondwb_ab_wg.RDATA')


## all occurrence
#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")
save(VP.test1_moth, TabMF_Moth, VP.test1_bf, TabMF_bf,
     VP.test1_rod, TabMF_rod, VP.test1_wg, TabMF_wg,
     VP.test1_bd, TabMF_bd, VP.test2_moth, VP.test2_bf,
     VP.test2_rod, VP.test2_wg,
     VP.test2_bd,file = 'VPCPcond_all_ab.RDATA')
