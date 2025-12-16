#
#             Testing VP - CP
#          Get Linear predictor beta*X
#   To calculate correlation between  predictors
#.....................................................................


library(dplyr)
library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used

# function to extract correlation between linear terms
setwd("D:/Helsinki/RECcII/CSC_HMSC")
source('Corr_functions.R')

#
#   Rodents
#
## -------  Occurrence data

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('Rodtest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('Rodtest2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('Rodtest2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('Rodtest2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


modelsrod = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcor_rod = Extract_corr(hM=modelsrod, group = group1,
                     groupnames = groupnames1)

dim(FRcor_rod) # predictor  predict species mcmc samples

#round(apply(FRcor_rod, c(1,2), mean), 3)


## -------  Abundance data

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('RodtestAB_ch1_models_thin_1000_samples_500_chains_1.Rdata')
modelsN1 = models
load('RodtestAB_ch2_models_thin_1000_samples_500_chains_1.Rdata')
modelsN2 = models
load('RodtestAB_ch3_models_thin_1000_samples_500_chains_1.Rdata')
modelsN3 = models
load('RodtestAB_ch4_models_thin_1000_samples_500_chains_1.Rdata')
modelsN4 = models


modelsrodAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       

FRcorAB_rod = Extract_corr(hM=modelsrodAB, group = group1,
                       groupnames = groupnames1)

dim(FRcorAB_rod)



#
#  Bird
#

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('Birdtest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('Birdtest_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('Birdtest_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('Birdtest_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsbd = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcor_bd = Extract_corr(hM=modelsbd, group = group1,
                     groupnames = groupnames1)


dim(FRcor_bd)


setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('BirdAB2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('BirdAB2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('BirdAB2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('BirdAB2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsbdAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcorAB_bd = Extract_corr(hM=modelsbdAB, group = group1,
                       groupnames = groupnames1)

dim(FRcorAB_bd)



#
#  moth
#

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('MothTest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothTest2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothTest2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothTest2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsmoth = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcor_moth = Extract_corr(hM=modelsmoth, group = group1,
                        groupnames = groupnames1)


dim(FRcor_moth)


setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('MothAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsmothAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcorAB_moth = Extract_corr(hM=modelsmothAB, group = group1,
                          groupnames = groupnames1)

dim(FRcorAB_moth)



#
#  bf
#

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('bftest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bftest_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bftest_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bftest_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsbf = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcor_bf = Extract_corr(hM=modelsbf, group = group1,
                          groupnames = groupnames1)


dim(FRcor_bf)


setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('bfAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bfAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bfAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bfAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsbfAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcorAB_bf = Extract_corr(hM=modelsbfAB, group = group1,
                            groupnames = groupnames1)

dim(FRcorAB_bf)



#
#  wg
#

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('WintGtest_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGtest_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGtest_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGtest_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models

modelswg = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcor_wg = Extract_corr(hM=modelswg, group = group1,
                          groupnames = groupnames1)


dim(FRcor_wg)


setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('WintGAB2_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGAB2_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGAB2_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGAB2_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models

modelswgAB = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

##.................................................................................

## grouping
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 4,4, 4,4, 4,4, 4,4, 4,4)
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')       


FRcorAB_wg = Extract_corr(hM=modelswgAB, group = group1,
                            groupnames = groupnames1)

dim(FRcorAB_wg)

## save them all

save(FRcor_bd, FRcor_bf, FRcor_moth, FRcor_rod, FRcor_wg,
     FRcorAB_bd, FRcorAB_bf, FRcorAB_moth, FRcorAB_rod, FRcorAB_wg,
     file = 'FRcorr_taxa1010.RDATA')
