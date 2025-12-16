#
#
#   .....   Sensitivity analysis 2: VP # covar
#
#
#............................................................
#...........................................................

library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)
library(tidyr)
library(dplyr)

setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_ABpure")


load('WintGAB2_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGAB2_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGAB2_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGAB2_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models


models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

#VP1 = computeVariancePartitioning(models)
#plotVariancePartitioning(models, VP1, las=2)

#VP.test1 = computeVarPartSummaries(hM = models, exclude = 2,
#                                          conditional = FALSE, marginal=FALSE)


group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test1g = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                    groupnames = groupnames1, conditional = FALSE, marginal=FALSE)

#...........................................................................
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_sensitivity/Sensitivity_VP")

load('WintGABtest2Ord_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGABtest2Ord_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGABtest2Ord_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGABtest2Ord_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models


modelsOrd = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

#VPord = computeVariancePartitioning(modelsOrd)

#VP.testOrd = computeVarPartSummaries(hM = modelsOrd, exclude = 2,
#                                          conditional = FALSE, marginal=FALSE)



## 
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1

VP.testOrdg = computeVarPartSummaries(hM = modelsOrd, exclude = 2, group = group1,
                                      groupnames = groupnames1, conditional = FALSE, marginal=FALSE)

#-........................................................................



load('WintGABtest3_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN1 = models
load('WintGABtest3_ch2_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN2 = models
load('WintGABtest3_ch3_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN3 = models
load('WintGABtest3_ch4_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelsN4 = models


models2 = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

#VP2 = computeVariancePartitioning(models2)

#VP.test2 = computeVarPartSummaries(hM = models2, exclude = 2,
#                                            conditional = FALSE, marginal=FALSE)


## 
group1 = c(1,1, 2,2, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test2g = computeVarPartSummaries(hM = models2, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

##
VP.test1g
VP.testOrdg
VP.test2g

Mean.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), mean)
Mean.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), mean)
Mean.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), mean)

Meanwg = as.data.frame(rbind(Mean.vp.hmsc1[1:2], Mean.vp.hmscOrd[1:2], Mean.vp.hmsc2[1:2]))
colnames(Meanwg)[1:2] = c('Climate', 'Habitat')
Meanwg$model = c('Normal', 'ordination', 'Proportion grouping')
Meanwg$model = factor(Meanwg$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Meanwglg <- Meanwg %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Meanwglg$Taxa = 'large mammals'





Up.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)

Upwg = as.data.frame(rbind(Up.vp.hmsc1[1:2], Up.vp.hmscOrd[1:2], Up.vp.hmsc2[1:2]))
colnames(Upwg)[1:2] = c('Climate', 'Habitat')
Upwg$model = c('Normal', 'ordination', 'Proportion grouping')

Upwg$model = factor(Upwg$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Upwglg <- Upwg %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Upwglg$Taxa = 'large mammals'





Low.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)

Lowwg = as.data.frame(rbind(Low.vp.hmsc1[1:2], Low.vp.hmscOrd[1:2], Low.vp.hmsc2[1:2]))
colnames(Lowwg)[1:2] = c('Climate', 'Habitat')
Lowwg$model = c('Normal', 'ordination', 'Proportion grouping')

Lowwg$model = factor(Lowwg$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Lowwglg <- Lowwg %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Lowwglg$Taxa = 'large mammals'




### butterflies
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('bfAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bfAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bfAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bfAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])




group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test1g = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#...........................................................................

setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_sensitivity/Sensitivity_VP")

load('bfABtestOrd_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bfABtestOrd_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bfABtestOrd_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bfABtestOrd_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


modelsOrd = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])


## 
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1

VP.testOrdg = computeVarPartSummaries(hM = modelsOrd, exclude = 2, group = group1,
                                      groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#-........................................................................



load('bfABtest3_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('bfABtest3_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('bfABtest3_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('bfABtest3_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models2 =  c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

## 
group1 = c(1,1, 2,2, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test2g = computeVarPartSummaries(hM = models2, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

##
VP.test1g
VP.testOrdg
VP.test2g


Mean.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), mean)
Mean.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), mean)
Mean.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), mean)

Meanbf = as.data.frame(rbind(Mean.vp.hmsc1[1:2], Mean.vp.hmscOrd[1:2], Mean.vp.hmsc2[1:2]))
colnames(Meanbf)[1:2] = c('Climate', 'Habitat')
Meanbf$model = c('Normal', 'ordination', 'Proportion grouping')

Meanbf$model = factor(Meanbf$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Meanbflg <- Meanbf %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Meanbflg$Taxa = 'butterflies'






Up.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)

Upbf = as.data.frame(rbind(Up.vp.hmsc1[1:2], Up.vp.hmscOrd[1:2], Up.vp.hmsc2[1:2]))
colnames(Upbf)[1:2] = c('Climate', 'Habitat')
Upbf$model = c('Normal', 'ordination', 'Proportion grouping')

Upbf$model = factor(Upbf$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Upbflg <- Upbf %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Upbflg$Taxa = 'butterflies'





Low.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)

Lowbf = as.data.frame(rbind(Low.vp.hmsc1[1:2], Low.vp.hmscOrd[1:2], Low.vp.hmsc2[1:2]))
colnames(Lowbf)[1:2] = c('Climate', 'Habitat')
Lowbf$model = c('Normal', 'ordination', 'Proportion grouping')

Lowbf$model = factor(Lowbf$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Lowbflg <- Lowbf %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Lowbflg$Taxa = 'butterflies'





### Moths
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('MothAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])


group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test1g = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#...........................................................................
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_sensitivity/Sensitivity_VP")

load('MothAB_test2Ord_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothAB_test2Ord_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothAB_test2Ord_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothAB_test2Ord_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


modelsOrd = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

#VPord = computeVariancePartitioning(modelsOrd)


VP.testOrd = computeVarPartSummaries(hM = modelsOrd, exclude = 2,
                                     conditional = FALSE, marginal=FALSE)



## 
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1

VP.testOrdg = computeVarPartSummaries(hM = modelsOrd, exclude = 2, group = group1,
                                      groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#-........................................................................



load('MothAB_test3_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('MothAB_test3_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('MothAB_test3_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('MothAB_test3_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models2 =  c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

## 
group1 = c(1,1, 2,2, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test2g = computeVarPartSummaries(hM = models2, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

##
VP.test1g
VP.testOrdg
VP.test2g


Mean.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), mean)
Mean.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), mean)
Mean.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), mean)

MeanMoth = as.data.frame(rbind(Mean.vp.hmsc1[1:2], Mean.vp.hmscOrd[1:2], Mean.vp.hmsc2[1:2]))
colnames(MeanMoth)[1:2] = c('Climate', 'Habitat')
MeanMoth$model = c('Normal', 'ordination', 'Proportion grouping')

MeanMoth$model = factor(MeanMoth$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

MeanMothlg <- MeanMoth %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

MeanMothlg$Taxa = 'moths'



Up.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)

UpMoth = as.data.frame(rbind(Up.vp.hmsc1[1:2], Up.vp.hmscOrd[1:2], Up.vp.hmsc2[1:2]))
colnames(UpMoth)[1:2] = c('Climate', 'Habitat')
UpMoth$model = c('Normal', 'ordination', 'Proportion grouping')

UpMoth$model = factor(UpMoth$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

UpMothlg <- UpMoth %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

UpMothlg$Taxa = 'moths'





Low.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)

LowMoth = as.data.frame(rbind(Low.vp.hmsc1[1:2], Low.vp.hmscOrd[1:2], Low.vp.hmsc2[1:2]))
colnames(LowMoth)[1:2] = c('Climate', 'Habitat')
LowMoth$model = c('Normal', 'ordination', 'Proportion grouping')

LowMoth$model = factor(LowMoth$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

LowMothlg <- LowMoth %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

LowMothlg$Taxa = 'moths'



### Birds
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('BirdAB2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('BirdAB2_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('BirdAB2_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('BirdAB2_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test1g = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#...........................................................................
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_sensitivity/Sensitivity_VP")

load('BirdAB_test2Ord_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('BirdAB_test2Ord_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('BirdAB_test2Ord_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('BirdAB_test2Ord_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


modelsOrd = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])


## 
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1

VP.testOrdg = computeVarPartSummaries(hM = modelsOrd, exclude = 2, group = group1,
                                      groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#-........................................................................



load('BirdAB_test3_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('BirdAB_test3_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('BirdAB_test3_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('BirdAB_test3_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models2 =  c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

group1 = c(1,1, 2,2, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test2g = computeVarPartSummaries(hM = models2, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

##
VP.test1g
VP.testOrdg
VP.test2g


Mean.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), mean)
Mean.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), mean)
Mean.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), mean)

MeanBird = as.data.frame(rbind(Mean.vp.hmsc1[1:2], Mean.vp.hmscOrd[1:2], Mean.vp.hmsc2[1:2]))
colnames(MeanBird)[1:2] = c('Climate', 'Habitat')
MeanBird$model = c('Normal', 'ordination', 'Proportion grouping')

MeanBird$model = factor(MeanBird$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

MeanBirdlg <- MeanBird %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

MeanBirdlg$Taxa = 'birds'




Up.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)

UpBird = as.data.frame(rbind(Up.vp.hmsc1[1:2], Up.vp.hmscOrd[1:2], Up.vp.hmsc2[1:2]))
colnames(UpBird)[1:2] = c('Climate', 'Habitat')
UpBird$model = c('Normal', 'ordination', 'Proportion grouping')

UpBird$model = factor(UpBird$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

UpBirdlg <- UpBird %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

UpBirdlg$Taxa = 'birds'





Low.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)

LowBird = as.data.frame(rbind(Low.vp.hmsc1[1:2], Low.vp.hmscOrd[1:2], Low.vp.hmsc2[1:2]))
colnames(LowBird)[1:2] = c('Climate', 'Habitat')
LowBird$model = c('Normal', 'ordination', 'Proportion grouping')

LowBird$model = factor(LowBird$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

LowBirdlg <- LowBird %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

LowBirdlg$Taxa = 'birds'


## Rodents
#..................................................

setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_ABpure")

load('RodtestAB_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('RodtestAB_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('RodtestAB_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('RodtestAB_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models


models = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test1g = computeVarPartSummaries(hM = models, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#...........................................................................
setwd("/home/local/guilbaul/tsclient/D/Helsinki/RECcII_project/CSC_HMSC/fullmodelABpure/Results_sensitivity/Sensitivity_VP")


load('RodtestABtest2Ord_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('RodtestABtest2Ord_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('RodtestABtest2Ord_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('RodtestABtest2Ord_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

modelsOrd = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])

VP.testOrd = computeVarPartSummaries(hM = modelsOrd, exclude = 2,
                                     conditional = FALSE, marginal=FALSE)



## 
group1 = c(1,1, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1

VP.testOrdg = computeVarPartSummaries(hM = modelsOrd, exclude = 2, group = group1,
                                      groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

#-........................................................................


load('RodtestABtest3_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelsN1 = models
load('RodtestABtest3_ch2_models_thin_1000_samples_250_chains_1.Rdata')
modelsN2 = models
load('RodtestABtest3_ch3_models_thin_1000_samples_250_chains_1.Rdata')
modelsN3 = models
load('RodtestABtest3_ch4_models_thin_1000_samples_250_chains_1.Rdata')
modelsN4 = models

models2 = c.Hmsc(modelsN1[[1]], modelsN2[[1]], modelsN3[[1]], modelsN4[[1]])


VP.test2 = computeVarPartSummaries(hM = models2, exclude = 2,
                                   conditional = FALSE, marginal=FALSE)


## 
group1 = c(1,1, 2,2, 2,2, 2,2, 2,2, 3,3, 3,3, 3,3)
groupnames1=c("Effort", "Climate", "Habitat")       

group = group1
groupnames = groupnames1


VP.test2g = computeVarPartSummaries(hM = models2, exclude = 2, group = group1,
                                    groupnames = groupnames1,conditional = FALSE, marginal=FALSE)

##
VP.test1g
VP.testOrdg
# VP.test2g
# 
Mean.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), mean)
Mean.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), mean)
Mean.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), mean)

Meanrod = as.data.frame(rbind(Mean.vp.hmsc1[1:2], Mean.vp.hmscOrd[1:2], Mean.vp.hmsc2[1:2]))
colnames(Meanrod)[1:2] = c('Climate', 'Habitat')
Meanrod$model = c('Normal', 'ordination', 'Proportion grouping')

Meanrod$model = factor(Meanrod$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Meanrodlg <- Meanrod %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Meanrodlg$Taxa = 'small mammals'


Up.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)
Up.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.975, na.rm=T)

Uprod = as.data.frame(rbind(Up.vp.hmsc1[1:2], Up.vp.hmscOrd[1:2], Up.vp.hmsc2[1:2]))
colnames(Uprod)[1:2] = c('Climate', 'Habitat')
Uprod$model = c('Normal', 'ordination', 'Proportion grouping')

Uprod$model = factor(Uprod$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Uprodlg <- Uprod %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Uprodlg$Taxa = 'small mammals'





Low.vp.hmsc1 = apply(VP.test1g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmsc2 = apply(VP.test2g$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)
Low.vp.hmscOrd = apply(VP.testOrdg$Vdiag, c(1,2), quantile, probs=0.025, na.rm=T)

Lowrod = as.data.frame(rbind(Low.vp.hmsc1[1:2], Low.vp.hmscOrd[1:2], Low.vp.hmsc2[1:2]))
colnames(Lowrod)[1:2] = c('Climate', 'Habitat')
Lowrod$model = c('Normal', 'ordination', 'Proportion grouping')

Lowrod$model = factor(Lowrod$model, levels = c('Normal', 'ordination', 'Proportion grouping'))

Lowrodlg <- Lowrod %>% 
  pivot_longer(
    cols = 'Climate':'Habitat', 
    names_to = "Driver",
    values_to = "VP"
  )

Lowrodlg$Taxa = 'small mammals'


###  combine all

Meanlg = rbind(Meanrodlg, Meanbflg, MeanBirdlg, MeanMothlg, Meanwglg)
Meanlg$Taxa = factor(Meanlg$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
Meanlg$stat = 'Mean'

Uplg = rbind(Uprodlg, Upbflg, UpBirdlg, UpMothlg, Upwglg)
Uplg$Taxa = factor(Uplg$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
Uplg$stat = 'Up'

Lowlg = rbind(Lowrodlg, Lowbflg, LowBirdlg, LowMothlg, Lowwglg)
Lowlg$Taxa = factor(Lowlg$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
Lowlg$stat = 'Low'

## New plots
library(ggplot2)

Meanlg %>%mutate( 
  VP = ifelse(Driver=="Climate", VP*(-1), 
              VP*1))%>% 
  ggplot(aes(x = model, y = VP,fill = Driver)) +  
  geom_bar(stat = "identity") + 
  facet_wrap(.~ Taxa, ncol=1) +
  scale_fill_manual(values = c("#66C2A5", 'orange3')) +
  scale_y_continuous(limits = c(-0.6, 0.6)) + coord_flip()+ 
  labs(title = "", x = "Model",  
       y = "VP") + theme_classic() +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) #change font size of legend title  


Alllg = rbind(Meanlg, Uplg, Lowlg)

levels(Alllg$model)[1] = 'Default'
levels(Meanlg$model)[1] = 'Default'
levels(Uplg$model)[1] = 'Default'
levels(Lowlg$model)[1] = 'Default'

levels(Meanlg$model)[2] = 'Ordination'
levels(Uplg$model)[2] = 'Ordination'
levels(Lowlg$model)[2] = 'Ordination'

library(ggplot2)

ggplot()+  
  # Mean
  geom_col(
    data = Meanlg %>% mutate( 
      VP = ifelse(Driver=="Climate", VP*(-1), 
                  VP*1)),
    mapping = aes(
      x = model,
      y = VP,
      fill = Driver),
    #colour = "black",                               # black color around bars
    alpha = 1,                                    # more transparent
    width = 0.6)+                                     # full width
  
  # up
  geom_col(
    data = Uplg %>% mutate( 
      VP = ifelse(Driver=="Climate", VP*(-1), 
                  VP*1)),
    mapping = aes(
      x = model,
      y = VP,
      fill = Driver),                           # fill of bars by gender
    colour = "black",                               # black color around bars
    alpha = 0.3,                                      # not transparent 
    width = 1)+                                   # half width
  # low
  geom_col(
    data = Lowlg %>% mutate( 
      VP = ifelse(Driver=="Climate", VP*(-1), 
                  VP*1)),
    mapping = aes(
      x = model,
      y = VP,
      fill = Driver),                           # fill of bars by gender
    colour = "black",                               # black color around bars
    alpha = 0.2,                                      # not transparent 
    width = 1)+                                   # half width
  
  facet_wrap(.~ Taxa, ncol=1) +
  scale_fill_manual(values = c("#66C2A5", 'orange3')) +
  scale_y_continuous(limits = c(-1.2, 1.2)) + coord_flip()+ 
  labs(title = "", x = "Model",  
       y = "Variance Partition") + theme_classic() +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) #change font size of legend title  
