#
#             Summary plots to compare taxa
#
#
# --------------------------------------------------------------------


library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)

setwd("~/Data/4.VPanalyses/VPCP") # to adapt depending on where it is saved

#
#   Occurrence
#

load('VPCP_rod.RDATA')
dim(VP.test3_rod$Vnorm)
dim(VP.test3_rod$Cnorm)

VP.test3.rod = array(NA, dim = dim(VP.test3_rod$Cnorm))

for (i in 1:dim(VP.test3_rod$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_rod$Vnorm)[3]){  ## species
    
    VP.test3.rod[,,j,i] = VP.test3_rod$Cnorm[,,j,i]
    diag(VP.test3.rod[,,j,i]) = VP.test3_rod$Vnorm[,,j,i]
    
    colnames( VP.test3.rod[,,j,i]) = colnames(VP.test3_rod$Cnorm[,,j,i])
    rownames( VP.test3.rod[,,j,i]) = rownames(VP.test3_rod$Cnorm[,,j,i])
    
  }
}


varr.rod <- aperm(array(TabMF_rod$TjurR2, dim = c(8L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.rod)

Mean.fix.sprod = apply(varr.rod*VP.test3.rod[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)



load('VPCP_bf.RDATA')
dim(VP.test3_bf$Vnorm)
dim(VP.test3_bf$Cnorm)

VP.test3.bf = array(NA, dim = dim(VP.test3_bf$Cnorm))

for (i in 1:dim(VP.test3_bf$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bf$Vnorm)[3]){  ## species
    
    VP.test3.bf[,,j,i] = VP.test3_bf$Cnorm[,,j,i]
    diag(VP.test3.bf[,,j,i]) = VP.test3_bf$Vnorm[,,j,i]
    
    colnames( VP.test3.bf[,,j,i]) = colnames(VP.test3_bf$Cnorm[,,j,i])
    rownames( VP.test3.bf[,,j,i]) = rownames(VP.test3_bf$Cnorm[,,j,i])
    
  }
}


varr.bf <- aperm(array(TabMF_bf$TjurR2, dim = c(57L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.bf)

Mean.fix.spbf = apply(varr.bf*VP.test3.bf[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)



load('VPCP_wg.RDATA')
dim(VP.test3_wg$Vnorm)
dim(VP.test3_wg$Cnorm)

VP.test3.wg = array(NA, dim = dim(VP.test3_wg$Cnorm))

for (i in 1:dim(VP.test3_wg$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_wg$Vnorm)[3]){  ## species
    
    VP.test3.wg[,,j,i] = VP.test3_wg$Cnorm[,,j,i]
    diag(VP.test3.wg[,,j,i]) = VP.test3_wg$Vnorm[,,j,i]
    
    colnames( VP.test3.wg[,,j,i]) = colnames(VP.test3_wg$Cnorm[,,j,i])
    rownames( VP.test3.wg[,,j,i]) = rownames(VP.test3_wg$Cnorm[,,j,i])
    
  }
}


varr.wg <- aperm(array(TabMF_wg$TjurR2, dim = c(15L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.wg)

Mean.fix.spwg = apply(varr.wg*VP.test3.wg[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)


load('VPCP_bd.RDATA')
dim(VP.test3_bd$Vnorm)
dim(VP.test3_bd$Cnorm)

VP.test3.bd = array(NA, dim = dim(VP.test3_bd$Cnorm))

for (i in 1:dim(VP.test3_bd$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bd$Vnorm)[3]){  ## species
    
    VP.test3.bd[,,j,i] = VP.test3_bd$Cnorm[,,j,i]
    diag(VP.test3.bd[,,j,i]) = VP.test3_bd$Vnorm[,,j,i]
    
    colnames( VP.test3.bd[,,j,i]) = colnames(VP.test3_bd$Cnorm[,,j,i])
    rownames( VP.test3.bd[,,j,i]) = rownames(VP.test3_bd$Cnorm[,,j,i])
    
  }
}


varr.bd <- aperm(array(TabMF_bd$TjurR2,dim = c(102L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))

Mean.fix.spbd = apply(varr.bd*VP.test3.bd[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)



load('VPCP_moth.RDATA')
dim(VP.test3_moth$Vnorm)
dim(VP.test3_moth$Cnorm)

VP.test3.moth = array(NA, dim = dim(VP.test3_moth$Cnorm))

for (i in 1:dim(VP.test3_moth$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_moth$Vnorm)[3]){  ## species
    
    VP.test3.moth[,,j,i] = VP.test3_moth$Cnorm[,,j,i]
    diag(VP.test3.moth[,,j,i]) = VP.test3_moth$Vnorm[,,j,i]
    
    colnames( VP.test3.moth[,,j,i]) = colnames(VP.test3_moth$Cnorm[,,j,i])
    rownames( VP.test3.moth[,,j,i]) = rownames(VP.test3_moth$Cnorm[,,j,i])
    
  }
}


varr.moth <- aperm(array(TabMF_Moth$TjurR2, dim = c(319L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.moth)


Mean.fix.spmoth = apply(varr.moth*VP.test3.moth[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)





####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

Mean.fix.sprod = Mean.fix.sprod*100
Mean.fix.spbf = Mean.fix.spbf*100
Mean.fix.spbd = Mean.fix.spbd*100
Mean.fix.spmoth = Mean.fix.spmoth*100
Mean.fix.spwg = Mean.fix.spwg*100


VP.sprod = data.frame(Climate = Mean.fix.sprod[1,1,], habprop = Mean.fix.sprod[3,3,], 
                      habConf = Mean.fix.sprod[2,2,], Site = Mean.fix.sprod[4,4,], 
                      year = Mean.fix.sprod[5,5,], region = Mean.fix.sprod[6,6,])

VP.sprodlg <- VP.sprod %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.sprodlg$taxa = 'small mammals'

VP.spmoth = data.frame(Climate = Mean.fix.spmoth[1,1,], habprop = Mean.fix.spmoth[3,3,], 
                      habConf = Mean.fix.spmoth[2,2,], Site = Mean.fix.spmoth[4,4,], 
                      year = Mean.fix.spmoth[5,5,], region = Mean.fix.spmoth[6,6,])

VP.spmothlg <- VP.spmoth %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spmothlg$taxa = 'moths'

VP.spbd = data.frame(Climate = Mean.fix.spbd[1,1,], habprop = Mean.fix.spbd[3,3,], 
                     habConf = Mean.fix.spbd[2,2,], Site = Mean.fix.spbd[4,4,], 
                     year = Mean.fix.spbd[5,5,], region = Mean.fix.spbd[6,6,])


VP.spbdlg <- VP.spbd %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbdlg$taxa = 'birds'

VP.spbf = data.frame(Climate = Mean.fix.spbf[1,1,], habprop = Mean.fix.spbf[3,3,], 
                     habConf = Mean.fix.spbf[2,2,], Site = Mean.fix.spbf[4,4,], 
                     year = Mean.fix.spbf[5,5,], region = Mean.fix.spbf[6,6,])


VP.spbflg <- VP.spbf %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbflg$taxa = 'butterflies'

VP.spwg = data.frame(Climate = Mean.fix.spwg[1,1,], habprop = Mean.fix.spwg[3,3,], 
                     habConf = Mean.fix.spwg[2,2,], Site = Mean.fix.spwg[4,4,], 
                     year = Mean.fix.spwg[5,5,], region = Mean.fix.spwg[6,6,])

VP.spwglg <- VP.spwg %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spwglg$taxa = 'large mammals'


VP.mcmc = rbind(VP.sprodlg, VP.spbflg, VP.spmothlg, VP.spbdlg, VP.spwglg)
VP.mcmc$VP = factor(VP.mcmc$VP, levels = c('Climate', 'habprop', 'habConf',
                                           'Site', 'year', 'region'))


ggplot(data = VP.mcmc, aes(x=value, y = taxa, group = taxa, fill = VP)) +
  geom_density_ridges(alpha=0.8) +
  scale_fill_manual(values = c("#66C2A5", 'orange3', "yellow3", "grey78", "blue4", 'green4')) +
  facet_wrap(.~ VP) + xlim(c(-15, 60)) +
  theme_classic() 


#### iqr


resrodmc.CIF.uppCI = apply(varr.rod*VP.test3.rod, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resrodmc.CIF.lowCI = apply(varr.rod*VP.test3.rod, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

resmothmc.CIF.uppCI = apply(varr.moth*VP.test3.moth, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resmothmc.CIF.lowCI = apply(varr.moth*VP.test3.moth, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

resbfmc.CIF.uppCI = apply(varr.bf*VP.test3.bf, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resbfmc.CIF.lowCI = apply(varr.bf*VP.test3.bf, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

reswgmc.CIF.uppCI = apply(varr.wg*VP.test3.wg, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
reswgmc.CIF.lowCI = apply(varr.wg*VP.test3.wg, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

resbdmc.CIF.uppCI = apply(varr.bd*VP.test3.bd, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resbdmc.CIF.lowCI = apply(varr.bd*VP.test3.bd, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

IQR_rod = 100*(resrodmc.CIF.uppCI - resrodmc.CIF.lowCI)
IQR_bf= 100*(resbfmc.CIF.uppCI - resbfmc.CIF.lowCI)
IQR_moth= 100*(resmothmc.CIF.uppCI - resmothmc.CIF.lowCI)
IQR_wg = 100*(reswgmc.CIF.uppCI - reswgmc.CIF.lowCI)
IQR_bd = 100*(resbdmc.CIF.uppCI - resbdmc.CIF.lowCI)



VP.sprod = data.frame(Climate = IQR_rod[1,1,], habprop = IQR_rod[3,3,], 
                      habConf = IQR_rod[2,2,], Site = IQR_rod[4,4,], 
                      year = IQR_rod[5,5,], region = IQR_rod[6,6,])

VP.sprodlg <- VP.sprod %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.sprodlg$taxa = 'small mammals'

VP.spbd = data.frame(Climate = IQR_bd[1,1,], habprop = IQR_bd[3,3,], 
                     habConf = IQR_bd[2,2,], Site = IQR_bd[4,4,], 
                     year = IQR_bd[5,5,], region = IQR_bd[6,6,])

VP.spbdlg <- VP.spbd %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbdlg$taxa = 'birds'

VP.spbf = data.frame(Climate = IQR_bf[1,1,], habprop = IQR_bf[3,3,], 
                     habConf = IQR_bf[2,2,], Site = IQR_bf[4,4,], 
                     year = IQR_bf[5,5,], region = IQR_bf[6,6,])

VP.spbflg <- VP.spbf %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbflg$taxa = 'butterflies'


VP.spmoth = data.frame(Climate = IQR_moth[1,1,], habprop = IQR_moth[3,3,], 
                     habConf = IQR_moth[2,2,], Site = IQR_moth[4,4,], 
                     year = IQR_moth[5,5,], region = IQR_moth[6,6,])

VP.spmothlg <- VP.spmoth %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spmothlg$taxa = 'moths'


VP.spwg = data.frame(Climate = IQR_wg[1,1,], habprop = IQR_wg[3,3,], 
                     habConf = IQR_wg[2,2,], Site = IQR_wg[4,4,], 
                     year = IQR_wg[5,5,], region = IQR_wg[6,6,])

VP.spwglg <- VP.spwg %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spwglg$taxa = 'large mammals'


VPiqr.sp = rbind(VP.sprodlg, VP.spbflg, VP.spbdlg, VP.spwglg, VP.spmothlg)
VPiqr.sp$VP = factor(VPiqr.sp$VP, levels = c('Climate', 'habprop', 'habConf',
                                       'Site', 'year', 'region'))


ggplot(data = VPiqr.sp, aes(x=value, y = taxa, group = taxa, fill = VP)) +
  geom_density_ridges(alpha=0.8) +
  scale_fill_manual(values = c("#66C2A5", 'orange3', "yellow3", "grey78", "blue4", 'green4')) +
  facet_wrap(.~ VP) + xlim(c(-10, 130)) +
  theme_classic() 


## both
##### mcmc

## both

head(VP.mcmc)

VPiqr.mcmc = VPiqr.sp
head(VPiqr.sp)

VP.mcmc$Measure = 'Mean'
VPiqr.mcmc$Measure = 'IQR'

VPmeasures.mcmc = rbind(VP.mcmc, VPiqr.mcmc)
VPmeasures.mcmc$Measure = factor(VPmeasures.mcmc$Measure, levels = c('Mean', 'IQR'))
VPmeasures.mcmc$taxa = factor(VPmeasures.mcmc$taxa, levels = rev(c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals')))


levels(VPmeasures.mcmc$VP)[c(2:3, 5:6)] = c("Habcomp", "LandConf", 'Year', 'Region')

ggplot(data = VPmeasures.mcmc, aes(x=value, y = taxa, group = taxa, fill = VP)) +
  geom_density_ridges(alpha=0.8, quantile_lines = TRUE, quantile=2) +
  scale_fill_manual(values = c("#66C2A5", 'orange3', "yellow3", "grey78", "blue4", 'green4')) +
  facet_grid(VP~ Measure) + xlim(c(-5, 75)) +
  labs(title = "", x = "Variance partition",  
       y = "Taxa") + theme_classic() +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) #change font size of legend title  




#
#   Abundance
#

load('VPCP_ab_rod.RDATA')
dim(VP.test3_rod$Vnorm)
dim(VP.test3_rod$Cnorm)

VP.test3.rod = array(NA, dim = dim(VP.test3_rod$Cnorm))

for (i in 1:dim(VP.test3_rod$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_rod$Vnorm)[3]){  ## species
    
    VP.test3.rod[,,j,i] = VP.test3_rod$Cnorm[,,j,i]
    diag(VP.test3.rod[,,j,i]) = VP.test3_rod$Vnorm[,,j,i]
    
    colnames( VP.test3.rod[,,j,i]) = colnames(VP.test3_rod$Cnorm[,,j,i])
    rownames( VP.test3.rod[,,j,i]) = rownames(VP.test3_rod$Cnorm[,,j,i])
    
  }
}


varr.rod <- aperm(array(TabMF_rod$SR2,dim = c(8L, 6L, 6L, 2000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.rod)

Mean.fix.sprod = apply(varr.rod*VP.test3.rod[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)



load('VPCP_ab_bf.RDATA')
dim(VP.test3_bf$Vnorm)
dim(VP.test3_bf$Cnorm)

VP.test3.bf = array(NA, dim = dim(VP.test3_bf$Cnorm))

for (i in 1:dim(VP.test3_bf$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bf$Vnorm)[3]){  ## species
    
    VP.test3.bf[,,j,i] = VP.test3_bf$Cnorm[,,j,i]
    diag(VP.test3.bf[,,j,i]) = VP.test3_bf$Vnorm[,,j,i]
    
    colnames( VP.test3.bf[,,j,i]) = colnames(VP.test3_bf$Cnorm[,,j,i])
    rownames( VP.test3.bf[,,j,i]) = rownames(VP.test3_bf$Cnorm[,,j,i])
    
  }
}


varr.bf <- aperm(array(TabMF_bf$SR2, dim = c(57L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.bf)

Mean.fix.spbf = apply(varr.bf*VP.test3.bf[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)



load('VPCP_ab_wg.RDATA')
dim(VP.test3_wg$Vnorm)
dim(VP.test3_wg$Cnorm)

VP.test3.wg = array(NA, dim = dim(VP.test3_wg$Cnorm))

for (i in 1:dim(VP.test3_wg$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_wg$Vnorm)[3]){  ## species
    
    VP.test3.wg[,,j,i] = VP.test3_wg$Cnorm[,,j,i]
    diag(VP.test3.wg[,,j,i]) = VP.test3_wg$Vnorm[,,j,i]
    
    colnames( VP.test3.wg[,,j,i]) = colnames(VP.test3_wg$Cnorm[,,j,i])
    rownames( VP.test3.wg[,,j,i]) = rownames(VP.test3_wg$Cnorm[,,j,i])
    
  }
}


varr.wg <- aperm(array(TabMF_wg$SR2, dim = c(15L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.wg)

Mean.fix.spwg = apply(varr.wg*VP.test3.wg[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)


load('VPCP_ab_bd.RDATA')
dim(VP.test3_bd$Vnorm)
dim(VP.test3_bd$Cnorm)

VP.test3.bd = array(NA, dim = dim(VP.test3_bd$Cnorm))

for (i in 1:dim(VP.test3_bd$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bd$Vnorm)[3]){  ## species
    
    VP.test3.bd[,,j,i] = VP.test3_bd$Cnorm[,,j,i]
    diag(VP.test3.bd[,,j,i]) = VP.test3_bd$Vnorm[,,j,i]
    
    colnames( VP.test3.bd[,,j,i]) = colnames(VP.test3_bd$Cnorm[,,j,i])
    rownames( VP.test3.bd[,,j,i]) = rownames(VP.test3_bd$Cnorm[,,j,i])
    
  }
}


varr.bd <- aperm(array(TabMF_bd$SR2, dim = c(102L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))

Mean.fix.spbd = apply(varr.bd*VP.test3.bd[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)



load('VPCP_ab_moth.RDATA')
dim(VP.test3_moth$Vnorm)
dim(VP.test3_moth$Cnorm)

VP.test3.moth = array(NA, dim = dim(VP.test3_moth$Cnorm))

for (i in 1:dim(VP.test3_moth$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_moth$Vnorm)[3]){  ## species
    
    VP.test3.moth[,,j,i] = VP.test3_moth$Cnorm[,,j,i]
    diag(VP.test3.moth[,,j,i]) = VP.test3_moth$Vnorm[,,j,i]
    
    colnames( VP.test3.moth[,,j,i]) = colnames(VP.test3_moth$Cnorm[,,j,i])
    rownames( VP.test3.moth[,,j,i]) = rownames(VP.test3_moth$Cnorm[,,j,i])
    
  }
}


varr.moth <- aperm(array(TabMF_Moth$SR2,dim = c(319L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.moth)


Mean.fix.spmoth = apply(varr.moth*VP.test3.moth[1:6, 1:6, ,], c(1, 2, 4), mean, na.rm=T)





####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

Mean.fix.sprod = Mean.fix.sprod*100
Mean.fix.spbf = Mean.fix.spbf*100
Mean.fix.spbd = Mean.fix.spbd*100
Mean.fix.spmoth = Mean.fix.spmoth*100
Mean.fix.spwg = Mean.fix.spwg*100


VP.sprod = data.frame(Climate = Mean.fix.sprod[1,1,], habprop = Mean.fix.sprod[3,3,], 
                      habConf = Mean.fix.sprod[2,2,], Site = Mean.fix.sprod[4,4,], 
                      year = Mean.fix.sprod[5,5,], region = Mean.fix.sprod[6,6,])

VP.sprodlg <- VP.sprod %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.sprodlg$taxa = 'small mammals'

VP.spmoth = data.frame(Climate = Mean.fix.spmoth[1,1,], habprop = Mean.fix.spmoth[3,3,], 
                       habConf = Mean.fix.spmoth[2,2,], Site = Mean.fix.spmoth[4,4,], 
                       year = Mean.fix.spmoth[5,5,], region = Mean.fix.spmoth[6,6,])

VP.spmothlg <- VP.spmoth %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spmothlg$taxa = 'moths'

VP.spbd = data.frame(Climate = Mean.fix.spbd[1,1,], habprop = Mean.fix.spbd[3,3,], 
                     habConf = Mean.fix.spbd[2,2,], Site = Mean.fix.spbd[4,4,], 
                     year = Mean.fix.spbd[5,5,], region = Mean.fix.spbd[6,6,])


VP.spbdlg <- VP.spbd %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbdlg$taxa = 'birds'

VP.spbf = data.frame(Climate = Mean.fix.spbf[1,1,], habprop = Mean.fix.spbf[3,3,], 
                     habConf = Mean.fix.spbf[2,2,], Site = Mean.fix.spbf[4,4,], 
                     year = Mean.fix.spbf[5,5,], region = Mean.fix.spbf[6,6,])


VP.spbflg <- VP.spbf %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbflg$taxa = 'butterflies'

VP.spwg = data.frame(Climate = Mean.fix.spwg[1,1,], habprop = Mean.fix.spwg[3,3,], 
                     habConf = Mean.fix.spwg[2,2,], Site = Mean.fix.spwg[4,4,], 
                     year = Mean.fix.spwg[5,5,], region = Mean.fix.spwg[6,6,])

VP.spwglg <- VP.spwg %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spwglg$taxa = 'large mammals'


VP.mcmc = rbind(VP.sprodlg, VP.spbflg, VP.spmothlg, VP.spbdlg, VP.spwglg)
VP.mcmc$VP = factor(VP.mcmc$VP, levels = c('Climate', 'habprop', 'habConf',
                                           'Site', 'year', 'region'))


ggplot(data = VP.mcmc, aes(x=value, y = taxa, group = taxa, fill = VP)) +
  geom_density_ridges(alpha=0.8) +
  scale_fill_manual(values = c("#66C2A5", 'orange3', "yellow3", "grey78", "blue4", 'green4')) +
  facet_wrap(.~ VP) + xlim(c(-15, 60)) +
  theme_classic() 


#### iqr


resrodmc.CIF.uppCI = apply(varr.rod*VP.test3.rod, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resrodmc.CIF.lowCI = apply(varr.rod*VP.test3.rod, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

resmothmc.CIF.uppCI = apply(varr.moth*VP.test3.moth, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resmothmc.CIF.lowCI = apply(varr.moth*VP.test3.moth, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

resbfmc.CIF.uppCI = apply(varr.bf*VP.test3.bf, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resbfmc.CIF.lowCI = apply(varr.bf*VP.test3.bf, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

reswgmc.CIF.uppCI = apply(varr.wg*VP.test3.wg, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
reswgmc.CIF.lowCI = apply(varr.wg*VP.test3.wg, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

resbdmc.CIF.uppCI = apply(varr.bd*VP.test3.bd, c(1, 2, 4), quantile, probs=0.975, na.rm=T)
resbdmc.CIF.lowCI = apply(varr.bd*VP.test3.bd, c(1, 2, 4), quantile, probs=0.025, na.rm=T)

IQR_rod = 100*(resrodmc.CIF.uppCI - resrodmc.CIF.lowCI)
IQR_bf= 100*(resbfmc.CIF.uppCI - resbfmc.CIF.lowCI)
IQR_moth= 100*(resmothmc.CIF.uppCI - resmothmc.CIF.lowCI)
IQR_wg = 100*(reswgmc.CIF.uppCI - reswgmc.CIF.lowCI)
IQR_bd = 100*(resbdmc.CIF.uppCI - resbdmc.CIF.lowCI)



VP.sprod = data.frame(Climate = IQR_rod[1,1,], habprop = IQR_rod[3,3,], 
                      habConf = IQR_rod[2,2,], Site = IQR_rod[4,4,], 
                      year = IQR_rod[5,5,], region = IQR_rod[6,6,])

VP.sprodlg <- VP.sprod %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.sprodlg$taxa = 'small mammals'

VP.spbd = data.frame(Climate = IQR_bd[1,1,], habprop = IQR_bd[3,3,], 
                     habConf = IQR_bd[2,2,], Site = IQR_bd[4,4,], 
                     year = IQR_bd[5,5,], region = IQR_bd[6,6,])

VP.spbdlg <- VP.spbd %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbdlg$taxa = 'birds'

VP.spbf = data.frame(Climate = IQR_bf[1,1,], habprop = IQR_bf[3,3,], 
                     habConf = IQR_bf[2,2,], Site = IQR_bf[4,4,], 
                     year = IQR_bf[5,5,], region = IQR_bf[6,6,])

VP.spbflg <- VP.spbf %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spbflg$taxa = 'butterflies'


VP.spmoth = data.frame(Climate = IQR_moth[1,1,], habprop = IQR_moth[3,3,], 
                       habConf = IQR_moth[2,2,], Site = IQR_moth[4,4,], 
                       year = IQR_moth[5,5,], region = IQR_moth[6,6,])

VP.spmothlg <- VP.spmoth %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spmothlg$taxa = 'moths'


VP.spwg = data.frame(Climate = IQR_wg[1,1,], habprop = IQR_wg[3,3,], 
                     habConf = IQR_wg[2,2,], Site = IQR_wg[4,4,], 
                     year = IQR_wg[5,5,], region = IQR_wg[6,6,])

VP.spwglg <- VP.spwg %>% 
  pivot_longer(
    cols = 'Climate':'region', 
    names_to = "VP",
    values_to = "value"
  )
VP.spwglg$taxa = 'large mammals'


VPiqr.sp = rbind(VP.sprodlg, VP.spbflg, VP.spbdlg, VP.spwglg, VP.spmothlg)
VPiqr.sp$VP = factor(VPiqr.sp$VP, levels = c('Climate', 'habprop', 'habConf',
                                          'Site', 'year', 'region'))


ggplot(data = VPiqr.sp, aes(x=value, y = taxa, group = taxa, fill = VP)) +
  geom_density_ridges(alpha=0.8) +
  scale_fill_manual(values = c("#66C2A5", 'orange3', "yellow3", "grey78", "blue4", 'green4')) +
  facet_wrap(.~ VP) + xlim(c(-10, 130)) +
  theme_classic() 


## both
##### mcmc

## both

head(VP.mcmc)

VPiqr.mcmc = VPiqr.sp
head(VPiqr.mcmc)

VP.mcmc$Measure = 'Mean'
VPiqr.mcmc$Measure = 'IQR'

VPmeasures.mcmc = rbind(VP.mcmc, VPiqr.mcmc)
VPmeasures.mcmc$Measure = factor(VPmeasures.mcmc$Measure, levels = c('Mean', 'IQR'))
VPmeasures.mcmc$taxa = factor(VPmeasures.mcmc$taxa, levels = rev(c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals')))

levels(VPmeasures.mcmc$VP)[c(2:3, 5:6)] = c("Habcomp", "LandConf", 'Year', 'Region')

ggplot(data = VPmeasures.mcmc, aes(x=value, y = taxa, group = taxa, fill = VP)) +
  geom_density_ridges(alpha=0.8, quantile_lines = TRUE, quantile=2) +
  scale_fill_manual(values = c("#66C2A5", 'orange3', "yellow3", "grey78", "blue4", 'green4')) +
  facet_grid(VP~ Measure) + xlim(c(-5, 100)) +
  labs(title = "", x = "Variance partition",  
       y = "Taxa") + theme_classic() +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) #change font size of legend title  

