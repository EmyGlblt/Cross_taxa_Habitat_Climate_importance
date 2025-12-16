#
#         Ternplots plots : Global
#
#
# ..............................................................................................

library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)
library(ggplot2)
library(ggExtra)
library(cowplot)
library(patchwork)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggridges)

library(ggtern)


# colors
library(wesanderson)
names(wes_palettes)

#remotes::install_github("wilkelab/cowplot")
#install.packages("colorspace", repos = "http://R-Forge.R-project.org")
#remotes::install_github("clauswilke/colorblindr")
library(colorblindr)

setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")

#
#   Occurrence
#

## Load species specific variance partition summarises: focusing on the normalized VP
##..........................................................................

load('VPCP_rod.RDATA')
dim(VP.test3_rod$Vdiag)
dim(VP.test3_rod$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_rod$Cnorm))

for (i in 1:dim(VP.test3_rod$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_rod$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_rod$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_rod$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_rod$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_rod$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_rod$TjurR2, dim = c(1000L, 8L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)

Mean.fix.sprod = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_bf.RDATA')
dim(VP.test3_bf$Vdiag)
dim(VP.test3_bf$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bf$Cnorm))

for (i in 1:dim(VP.test3_bf$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bf$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bf$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bf$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bf$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bf$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_bf$TjurR2, dim = c(1000L, 57L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)


Mean.fix.spbf = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_wg.RDATA')
dim(VP.test3_wg$Vdiag)
dim(VP.test3_wg$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_wg$Cnorm))

for (i in 1:dim(VP.test3_wg$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_wg$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_wg$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_wg$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_wg$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_wg$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_wg$TjurR2, dim = c(1000L, 15L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)


Mean.fix.spwg = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_bd.RDATA')
dim(VP.test3_bd$Vdiag)
dim(VP.test3_bd$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bd$Cnorm))

for (i in 1:dim(VP.test3_bd$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bd$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bd$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bd$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bd$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bd$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_bd$TjurR2, dim = c(1000L, 102L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))

Mean.fix.spbd = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_moth.RDATA')
dim(VP.test3_moth$Vdiag)
dim(VP.test3_moth$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_moth$Cnorm))

for (i in 1:dim(VP.test3_moth$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_moth$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_moth$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_moth$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_moth$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_moth$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_Moth$TjurR2, dim = c(1000L, 319L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)


Mean.fix.spmoth = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)




## group information

Mean.fix.sprod = Mean.fix.sprod*100
Mean.fix.spbf = Mean.fix.spbf*100
Mean.fix.spbd = Mean.fix.spbd*100
Mean.fix.spmoth = Mean.fix.spmoth*100
Mean.fix.spwg = Mean.fix.spwg*100


VP.sprod = data.frame(Climate = Mean.fix.sprod[1,1,], habprop = Mean.fix.sprod[3,3,], 
                      habConf = Mean.fix.sprod[2,2,])

VP.sprod$species = TabMF_rod$species
VP.sprod$taxa = 'rodents'

VP.spmoth = data.frame(Climate = Mean.fix.spmoth[1,1,], habprop = Mean.fix.spmoth[3,3,], 
                       habConf = Mean.fix.spmoth[2,2,])

VP.spmoth$species = TabMF_Moth$species
VP.spmoth$taxa = 'moths'

VP.spbd = data.frame(Climate = Mean.fix.spbd[1,1,], habprop = Mean.fix.spbd[3,3,], 
                     habConf = Mean.fix.spbd[2,2,])

VP.spbd$species = TabMF_bd$species
VP.spbd$taxa = 'birds'

VP.spbf = data.frame(Climate = Mean.fix.spbf[1,1,], habprop = Mean.fix.spbf[3,3,], 
                     habConf = Mean.fix.spbf[2,2,])

VP.spbf$species = TabMF_bf$species
VP.spbf$taxa = 'butterflies'

VP.spwg = data.frame(Climate = Mean.fix.spwg[1,1,], habprop = Mean.fix.spwg[3,3,], 
                     habConf = Mean.fix.spwg[2,2,])

VP.spwg$species = TabMF_wg$species
VP.spwg$taxa = 'large mammals'


VP.sp = rbind(VP.sprod, VP.spbf, VP.spmoth, VP.spbd, VP.spwg)
VP.sp$taxa = as.factor(VP.sp$taxa)
VP.sp$taxa = factor(VP.sp$taxa, levels = c('birds', 'butterflies', 'moths',
                                              'rodents', 'large mammals'))
summary(VP.sp)
VP.sp2 = VP.sp[, -4]

#  if rescale already to sum up to 100
## it is the same

VP.sp2$sum = VP.sp2$Climate + VP.sp2$habprop + VP.sp2$habConf

VP.sp2$Climate = VP.sp2$Climate*100/VP.sp2$sum
VP.sp2$habprop = VP.sp2$habprop*100/VP.sp2$sum
VP.sp2$habConf = VP.sp2$habConf*100/VP.sp2$sum

## add group means
VP.spmean = VP.sp2 %>%
  group_by(taxa) %>%
  summarise_all(list(mean))

Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")

Tern_occ = ggtern(VP.sp2, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
                    #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  geom_point(aes(color = taxa), shape = 4, size = 0.7) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmean, lty = 2, size = 2) +  # add mean information and lines
  geom_point(data = VP.spmean, aes(color = taxa), show.legend = F) + 
  #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
  #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
  theme_showarrows() +
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)) + 
  guides(alpha='none', colour='none') + 
  labs(x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave( Tern_occ, file=paste0("VTern_occ",".png", sep = ""))

## checking if this colour setting is colourblind friendly
cvd_grid(Tern_occ)

#ggglobal
leg <- get_legend(Tern_occ)

# Convert to a ggplot and print
legend = as_ggplot(leg)
ggsave( as_ggplot(leg), file=paste0("VPtern_legend",".png", sep = ""))

# Marginal density plot of x (top panel) and y (right panel)
climplot <- ggdensity(VP.sp2, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave(climplot, file=paste0("climplot_occ",".png", sep = ""))

habPplot <- ggdensity(VP.sp2, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave(habPplot, file=paste0("habPplot_occ",".png", sep = ""))

LandCplot <- ggdensity(VP.sp2, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave(LandCplot, file=paste0("LandCplot_occ",".png", sep = ""))



##..........................................................................
##..........................................................................



#
#   Abundance
#

## Load species specific variance partition summarises: focusing on the normalized VP
##..........................................................................

load('VPCP_ab_rod.RDATA')
dim(VP.test3_rod$Vdiag)
dim(VP.test3_rod$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_rod$Cnorm))

for (i in 1:dim(VP.test3_rod$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_rod$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_rod$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_rod$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_rod$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_rod$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_rod$SR2, dim = c(2000L, 8L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)

Mean.fix.sprod = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_ab_bf.RDATA')
dim(VP.test3_bf$Vdiag)
dim(VP.test3_bf$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bf$Cnorm))

for (i in 1:dim(VP.test3_bf$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bf$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bf$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bf$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bf$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bf$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_bf$SR2, dim = c(1000L, 57L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)


Mean.fix.spbf = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_ab_wg.RDATA')
dim(VP.test3_wg$Vdiag)
dim(VP.test3_wg$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_wg$Cnorm))

for (i in 1:dim(VP.test3_wg$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_wg$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_wg$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_wg$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_wg$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_wg$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_wg$SR2, dim = c(1000L, 15L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)


Mean.fix.spwg = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_ab_bd.RDATA')
dim(VP.test3_bd$Vdiag)
dim(VP.test3_bd$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bd$Cnorm))

for (i in 1:dim(VP.test3_bd$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bd$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bd$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bd$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bd$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bd$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_bd$SR2, dim = c(1000L, 102L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))

Mean.fix.spbd = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)



load('VPCP_ab_moth.RDATA')
dim(VP.test3_moth$Vdiag)
dim(VP.test3_moth$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_moth$Cnorm))

for (i in 1:dim(VP.test3_moth$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_moth$Vdiag)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_moth$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_moth$Vdiag[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_moth$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_moth$Cnorm[,,j,i])
    
  }
}


varr <- aperm(array(TabMF_Moth$SR2, dim = c(1000L, 319L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)


Mean.fix.spmoth = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2, 3), mean, na.rm=T)




## group information

Mean.fix.sprod = Mean.fix.sprod*100
Mean.fix.spbf = Mean.fix.spbf*100
Mean.fix.spbd = Mean.fix.spbd*100
Mean.fix.spmoth = Mean.fix.spmoth*100
Mean.fix.spwg = Mean.fix.spwg*100


VP.sprod = data.frame(Climate = Mean.fix.sprod[1,1,], habprop = Mean.fix.sprod[3,3,], 
                      habConf = Mean.fix.sprod[2,2,])

VP.sprod$species = TabMF_rod$species
VP.sprod$taxa = 'rodents'

VP.spmoth = data.frame(Climate = Mean.fix.spmoth[1,1,], habprop = Mean.fix.spmoth[3,3,], 
                       habConf = Mean.fix.spmoth[2,2,])

VP.spmoth$species = TabMF_Moth$species
VP.spmoth$taxa = 'moths'

VP.spbd = data.frame(Climate = Mean.fix.spbd[1,1,], habprop = Mean.fix.spbd[3,3,], 
                     habConf = Mean.fix.spbd[2,2,])

VP.spbd$species = TabMF_bd$species
VP.spbd$taxa = 'birds'

VP.spbf = data.frame(Climate = Mean.fix.spbf[1,1,], habprop = Mean.fix.spbf[3,3,], 
                     habConf = Mean.fix.spbf[2,2,])

VP.spbf$species = TabMF_bf$species
VP.spbf$taxa = 'butterflies'

VP.spwg = data.frame(Climate = Mean.fix.spwg[1,1,], habprop = Mean.fix.spwg[3,3,], 
                     habConf = Mean.fix.spwg[2,2,])

VP.spwg$species = TabMF_wg$species
VP.spwg$taxa = 'large mammals'


VP.sp = rbind(VP.sprod, VP.spbf, VP.spmoth, VP.spbd, VP.spwg)
VP.sp$taxa = as.factor(VP.sp$taxa)
VP.sp$taxa = factor(VP.sp$taxa, levels = c('birds', 'butterflies', 'moths',
                                           'rodents', 'large mammals'))
summary(VP.sp)
VP.sp2 = VP.sp[, -4]

#  if rescale already to sum up to 100
## it is the same

VP.sp2$sum = VP.sp2$Climate + VP.sp2$habprop + VP.sp2$habConf

VP.sp2$Climate = VP.sp2$Climate*100/VP.sp2$sum
VP.sp2$habprop = VP.sp2$habprop*100/VP.sp2$sum
VP.sp2$habConf = VP.sp2$habConf*100/VP.sp2$sum

## add group means
VP.spmean = VP.sp2 %>%
  group_by(taxa) %>%
  summarise_all(list(mean))

Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")

Tern_ab = ggtern(VP.sp2, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  geom_point(aes(color = taxa), shape = 4, size = 0.7) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmean, lty = 2, size = 2) +  # add mean information and lines
  geom_point(data = VP.spmean, aes(color = taxa), show.legend = F) + 
  #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
  #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
  theme_showarrows() +
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)) + 
  guides(alpha='none', colour='none') + 
  labs(x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")


ggsave(Tern_ab, file=paste0("VTern_ab",".png", sep = ""))

## checking if this colour setting is colourblind friendly
cvd_grid(Tern_ab)

# Marginal density plot of x (top panel) and y (right panel)
climplot <- ggdensity(VP.sp2, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave(climplot, file=paste0("climplot_ab",".png", sep = ""))

habPplot <- ggdensity(VP.sp2, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave(habPplot, file=paste0("habPplot_ab",".png", sep = ""))

LandCplot <- ggdensity(VP.sp2, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave(LandCplot, file=paste0("LandCplot_ab",".png", sep = ""))




