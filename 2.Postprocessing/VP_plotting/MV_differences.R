
library(dplyr)
library(tidyr)
library(abind)
#library(quantmod)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggridges)
library(tidybayes)


## load all VCNMP variance partitions
setwd("~/Data/4.VPanalyses/VPCP") # to adapt depending on where it is saved
load('VPCP_bd.RDATA')
load('VPCP_bf.RDATA')
load('VPCP_rod.RDATA')
load('VPCP_wg.RDATA')
load('VPCP_moth.RDATA')

## load all linear predictor correlation
#load('FRcorr_taxa1010.RDATA')

# example dimensions
dim(VP.test3_bf$Vnorm) #    1  /   predictors group  /  species /  mcmc samples
#dim(FRcor_bf) #    predictors group (+ effort in our case)  /  predictors group (+ effort in our case)  /  species /  mcmc samples
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')      # + 3 random effect so 7 in total


# Scaling by R2 all VP summaries
# 
# rod
varr <- aperm(array(TabMF_rod$TjurR2, dim = c(8L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_rod$Vnorm = varr * VP.test3_rod$Vnorm
VP.test3_rod$Vdiag = varr * VP.test3_rod$Vdiag


## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_rod$TjurR2, dim = c(8L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_rod$Cnorm = Carr * VP.test3_rod$Cnorm
VP.test3_rod$Vmarg = Carr * VP.test3_rod$Vmarg
VP.test3_rod$Vpart = Carr * VP.test3_rod$Vpart

# bf
varr <- aperm(array(TabMF_bf$TjurR2, dim = c(57L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_bf$Vnorm = varr * VP.test3_bf$Vnorm
VP.test3_bf$Vdiag = varr * VP.test3_bf$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_bf$TjurR2, dim = c(57L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_bf$Cnorm = Carr * VP.test3_bf$Cnorm
VP.test3_bf$Vmarg = Carr * VP.test3_bf$Vmarg
VP.test3_bf$Vpart = Carr * VP.test3_bf$Vpart

# moth
varr <- aperm(array(TabMF_Moth$TjurR2, dim = c(319L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_moth$Vnorm = varr * VP.test3_moth$Vnorm
VP.test3_moth$Vdiag = varr * VP.test3_moth$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_Moth$TjurR2, dim = c(319L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_moth$Cnorm = Carr * VP.test3_moth$Cnorm
VP.test3_moth$Vmarg = Carr * VP.test3_moth$Vmarg
VP.test3_moth$Vpart = Carr * VP.test3_moth$Vpart

# wg
varr <- aperm(array(TabMF_wg$TjurR2, dim = c(15L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_wg$Vnorm = varr * VP.test3_wg$Vnorm
VP.test3_wg$Vdiag = varr * VP.test3_wg$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_wg$TjurR2, dim = c(15L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_wg$Cnorm = Carr * VP.test3_wg$Cnorm
VP.test3_wg$Vmarg = Carr * VP.test3_wg$Vmarg
VP.test3_wg$Vpart = Carr * VP.test3_wg$Vpart


# bird
varr <- aperm(array(TabMF_bd$TjurR2, dim = c(102L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_bd$Vnorm = varr * VP.test3_bd$Vnorm
VP.test3_bd$Vdiag = varr * VP.test3_bd$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_bd$TjurR2, dim = c(102L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_bd$Cnorm = Carr * VP.test3_bd$Cnorm
VP.test3_bd$Vmarg = Carr * VP.test3_bd$Vmarg
VP.test3_bd$Vpart = Carr * VP.test3_bd$Vpart

## Posterior differences between M and V
VP.test3_bd$Vmarg[,,1,1]
VP.test3_bd$Vnorm[,,1,1]

mat.diff.bd = array(data=NA, dim = dim(VP.test3_bd$Vmarg))
for (i in 1:dim(VP.test3_bd$Vmarg)[4]){
  for (j in 1:dim(VP.test3_bd$Vmarg)[3]) {
    
    mat.diff.bd[,,j,i] = sweep(VP.test3_bd$Vmarg[,,j,i], 2, VP.test3_bd$Vnorm[,,j,i], FUN = "-")

  }
}
dim(mat.diff.bd)
dimnames(mat.diff.bd) = dimnames(VP.test3_bd$Vmarg)


mat.diff.bf = array(data=NA, dim = dim(VP.test3_bf$Vmarg))
for (i in 1:dim(VP.test3_bf$Vmarg)[4]){
  for (j in 1:dim(VP.test3_bf$Vmarg)[3]) {
    
    mat.diff.bf[,,j,i] = sweep(VP.test3_bf$Vmarg[,,j,i], 2, VP.test3_bf$Vnorm[,,j,i], FUN = "-")
    
  }
}
dimnames(mat.diff.bf) = dimnames(VP.test3_bf$Vmarg)



mat.diff.rod = array(data=NA, dim = dim(VP.test3_rod$Vmarg))
for (i in 1:dim(VP.test3_rod$Vmarg)[4]){
  for (j in 1:dim(VP.test3_rod$Vmarg)[3]) {
    
    mat.diff.rod[,,j,i] = sweep(VP.test3_rod$Vmarg[,,j,i], 2, VP.test3_rod$Vnorm[,,j,i], FUN = "-")
    
  }
}
dimnames(mat.diff.rod) = dimnames(VP.test3_rod$Vmarg)


mat.diff.moth = array(data=NA, dim = dim(VP.test3_moth$Vmarg))
for (i in 1:dim(VP.test3_moth$Vmarg)[4]){
  for (j in 1:dim(VP.test3_moth$Vmarg)[3]) {
    
    mat.diff.moth[,,j,i] = sweep(VP.test3_moth$Vmarg[,,j,i], 2, VP.test3_moth$Vnorm[,,j,i], FUN = "-")
    
  }
}
dimnames(mat.diff.moth) = dimnames(VP.test3_moth$Vmarg)


mat.diff.wg = array(data=NA, dim = dim(VP.test3_wg$Vmarg))
for (i in 1:dim(VP.test3_wg$Vmarg)[4]){
  for (j in 1:dim(VP.test3_wg$Vmarg)[3]) {
    
    mat.diff.wg[,,j,i] = sweep(VP.test3_wg$Vmarg[,,j,i], 2, VP.test3_wg$Vnorm[,,j,i], FUN = "-")
    
  }
}
dimnames(mat.diff.wg) = dimnames(VP.test3_wg$Vmarg)


## posterior mean

DiffMV_bd = melt(apply(mat.diff.bd, c(1,2,3), mean)) 
DiffMV_bd$Taxa = 'birds'

DiffMV_bf = melt(apply(mat.diff.bf, c(1,2,3), mean)) 
DiffMV_bf$Taxa = 'butterflies'

DiffMV_moth = melt(apply(mat.diff.moth, c(1,2,3), mean)) 
DiffMV_moth$Taxa = 'moths'

DiffMV_rod = melt(apply(mat.diff.rod, c(1,2,3), mean)) 
DiffMV_rod$Taxa = 'small mammals'

DiffMV_wg = melt(apply(mat.diff.wg, c(1,2,3), mean)) 
DiffMV_wg$Taxa = 'large mammals'


## combine all

DiffMV = rbind(DiffMV_bd, DiffMV_bf, DiffMV_moth, DiffMV_rod, DiffMV_wg)
DiffMV = DiffMV[!(DiffMV$Value == DiffMV$element) ,]
DiffMV = DiffMV[!(DiffMV$Value %in% c('site', 'year', 'bg') | DiffMV$element %in% c('site', 'year', 'bg')) ,]

DiffMV$Taxa = factor(DiffMV$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
DiffMV$Value = factor(DiffMV$Value, levels = c("Climate", "HabitatConf","Habitatprop"),
                      labels = c("Climate", "LandscapeConf","Habitatcomp"))
DiffMV$element = factor(DiffMV$element, levels = c("Climate", "HabitatConf","Habitatprop"),
                        labels = c("Climate", "LandscapeConf","Habitatcomp"))


## plotting

ggplot(data = DiffMV, aes(y = element, x=value, fill=Value)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75,
                      quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = c("#66C2A5", "yellow3", 'orange3')) +
  facet_grid(Taxa ~ Value) +
  theme_classic() 


ggplot(data = DiffMV, aes(y = element, x=value, fill=Value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = T) +
  geom_vline(xintercept = 0, colour='red', linetype = "dashed") +
  scale_fill_manual(values = c("#66C2A5", "yellow3", 'orange3')) +
  facet_grid(Taxa ~ Value) +
  theme_classic() 


plot = ggplot(data = DiffMV, aes(y = element, x=value, fill=element)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = T) +
  geom_vline(xintercept = 0, colour='red', linetype = "dashed") +
  scale_fill_manual(values = c("#66C2A5", "yellow3", 'orange3')) + 
  coord_flip() + ylab('') + xlab('M-V') +
  facet_grid(Value ~ Taxa) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

g <- ggplot_gtable(ggplot_build(plot))
strips <- which(grepl('strip-', g$layout$name))

col_strip = c(rep("white",5), "#66C2A5", "yellow3", 'orange3')
col_text = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", rep("black",3))

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_strip[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$col <- col_text[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- col_text[i]
  
}

plot(g)
