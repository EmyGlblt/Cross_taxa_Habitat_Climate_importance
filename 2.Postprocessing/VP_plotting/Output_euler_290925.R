#
#             Summary plots to compare taxa
#
#
# --------------------------------------------------------------------

library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)
library(grid)
library(eulerr)
library(gridExtra)
library(ggplot2)


setwd("D:/Helsinki/RECcII/CSC_HMSC/VPres")

#....................................................................................................
# simplified version .....................................................................................................


#
#  Occurrences
#
## .............................

## Posterior Mean --------------------------------------------------------

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

Mean.fix.ROD = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)

dim(Mean.fix.ROD) 

Mean.fix.rod = round(abs(Mean.fix.ROD)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.rod[1,1], "Habitat composition" = Mean.fix.rod[3,3], "Landscape configuration" = Mean.fix.rod[2,2], 
                    "Climate&Habitat composition" = Mean.fix.rod[1,3], "Habitat composition&Landscape configuration" = Mean.fix.rod[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.rod[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_rod = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                 fill=list(fill=c("#66C2A5", "orange3", 'yellow3', 'gray50', 'gray50', 'gray50')), 
                 counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                 labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_rod


# Change quanity labels manually:
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label =  'climate'
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label =  'Habitat composition'
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.rod[1,1], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.rod[3,3], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.rod[2,2], ' %', sep=' ')

Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste('-', Mean.fix.rod[1,3], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste('-', Mean.fix.rod[2,3], ' %', sep=' ')
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste('-', Mean.fix.rod[2,3], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'


Vplot_rod


####
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


Mean.fix.bf = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.bf) 

Mean.fix.bf = round(abs(Mean.fix.bf)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.bf[1,1], "Habitat composition" = Mean.fix.bf[3,3], "Landscape configuration" = Mean.fix.bf[2,2], 
                    "Climate&Habitat composition" = Mean.fix.bf[1,3], "Habitat composition&Landscape configuration" = Mean.fix.bf[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.bf[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_bf = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                fill=list(fill=c("#66C2A5", "orange3", 'yellow3')), 
                counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_bf


# Change quanity labels manually:

Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.bf[1,1], ' %', sep=' ')
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.bf[3,3], ' %', sep=' ')
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.bf[2,2], ' %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.bf[1,3], ' %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste(Mean.fix.bf[1,2], ' %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.bf[2,3], " [", low.fix.bf[2,3],";", Up.fix.bf[2,3], '] %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_bf



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


Mean.fix.wg = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.wg) 

Mean.fix.wg = round(abs(Mean.fix.wg)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.wg[1,1], "Habitat composition" = Mean.fix.wg[3,3], "Landscape configuration" = Mean.fix.wg[2,2], 
                    "Climate&Habitat composition" = Mean.fix.wg[1,3], "Habitat composition&Landscape configuration" = Mean.fix.wg[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.wg[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_wg = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                fill=list(fill=c("#66C2A5", "orange3", 'yellow3')), 
                counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_wg



# Change quanity labels manually:

Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.wg[1,1], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.wg[3,3], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.wg[2,2], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.wg[1,3], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste(Mean.fix.wg[1,2], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.wg[2,3], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_wg



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
dim(varr)


Mean.fix.bd = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.bd) 

Mean.fix.bd = round(abs(Mean.fix.bd)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.bd[1,1], "Habitat composition" = Mean.fix.bd[3,3], "Landscape configuration" = Mean.fix.bd[2,2], 
                    "Climate&Habitat composition" = Mean.fix.bd[1,3], "Habitat composition&Landscape configuration" = Mean.fix.bd[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.bd[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))

mix_colors <- function(rcol_in) {
  rgb_in <- t(grDevices::col2rgb(rcol_in))
  lab_in <- grDevices::convertColor(rgb_in, "sRGB", "Lab", scale.in = 255)
  mean_col <- colMeans(lab_in)
  rgb_out <- grDevices::convertColor(mean_col, "Lab", "sRGB", scale.out = 1)
  grDevices::rgb(rgb_out)
}

mix_colors(c("yellow3", "orange3"))


Vplot_bird = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                fill=list(fill=c("#66C2A5", "orange3", 'yellow3', '#ACA560', 'gray50', '#CFA900')), 
                counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_bird


# Change quanity labels manually:

Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.bd[1,1], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.bd[3,3], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.bd[2,2], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.bd[1,3], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste(Mean.fix.bd[2,3], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste('-', Mean.fix.bd[1,2], ' %', sep=' ')
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_bird





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


Mean.fix.moth = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.moth) 

Mean.fix.moth = round(abs(Mean.fix.moth)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.moth[1,1], "Habitat composition" = Mean.fix.moth[3,3], "Landscape configuration" = Mean.fix.moth[2,2], 
                    "Climate&Habitat composition" = Mean.fix.moth[1,3], "Habitat composition&Landscape configuration" = Mean.fix.moth[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.moth[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))

mix_colors <- function(rcol_in) {
  rgb_in <- t(grDevices::col2rgb(rcol_in))
  lab_in <- grDevices::convertColor(rgb_in, "sRGB", "Lab", scale.in = 255)
  mean_col <- colMeans(lab_in)
  rgb_out <- grDevices::convertColor(mean_col, "Lab", "sRGB", scale.out = 1)
  grDevices::rgb(rgb_out)
}

mix_colors(c("#66C2A5", "yellow3"))

Vplot_moth = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                fill=list(fill=c("#66C2A5", "orange3", 'yellow3', '#ACA560', '#A6C769', 'gray50')), 
                counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_moth


# Change quanity labels manually:

Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.moth[1,1], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.moth[3,3], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.moth[2,2], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.moth[1,3], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste('-', Mean.fix.moth[2,3], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.moth[1,2], ' %', sep=' ')
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_moth

## combine all


grid.arrange(Vplot_rod, Vplot_bf, Vplot_wg, Vplot_bird, Vplot_moth, ncol = 3, nrow = 2)

##
grid.arrange(Vplot_rod, Vplot_wg, 
             Vplot_bf, Vplot_moth, 
             Vplot_bird, ncol = 2, nrow = 3)

## what is explained per species
Srod = sum(Mean.fix.rod[, -4], na.rm = T) # 7.7
Swg = sum(Mean.fix.wg[, -4], na.rm = T) # 13.1
Smoth = sum(Mean.fix.moth[, -4], na.rm = T) # 8.8
Sbf = sum(Mean.fix.bf[, -4], na.rm = T) # 9.9
Sbird = sum(Mean.fix.bd[, -4], na.rm = T) # 18.4

# Smax = 20 #max(Sbird, Sbf, Smoth, Swg, Srod)
# ## max from the 75% CI
# 
# ### change sizes
# Vplot_rod$vp$height = Vplot_rod$vp$height*(Srod/Smax)
# Vplot_rod$vp$width = Vplot_rod$vp$width*(Srod/Smax)
# Vplot_rod
# 
# Vplot_bf$vp$height = Vplot_bf$vp$height*(Sbf/Smax)
# Vplot_bf$vp$width = Vplot_bf$vp$width*(Sbf/Smax)
# Vplot_bf
# 
# Vplot_moth$vp$height = Vplot_moth$vp$height*(Smoth/Smax)
# Vplot_moth$vp$width = Vplot_moth$vp$width*(Smoth/Smax)
# Vplot_moth
# 
# Vplot_wg$vp$height = Vplot_wg$vp$height*(Swg/Smax)
# Vplot_wg$vp$width = Vplot_wg$vp$width*(Swg/Smax)
# Vplot_wg
# 
# Vplot_bird$vp$height = Vplot_bird$vp$height*(Sbird/Smax)
# Vplot_bird$vp$width = Vplot_bird$vp$width*(Sbird/Smax)
# Vplot_bird
# 
# Empty = ggplot() + theme_void()
# 
# grid.arrange(Vplot_rod, Vplot_wg, 
#              Vplot_moth, Vplot_bf, 
#              Vplot_bird, ncol = 2, nrow = 3)


#
#  Abundances
#
## .............................



#setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")

#....................................................................................................
# simplified version .....................................................................................................


## Mean --------------------------------------------------------

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

Mean.fix.ROD = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)

dim(Mean.fix.ROD) 

Mean.fix.rod = round(abs(Mean.fix.ROD)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.rod[1,1], "Habitat composition" = Mean.fix.rod[3,3], "Landscape configuration" = Mean.fix.rod[2,2], 
                    "Climate&Habitat composition" = Mean.fix.rod[1,3], "Habitat composition&Landscape configuration" = Mean.fix.rod[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.rod[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_rod = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                 fill=list(fill=c("#66C2A5", "orange3", 'yellow3', 'gray50', '#A6C769', 'gray50')), 
                 counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                 labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_rod


# Change quanity labels manually:
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label =  'climate'
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label =  'Habitat composition'
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.rod[1,1], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.rod[3,3], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.rod[2,2], ' %', sep=' ')

Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste('-', Mean.fix.rod[1,3], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste('-', Mean.fix.rod[2,3], ' %', sep=' ')
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste('-', Mean.fix.rod[2,3], ' %', sep=' ')
Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_rod$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'


Vplot_rod


####
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


Mean.fix.bf = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.bf) 

Mean.fix.bf = round(abs(Mean.fix.bf)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.bf[1,1], "Habitat composition" = Mean.fix.bf[3,3], "Landscape configuration" = Mean.fix.bf[2,2], 
                    "Climate&Habitat composition" = Mean.fix.bf[1,3], "Habitat composition&Landscape configuration" = Mean.fix.bf[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.bf[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_bf = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                fill=list(fill=c("#66C2A5", "orange3", 'yellow3', '#ACA560', '#A6C769', 'gray50')), 
                counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_bf


# Change quanity labels manually:

Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.bf[1,1], ' %', sep=' ')
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.bf[3,3], ' %', sep=' ')
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.bf[2,2], ' %', sep=' ')
Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste('-', Mean.fix.bf[2,3], ' %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste(Mean.fix.bf[1,2], ' %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.bf[2,3], " [", low.fix.bf[2,3],";", Up.fix.bf[2,3], '] %', sep=' ')
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_bf$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_bf



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


Mean.fix.wg = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.wg) 

Mean.fix.wg = round(abs(Mean.fix.wg)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.wg[1,1], "Habitat composition" = Mean.fix.wg[3,3], "Landscape configuration" = Mean.fix.wg[2,2], 
                    "Climate&Habitat composition" = Mean.fix.wg[1,3], "Habitat composition&Landscape configuration" = Mean.fix.wg[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.wg[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_wg = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                fill=list(fill=c("#66C2A5", "orange3", 'yellow3', '#ACA560', '#A6C769', 'gray50')), 
                counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_wg



# Change quanity labels manually:

Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.wg[1,1], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.wg[3,3], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.wg[2,2], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.wg[1,3], ' %', sep=' ')
Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste(Mean.fix.wg[1,2], ' %', sep=' ')
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.wg[2,3], " [", low.fix.wg[2,3],";", Up.fix.wg[2,3], '] %', sep=' ')
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_wg$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_wg



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
dim(varr)


Mean.fix.bd = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.bd) 

Mean.fix.bd = round(abs(Mean.fix.bd)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.bd[1,1], "Habitat composition" = Mean.fix.bd[3,3], "Landscape configuration" = Mean.fix.bd[2,2], 
                    "Climate&Habitat composition" = Mean.fix.bd[1,3], "Habitat composition&Landscape configuration" = Mean.fix.bd[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.bd[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))


Vplot_bird = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                  fill=list(fill=c("#66C2A5", "orange3", 'yellow3')), 
                  counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                  labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_bird


# Change quanity labels manually:

Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.bd[1,1], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.bd[3,3], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.bd[2,2], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.bd[1,2], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste(Mean.fix.bd[1,3], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.bd[2,3], ' %', sep=' ')
Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_bird$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_bird





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


Mean.fix.moth = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2), mean, na.rm=T)
dim(Mean.fix.moth) 

Mean.fix.moth = round(abs(Mean.fix.moth)*100,1)


VennDiag <- euler(c("Climate" = Mean.fix.moth[1,1], "Habitat composition" = Mean.fix.moth[3,3], "Landscape configuration" = Mean.fix.moth[2,2], 
                    "Climate&Habitat composition" = Mean.fix.moth[1,3], "Habitat composition&Landscape configuration" = Mean.fix.moth[2,3], 
                    "Climate&Landscape configuration" = Mean.fix.moth[1,2], "Climate&Habitat composition&Landscape configuration" = 0),
                  control = list(extraopt = FALSE))

mix_colors <- function(rcol_in) {
  rgb_in <- t(grDevices::col2rgb(rcol_in))
  lab_in <- grDevices::convertColor(rgb_in, "sRGB", "Lab", scale.in = 255)
  mean_col <- colMeans(lab_in)
  rgb_out <- grDevices::convertColor(mean_col, "Lab", "sRGB", scale.out = 1)
  grDevices::rgb(rgb_out)
}

mix_colors(c("#66C2A5", "yellow3"))

Vplot_moth = plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
                  fill=list(fill=c("#66C2A5", "orange3", 'yellow3', '#ACA560', '#A6C769', 'gray50')), 
                  counts = list(cex=3), shape ='ellipse', input = "union", key = TRUE,#factor_names = TRUE,
                  labels = list(font = 4, cex=1), quantities = list(fontsize = 12), main = '')
Vplot_moth


# Change quanity labels manually:

Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.label.1$label = 'climate'
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.label.2$label = 'Habitat composition'
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.label.3$label = 'Landscape configuration'

Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label = paste(Mean.fix.moth[1,1], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label = paste(Mean.fix.moth[3,3], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label = paste(Mean.fix.moth[2,2], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.4$children$tag.quantity.4$label = paste(Mean.fix.moth[1,3], ' %', sep=' ')
Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.5$children$tag.quantity.5$label = paste('-', Mean.fix.moth[2,3], ' %', sep=' ')
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.6$children$tag.quantity.6$label = paste(Mean.fix.moth[1,2], ' %', sep=' ')
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.7$children$tag.quantity.7$label = ''

#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$vjust = 2
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$vjust = 2
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$vjust = 2

#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$fills.grob.4$gp$fill
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$fills.grob.5$gp$fill = 'white'
#Vplot_moth$children$canvas.grob$children$diagram.grob.1$children$fills.grob.6$gp$fill = 'white'

Vplot_moth

## combine all


##
grid.arrange(Vplot_rod, Vplot_wg, 
             Vplot_bf, Vplot_moth, 
             Vplot_bird, ncol = 2, nrow = 3)
