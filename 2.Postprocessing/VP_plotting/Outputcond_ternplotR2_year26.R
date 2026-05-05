#
#             Ternplot conditional variance
#
#
# --------------------------------------------------------------------

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

setwd("~/Data/4.VPanalyses/VP_cond_trait") # to adapt depending on where it is saved

#
#   Occurrence
#

## Load species specific variance partition summarises: focusing on the normalized VP
##..........................................................................

load('VPcondwb_temp_Rod.RDATA')
dim(VP.test1_rod$Vdiag)
dim(VP.test2_rod$Vdiag)

## conditional
VP.test_cond = array(NA, dim = dim(VP.test1_rod$Cnorm), dimnames = dimnames(VP.test1_rod$Cnorm))

for (i in 1:dim(VP.test1_rod$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_rod$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_rod$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_rod$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_rod$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_rod$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_rod$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.rod <- aperm(array(TabMF_rod$TjurR2, dim = c(8L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.rod)


Mean_fix_rod = array(NA, dim = c(3,3,dim(VP.test1_rod$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_rod$Vdiag)[5]){
  Mean_fix_rod[,,, k] = apply(varr.rod*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_rod$Cnorm), dimnames = dimnames(VP.test2_rod$Cnorm)) 

for (i in 1:dim(VP.test2_rod$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_rod$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_rod$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_rod$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_rod$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_rod$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_rod$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.rodbtw = apply(varr.rod*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.rodwit = apply(varr.rod*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


load('VPcondwb_temp_bf.RDATA')
dim(VP.test1_bf$Vdiag)
dim(VP.test2_bf$Vdiag)

# conditional
VP.test_cond = array(NA, dim = dim(VP.test1_bf$Cnorm), dimnames = dimnames(VP.test1_bf$Cnorm))

for (i in 1:dim(VP.test1_bf$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_bf$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_bf$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_bf$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_bf$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_bf$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_bf$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.bf <- aperm(array(TabMF_bf$TjurR2, dim = c(57L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.bf)


Mean_fix_bf = array(NA, dim = c(3,3,dim(VP.test1_bf$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_bf$Vdiag)[5]){
  Mean_fix_bf[,,, k] = apply(varr.bf*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_bf$Cnorm), dimnames = dimnames(VP.test2_bf$Cnorm)) 

for (i in 1:dim(VP.test2_bf$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_bf$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_bf$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_bf$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_bf$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_bf$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_bf$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.bfbtw = apply(varr.bf*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.bfwit = apply(varr.bf*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


load('VPcondwb_temp_moth.RDATA')
dim(VP.test1_moth$Vdiag)
dim(VP.test2_moth$Vdiag)

VP.test_cond = array(NA, dim = dim(VP.test1_moth$Cnorm), dimnames = dimnames(VP.test1_moth$Cnorm))

for (i in 1:dim(VP.test1_moth$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_moth$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_moth$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_moth$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_moth$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_moth$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_moth$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.moth <- aperm(array(TabMF_Moth$TjurR2, dim = c(319L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))


Mean_fix_moth = array(NA, dim = c(3,3,dim(VP.test1_moth$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_moth$Vdiag)[5]){
  Mean_fix_moth[,,, k] = apply(varr.moth*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_moth$Cnorm), dimnames = dimnames(VP.test2_moth$Cnorm)) 

for (i in 1:dim(VP.test2_moth$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_moth$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_moth$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_moth$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_moth$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_moth$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_moth$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.mothbtw = apply(varr.moth*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.mothwit = apply(varr.moth*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 



load('VPcondwb_temp_wg.RDATA')
dim(VP.test1_wg$Vdiag)
dim(VP.test2_wg$Vdiag)

# conditional
VP.test_cond = array(NA, dim = dim(VP.test1_wg$Cnorm), dimnames = dimnames(VP.test1_wg$Cnorm))

for (i in 1:dim(VP.test1_wg$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_wg$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_wg$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_wg$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_wg$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_wg$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_wg$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.wg <- aperm(array(TabMF_wg$TjurR2, dim = c(15L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.wg)


Mean_fix_wg = array(NA, dim = c(3,3,dim(VP.test1_wg$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_wg$Vdiag)[5]){
  Mean_fix_wg[,,, k] = apply(varr.wg*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_wg$Cnorm), dimnames = dimnames(VP.test2_wg$Cnorm)) 

for (i in 1:dim(VP.test2_wg$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_wg$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_wg$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_wg$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_wg$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_wg$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_wg$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.wgbtw = apply(varr.wg*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.wgwit = apply(varr.wg*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


load('VPcondwb_temp_bd.RDATA')
dim(VP.test1_bd$Vdiag)
dim(VP.test2_bd$Vdiag)

# conditional
VP.test_cond = array(NA, dim = dim(VP.test1_bd$Cnorm), dimnames = dimnames(VP.test1_bd$Cnorm))

for (i in 1:dim(VP.test1_bd$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_bd$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_bd$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_bd$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_bd$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_bd$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_bd$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.bd <- aperm(array(TabMF_bd$TjurR2, dim = c(102L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.bd)


Mean_fix_bd = array(NA, dim = c(3,3,dim(VP.test1_bd$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_bd$Vdiag)[5]){
  Mean_fix_bd[,,, k] = apply(varr.bd*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}



# btw wit
## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_bd$Cnorm), dimnames = dimnames(VP.test2_bd$Cnorm)) 

for (i in 1:dim(VP.test2_bd$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_bd$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_bd$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_bd$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_bd$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_bd$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_bd$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.bdbtw = apply(varr.bd*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.bdwit = apply(varr.bd*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggpubr)

VP.year = list()

for (k in 1:dim(VP.test1_bd$Vdiag)[5]){ # Same number of conditions year for all taxa
  Mean_fix_bd[,,, k] = Mean_fix_bd[,,, k]*100
  Mean_fix_bf[,,, k] = Mean_fix_bf[,,, k]*100
  Mean_fix_moth[,,, k] = Mean_fix_moth[,,, k]*100
  Mean_fix_rod[,,, k] = Mean_fix_rod[,,, k]*100
  Mean_fix_wg[,,, k] = Mean_fix_wg[,,, k]*100
  
  VP.rod = data.frame(Climate = Mean_fix_rod[1,1,, k], habprop = Mean_fix_rod[3,3,, k], 
                        habConf = Mean_fix_rod[2,2,, k])
  VP.rod$taxa = 'rodent'
  VP.rod$BG = dimnames(VP.test1_rod$Vdiag)$group[k]
  
  VP.bd = data.frame(Climate = Mean_fix_bd[1,1,, k], habprop = Mean_fix_bd[3,3,, k], 
                       habConf = Mean_fix_bd[2,2,, k])
  VP.bd$taxa = 'bird'
  VP.bd$BG = dimnames(VP.test1_bd$Vdiag)$group[k]
  
  VP.bf = data.frame(Climate = Mean_fix_bf[1,1,, k], habprop = Mean_fix_bf[3,3,, k], 
                       habConf = Mean_fix_bf[2,2,, k])
  VP.bf$taxa = 'butterfly'
  VP.bf$BG = dimnames(VP.test1_bf$Vdiag)$group[k]
  
  VP.moth = data.frame(Climate = Mean_fix_moth[1,1,, k], habprop = Mean_fix_moth[3,3,, k], 
                         habConf = Mean_fix_moth[2,2,, k])
  VP.moth$taxa = 'moth'
  VP.moth$BG = dimnames(VP.test1_moth$Vdiag)$group[k]
  
  
  VP.wg = data.frame(Climate = Mean_fix_wg[1,1,, k], habprop = Mean_fix_wg[3,3,, k], 
                       habConf = Mean_fix_wg[2,2,, k])
  VP.wg$taxa = 'large mammals'
  VP.wg$BG = dimnames(VP.test1_wg$Vdiag)$group[k]
  
  
  VP.year[[k]] = rbind(VP.rod, VP.bf, VP.moth, VP.bd, VP.wg)
}




Mean.fix.rodbtw = Mean.fix.rodbtw*100
Mean.fix.bfbtw = Mean.fix.bfbtw*100
Mean.fix.mothbtw = Mean.fix.mothbtw*100
Mean.fix.bdbtw = Mean.fix.bdbtw*100
Mean.fix.wgbtw = Mean.fix.wgbtw*100


VP.rodbtw = data.frame(Climate = Mean.fix.rodbtw[1,1,], habprop = Mean.fix.rodbtw[3,3,], 
                         habConf = Mean.fix.rodbtw[2,2,])

VP.rodbtw$taxa = 'rodent'

VP.bdbtw = data.frame(Climate = Mean.fix.bdbtw[1,1,], habprop = Mean.fix.bdbtw[3,3,], 
                        habConf = Mean.fix.bdbtw[2,2,])

VP.bdbtw$taxa = 'bird'

VP.bfbtw = data.frame(Climate = Mean.fix.bfbtw[1,1,], habprop = Mean.fix.bfbtw[3,3,], 
                        habConf = Mean.fix.bfbtw[2,2,])

VP.bfbtw$taxa = 'butterfly'

VP.mothbtw = data.frame(Climate = Mean.fix.mothbtw[1,1,], habprop = Mean.fix.mothbtw[3,3,], 
                          habConf = Mean.fix.mothbtw[2,2,])

VP.mothbtw$taxa = 'moth'

VP.wgbtw = data.frame(Climate = Mean.fix.wgbtw[1,1,], habprop = Mean.fix.wgbtw[3,3,], 
                        habConf = Mean.fix.wgbtw[2,2,])

VP.wgbtw$taxa = 'large mammals'


VP.btw = rbind(VP.rodbtw, VP.bfbtw, VP.mothbtw, VP.bdbtw, VP.wgbtw)
VP.btw$BG = 'btw'



Mean.fix.rodwit = Mean.fix.rodwit*100
Mean.fix.bfwit = Mean.fix.bfwit*100
Mean.fix.mothwit = Mean.fix.mothwit*100
Mean.fix.bdwit = Mean.fix.bdwit*100
Mean.fix.wgwit = Mean.fix.wgwit*100

VP.rodwit = data.frame(Climate = Mean.fix.rodwit[1,1,], habprop = Mean.fix.rodwit[3,3,], 
                         habConf = Mean.fix.rodwit[2,2,])

VP.rodwit$taxa = 'rodent'

VP.bdwit = data.frame(Climate = Mean.fix.bdwit[1,1,], habprop = Mean.fix.bdwit[3,3,], 
                        habConf = Mean.fix.bdwit[2,2,])

VP.bdwit$taxa = 'bird'

VP.bfwit = data.frame(Climate = Mean.fix.bfwit[1,1,], habprop = Mean.fix.bfwit[3,3,], 
                        habConf = Mean.fix.bfwit[2,2,])

VP.bfwit$taxa = 'butterfly'


VP.mothwit = data.frame(Climate = Mean.fix.mothwit[1,1,], habprop = Mean.fix.mothwit[3,3,], 
                          habConf = Mean.fix.mothwit[2,2,])

VP.mothwit$taxa = 'moth'

VP.wgwit = data.frame(Climate = Mean.fix.wgwit[1,1,], habprop = Mean.fix.wgwit[3,3,], 
                        habConf = Mean.fix.wgwit[2,2,])

VP.wgwit$taxa = 'large mammals'


VP.wit = rbind(VP.rodwit, VP.bfwit, VP.mothwit, VP.bdwit, VP.wgwit)
VP.wit$BG = 'wit'


VP.year.array = do.call(rbind, VP.year)

VP.sp = rbind(VP.btw, VP.wit, VP.year.array)
VP.sp$BG = factor(VP.sp$BG, levels = c('1999', '2000', '2001', '2002', '2003', '2004', '2005','2006',
                                       '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', 
                                       '2015', '2016', '2017', '2018', '2019','btw', 'wit'))
VP.sp$taxa = factor(VP.sp$taxa, levels = c('bird', 'butterfly', 'moth', 'rodent', 'large mammals'))


# colors
library(wesanderson)
names(wes_palettes)


library(grid)
head(VP.sp)
summary(VP.sp)

# 
# VP.sp2 = VP.sp[VP.sp$Climate<100 & VP.sp$Climate>0,]
# VP.sp3 = VP.sp2[VP.sp2$habprop<100 & VP.sp2$habprop>0,]
# VP.sp = VP.sp3[VP.sp3$habConf<100 & VP.sp3$habConf>0,]
# 

VP.sp$sum = VP.sp$Climate + VP.sp$habprop + VP.sp$habConf

VP.sp$Climate = VP.sp$Climate*100/VP.sp$sum
VP.sp$habprop = VP.sp$habprop*100/VP.sp$sum
VP.sp$habConf = VP.sp$habConf*100/VP.sp$sum

summary(VP.sp)


Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")


VP.sp.btw = VP.sp[VP.sp$BG == 'btw',]

VP.spmeanbtw = VP.sp.btw[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanbtw)[2:4] = c('Climate', 'habprop', 'habConf')


ternbtw = ggtern(VP.sp.btw, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanbtw, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanbtw, aes(color = taxa), show.legend = F) + 
  #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
  #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
  #theme_showarrows() +
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 28)) + 
  guides(alpha='none', colour='none') + 
  labs(title = 'Between group', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternbtw, file=paste0("VTernBTW_occ",".png", sep = ""))

climplot <- ggdensity(VP.sp.btw, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_btw.png")

habPplot <- ggdensity(VP.sp.btw, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_btw.png")

habCplot <- ggdensity(VP.sp.btw, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_btw.png")



VP.sp.wit = VP.sp[VP.sp$BG == 'wit',]

VP.spmeanwit = VP.sp.wit[, c(1:6)] %>%
  group_by(taxa) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanwit)[2:4] = c('Climate', 'habprop', 'habConf')

ternwit = ggtern(VP.sp.wit, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", 
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanwit, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanwit, aes(color = taxa), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 28)) + 
  guides(alpha='none', colour='none') + 
  labs(title = 'Within group', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternwit, file=paste0("VTernWIT_occ",".png", sep = ""))


climplot <- ggdensity(VP.sp.wit, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_wit.png")

habPplot <- ggdensity(VP.sp.wit, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_wit.png")

habCplot <- ggdensity(VP.sp.wit, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_wit.png")


# or do iteratively
# 

year_cond = levels(VP.sp$BG)

for (i in 1:length(year_cond)) {
  VP.sp.cond = VP.sp[VP.sp$BG == year_cond[i],]
  
  VP.spmean = VP.sp.cond[, c(1:4)] %>%
    group_by(taxa) %>%
    #summarise(across(everything(), list(mean)))
    summarise(across(everything(), list(mean)))
  colnames(VP.spmean)[2:4] = c('Climate', 'habprop', 'habConf')
  
  
  
  tern = ggtern(VP.sp.cond, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
    #geom_point()+
    stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                      position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
    #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
    scale_alpha_continuous(range = c(0.1, 0.3)) +
    scale_fill_manual(values = Newcol) +
    #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
    scale_color_manual(values = Newcol) +
    geom_crosshair_tern(data = VP.spmean, lty = 2, linewidth = 3.2) +  # add mean information and lines
    geom_point(data = VP.spmean, aes(color = taxa), show.legend = F) + 
    #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
    #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
    #theme_showarrows() +
    theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
          tern.axis.arrow.T = element_line(size=0, color="white"),
          tern.axis.arrow.R = element_line(size=0, color="white"),
          tern.axis.line.L = element_line(color='orange3',size=2),
          tern.axis.line.T = element_line(color='#66C2A5',size=2),
          tern.axis.line.R = element_line(color='yellow3',size=2),
          axis.text = element_text(size = 28),
          axis.title = element_text(size = 28)) + 
    guides(alpha='none', colour='none') + 
    labs(title = year_cond[i], x='', y='', z='', 
         xarrow = "Habitat proportions",
         yarrow = "Climate",
         zarrow = "Habitat configuration")
  ggsave(tern, file=paste0("VTern_occ", year_cond[i],".png", sep = ""))
  
  climplot <- ggdensity(VP.sp.cond, 'Climate', fill = 'taxa') +
    scale_fill_manual(values = Newcol) +
    scale_y_reverse()+ scale_x_reverse() +
    xlim(100,0) + guides(fill='none') + clean_theme() 
  ggsave(paste0("climglobalocc", year_cond[i],".png", sep = ""))
  
  habPplot <- ggdensity(VP.sp.cond, "habprop",  fill = 'taxa') +
    scale_fill_manual(values = Newcol) +
    xlim(0,100) + guides(fill='none') + clean_theme() 
  ggsave(paste0("habPglobalocc", year_cond[i],".png", sep = ""))
  
  habCplot <- ggdensity(VP.sp.cond, "habConf",  fill = 'taxa') +
    scale_fill_manual(values = Newcol) +
    xlim(0,100) + guides(fill='none') + clean_theme() 
  ggsave(paste0("landCglobal", year_cond[i],".png", sep = ""))
  
}




#..........................................................................................

#
#   Abundance
#

## Load species specific variance partition summarises: focusing on the normalized VP
##..........................................................................

load('VPcondwb_temp_ab_Rod.RDATA')
dim(VP.test1_rod$Vdiag)
dim(VP.test2_rod$Vdiag)

## conditional
VP.test_cond = array(NA, dim = dim(VP.test1_rod$Cnorm), dimnames = dimnames(VP.test1_rod$Cnorm))

for (i in 1:dim(VP.test1_rod$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_rod$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_rod$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_rod$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_rod$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_rod$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_rod$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.rod <- aperm(array(TabMF_rod$SR2, dim = c(8L, 3L, 3L, 2000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.rod)


Mean_fix_rod = array(NA, dim = c(3,3,dim(VP.test1_rod$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_rod$Vdiag)[5]){
  Mean_fix_rod[,,, k] = apply(varr.rod*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_rod$Cnorm), dimnames = dimnames(VP.test2_rod$Cnorm)) 

for (i in 1:dim(VP.test2_rod$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_rod$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_rod$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_rod$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_rod$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_rod$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_rod$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.rodbtw = apply(varr.rod*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.rodwit = apply(varr.rod*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


load('VPcondwb_temp_ab_bf.RDATA')
dim(VP.test1_bf$Vdiag)
dim(VP.test2_bf$Vdiag)

# conditional
VP.test_cond = array(NA, dim = dim(VP.test1_bf$Cnorm), dimnames = dimnames(VP.test1_bf$Cnorm))

for (i in 1:dim(VP.test1_bf$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_bf$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_bf$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_bf$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_bf$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_bf$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_bf$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.bf <- aperm(array(TabMF_bf$SR2, dim = c(57L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.bf)


Mean_fix_bf = array(NA, dim = c(3,3,dim(VP.test1_bf$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_bf$Vdiag)[5]){
  Mean_fix_bf[,,, k] = apply(varr.bf*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_bf$Cnorm), dimnames = dimnames(VP.test2_bf$Cnorm)) 

for (i in 1:dim(VP.test2_bf$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_bf$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_bf$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_bf$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_bf$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_bf$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_bf$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.bfbtw = apply(varr.bf*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.bfwit = apply(varr.bf*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


load('VPcondwb_temp_ab_moth.RDATA')
dim(VP.test1_moth$Vdiag)
dim(VP.test2_moth$Vdiag)

VP.test_cond = array(NA, dim = dim(VP.test1_moth$Cnorm), dimnames = dimnames(VP.test1_moth$Cnorm))

for (i in 1:dim(VP.test1_moth$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_moth$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_moth$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_moth$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_moth$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_moth$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_moth$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.moth <- aperm(array(TabMF_Moth$SR2, dim = c(319L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.moth)


Mean_fix_moth = array(NA, dim = c(3,3,dim(VP.test1_moth$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_moth$Vdiag)[5]){
  Mean_fix_moth[,,, k] = apply(varr.moth*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_moth$Cnorm), dimnames = dimnames(VP.test2_moth$Cnorm)) 

for (i in 1:dim(VP.test2_moth$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_moth$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_moth$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_moth$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_moth$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_moth$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_moth$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.mothbtw = apply(varr.moth*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.mothwit = apply(varr.moth*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 



load('VPcondwb_temp_ab_wg.RDATA')
dim(VP.test1_wg$Vdiag)
dim(VP.test2_wg$Vdiag)

# conditional
VP.test_cond = array(NA, dim = dim(VP.test1_wg$Cnorm), dimnames = dimnames(VP.test1_wg$Cnorm))

for (i in 1:dim(VP.test1_wg$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_wg$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_wg$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_wg$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_wg$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_wg$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_wg$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.wg <- aperm(array(TabMF_wg$SR2, dim = c(15L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.wg)


Mean_fix_wg = array(NA, dim = c(3,3,dim(VP.test1_wg$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_wg$Vdiag)[5]){
  Mean_fix_wg[,,, k] = apply(varr.wg*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}


## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_wg$Cnorm), dimnames = dimnames(VP.test2_wg$Cnorm)) 

for (i in 1:dim(VP.test2_wg$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_wg$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_wg$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_wg$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_wg$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_wg$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_wg$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.wgbtw = apply(varr.wg*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.wgwit = apply(varr.wg*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


load('VPcondwb_temp_ab_bd.RDATA')
dim(VP.test1_bd$Vdiag)
dim(VP.test2_bd$Vdiag)

# conditional
VP.test_cond = array(NA, dim = dim(VP.test1_bd$Cnorm), dimnames = dimnames(VP.test1_bd$Cnorm))

for (i in 1:dim(VP.test1_bd$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test1_bd$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test1_bd$Vdiag)[5]){ ## condition
      VP.test_cond[,,j,i,k] = VP.test1_bd$Cnorm[,,j,i,k]
      diag(VP.test_cond[,,j,i,k]) = VP.test1_bd$Vdiag[,,j,i,k]
      
      colnames( VP.test_cond[,,j,i,k]) = colnames(VP.test1_bd$Cnorm[,,j,i,k])
      rownames( VP.test_cond[,,j,i,k]) = rownames(VP.test1_bd$Cnorm[,,j,i,k])
      
    }  
  }
}


varr.bd <- aperm(array(TabMF_bd$SR2, dim = c(102L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr.bd)


Mean_fix_bd = array(NA, dim = c(3,3,dim(VP.test1_bd$Cnorm)[c(3,5)]))
for (k in 1:dim(VP.test1_bd$Vdiag)[5]){
  Mean_fix_bd[,,, k] = apply(varr.bd*VP.test_cond[1:3, 1:3,,, k], c(1, 2, 3), mean, na.rm=T)
}



# btw wit
## btw within
VP.test_bw = array(NA, dim = dim(VP.test2_bd$Cnorm), dimnames = dimnames(VP.test2_bd$Cnorm)) 

for (i in 1:dim(VP.test2_bd$Vdiag)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test2_bd$Vdiag)[3]){  ## species
    for (k in 1:dim(VP.test2_bd$Vdiag)[5]){ ## condition
      VP.test_bw[,,j,i,k] = VP.test2_bd$Cnorm[,,j,i,k]
      diag(VP.test_bw[,,j,i,k]) = VP.test2_bd$Vdiag[,,j,i,k]
      
      colnames(VP.test_bw[,,j,i,k]) = colnames(VP.test2_bd$Cnorm[,,j,i,k])
      rownames(VP.test_bw[,,j,i,k]) = rownames(VP.test2_bd$Cnorm[,,j,i,k])
      
    }  
  }
}



Mean.fix.bdbtw = apply(varr.bd*VP.test_bw[1:3, 1:3,,,'btw'], c(1, 2, 3), mean, na.rm=T) # 
Mean.fix.bdwit = apply(varr.bd*VP.test_bw[1:3, 1:3,,,'wit'], c(1, 2, 3), mean, na.rm=T) # 


####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggpubr)

VP.year = list()

for (k in 1:dim(VP.test1_bd$Vdiag)[5]){ # Same number of conditions year for all taxa
  Mean_fix_bd[,,, k] = Mean_fix_bd[,,, k]*100
  Mean_fix_bf[,,, k] = Mean_fix_bf[,,, k]*100
  Mean_fix_moth[,,, k] = Mean_fix_moth[,,, k]*100
  Mean_fix_rod[,,, k] = Mean_fix_rod[,,, k]*100
  Mean_fix_wg[,,, k] = Mean_fix_wg[,,, k]*100
  
  VP.rod = data.frame(Climate = Mean_fix_rod[1,1,, k], habprop = Mean_fix_rod[3,3,, k], 
                      habConf = Mean_fix_rod[2,2,, k])
  VP.rod$taxa = 'rodent'
  VP.rod$BG = dimnames(VP.test1_rod$Vdiag)$group[k]
  
  VP.bd = data.frame(Climate = Mean_fix_bd[1,1,, k], habprop = Mean_fix_bd[3,3,, k], 
                     habConf = Mean_fix_bd[2,2,, k])
  VP.bd$taxa = 'bird'
  VP.bd$BG = dimnames(VP.test1_bd$Vdiag)$group[k]
  
  VP.bf = data.frame(Climate = Mean_fix_bf[1,1,, k], habprop = Mean_fix_bf[3,3,, k], 
                     habConf = Mean_fix_bf[2,2,, k])
  VP.bf$taxa = 'butterfly'
  VP.bf$BG = dimnames(VP.test1_bf$Vdiag)$group[k]
  
  VP.moth = data.frame(Climate = Mean_fix_moth[1,1,, k], habprop = Mean_fix_moth[3,3,, k], 
                       habConf = Mean_fix_moth[2,2,, k])
  VP.moth$taxa = 'moth'
  VP.moth$BG = dimnames(VP.test1_moth$Vdiag)$group[k]
  
  
  VP.wg = data.frame(Climate = Mean_fix_wg[1,1,, k], habprop = Mean_fix_wg[3,3,, k], 
                     habConf = Mean_fix_wg[2,2,, k])
  VP.wg$taxa = 'large mammals'
  VP.wg$BG = dimnames(VP.test1_wg$Vdiag)$group[k]
  
  
  VP.year[[k]] = rbind(VP.rod, VP.bf, VP.moth, VP.bd, VP.wg)
}




Mean.fix.rodbtw = Mean.fix.rodbtw*100
Mean.fix.bfbtw = Mean.fix.bfbtw*100
Mean.fix.mothbtw = Mean.fix.mothbtw*100
Mean.fix.bdbtw = Mean.fix.bdbtw*100
Mean.fix.wgbtw = Mean.fix.wgbtw*100


VP.rodbtw = data.frame(Climate = Mean.fix.rodbtw[1,1,], habprop = Mean.fix.rodbtw[3,3,], 
                       habConf = Mean.fix.rodbtw[2,2,])

VP.rodbtw$taxa = 'rodent'

VP.bdbtw = data.frame(Climate = Mean.fix.bdbtw[1,1,], habprop = Mean.fix.bdbtw[3,3,], 
                      habConf = Mean.fix.bdbtw[2,2,])

VP.bdbtw$taxa = 'bird'

VP.bfbtw = data.frame(Climate = Mean.fix.bfbtw[1,1,], habprop = Mean.fix.bfbtw[3,3,], 
                      habConf = Mean.fix.bfbtw[2,2,])

VP.bfbtw$taxa = 'butterfly'

VP.mothbtw = data.frame(Climate = Mean.fix.mothbtw[1,1,], habprop = Mean.fix.mothbtw[3,3,], 
                        habConf = Mean.fix.mothbtw[2,2,])

VP.mothbtw$taxa = 'moth'

VP.wgbtw = data.frame(Climate = Mean.fix.wgbtw[1,1,], habprop = Mean.fix.wgbtw[3,3,], 
                      habConf = Mean.fix.wgbtw[2,2,])

VP.wgbtw$taxa = 'large mammals'


VP.btw = rbind(VP.rodbtw, VP.bfbtw, VP.mothbtw, VP.bdbtw, VP.wgbtw)
VP.btw$BG = 'btw'



Mean.fix.rodwit = Mean.fix.rodwit*100
Mean.fix.bfwit = Mean.fix.bfwit*100
Mean.fix.mothwit = Mean.fix.mothwit*100
Mean.fix.bdwit = Mean.fix.bdwit*100
Mean.fix.wgwit = Mean.fix.wgwit*100

VP.rodwit = data.frame(Climate = Mean.fix.rodwit[1,1,], habprop = Mean.fix.rodwit[3,3,], 
                       habConf = Mean.fix.rodwit[2,2,])

VP.rodwit$taxa = 'rodent'

VP.bdwit = data.frame(Climate = Mean.fix.bdwit[1,1,], habprop = Mean.fix.bdwit[3,3,], 
                      habConf = Mean.fix.bdwit[2,2,])

VP.bdwit$taxa = 'bird'

VP.bfwit = data.frame(Climate = Mean.fix.bfwit[1,1,], habprop = Mean.fix.bfwit[3,3,], 
                      habConf = Mean.fix.bfwit[2,2,])

VP.bfwit$taxa = 'butterfly'


VP.mothwit = data.frame(Climate = Mean.fix.mothwit[1,1,], habprop = Mean.fix.mothwit[3,3,], 
                        habConf = Mean.fix.mothwit[2,2,])

VP.mothwit$taxa = 'moth'

VP.wgwit = data.frame(Climate = Mean.fix.wgwit[1,1,], habprop = Mean.fix.wgwit[3,3,], 
                      habConf = Mean.fix.wgwit[2,2,])

VP.wgwit$taxa = 'large mammals'


VP.wit = rbind(VP.rodwit, VP.bfwit, VP.mothwit, VP.bdwit, VP.wgwit)
VP.wit$BG = 'wit'


VP.year.array = do.call(rbind, VP.year)

VP.sp = rbind(VP.btw, VP.wit, VP.year.array)
VP.sp$BG = factor(VP.sp$BG, levels = c('1999', '2000', '2001', '2002', '2003', '2004', '2005','2006',
                                       '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', 
                                       '2015', '2016', '2017', '2018', '2019','btw', 'wit'))
VP.sp$taxa = factor(VP.sp$taxa, levels = c('bird', 'butterfly', 'moth', 'rodent', 'large mammals'))


# colors
library(wesanderson)
names(wes_palettes)


library(grid)
head(VP.sp)
summary(VP.sp)

# 
# VP.sp2 = VP.sp[VP.sp$Climate<100 & VP.sp$Climate>0,]
# VP.sp3 = VP.sp2[VP.sp2$habprop<100 & VP.sp2$habprop>0,]
# VP.sp = VP.sp3[VP.sp3$habConf<100 & VP.sp3$habConf>0,]
# 

VP.sp$sum = VP.sp$Climate + VP.sp$habprop + VP.sp$habConf

VP.sp$Climate = VP.sp$Climate*100/VP.sp$sum
VP.sp$habprop = VP.sp$habprop*100/VP.sp$sum
VP.sp$habConf = VP.sp$habConf*100/VP.sp$sum

summary(VP.sp)


Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")
# 

year_cond = levels(VP.sp$BG)

for (i in 1:length(year_cond)) {
  VP.sp.cond = VP.sp[VP.sp$BG == year_cond[i],]
  
  VP.spmean = VP.sp.cond[, c(1:4)] %>%
    group_by(taxa) %>%
    #summarise(across(everything(), list(mean)))
    summarise(across(everything(), list(mean)))
  colnames(VP.spmean)[2:4] = c('Climate', 'habprop', 'habConf')
  
  
  
  tern = ggtern(VP.sp.cond, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
    #geom_point()+
    stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                      position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
    #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
    scale_alpha_continuous(range = c(0.1, 0.3)) +
    scale_fill_manual(values = Newcol) +
    #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
    scale_color_manual(values = Newcol) +
    geom_crosshair_tern(data = VP.spmean, lty = 2, linewidth = 3.2) +  # add mean information and lines
    geom_point(data = VP.spmean, aes(color = taxa), show.legend = F) + 
    #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
    #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
    #theme_showarrows() +
    theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
          tern.axis.arrow.T = element_line(size=0, color="white"),
          tern.axis.arrow.R = element_line(size=0, color="white"),
          tern.axis.line.L = element_line(color='orange3',size=2),
          tern.axis.line.T = element_line(color='#66C2A5',size=2),
          tern.axis.line.R = element_line(color='yellow3',size=2),
          axis.text = element_text(size = 28),
          axis.title = element_text(size = 28)) + 
    guides(alpha='none', colour='none') + 
    labs(title = year_cond[i], x='', y='', z='', 
         xarrow = "Habitat proportions",
         yarrow = "Climate",
         zarrow = "Habitat configuration")
  ggsave(tern, file=paste0("VTern_ab", year_cond[i],".png", sep = ""))
  
  climplot <- ggdensity(VP.sp.cond, 'Climate', fill = 'taxa') +
    scale_fill_manual(values = Newcol) +
    scale_y_reverse()+ scale_x_reverse() +
    xlim(100,0) + guides(fill='none') + clean_theme() 
  ggsave(paste0("climglobalab", year_cond[i],".png", sep = ""))
  
  habPplot <- ggdensity(VP.sp.cond, "habprop",  fill = 'taxa') +
    scale_fill_manual(values = Newcol) +
    xlim(0,100) + guides(fill='none') + clean_theme() 
  ggsave(paste0("habPglobalab", year_cond[i],".png", sep = ""))
  
  habCplot <- ggdensity(VP.sp.cond, "habConf",  fill = 'taxa') +
    scale_fill_manual(values = Newcol) +
    xlim(0,100) + guides(fill='none') + clean_theme() 
  ggsave(paste0("landCglobalab", year_cond[i],".png", sep = ""))
  
}

