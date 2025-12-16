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

setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")

#
#   Occurrence
#

## Load species specific variance partition summarises: focusing on the normalized VP
##..........................................................................

load('VPcondwb_Rod.RDATA')
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


varr.rod <- aperm(array(TabMF_rod$TjurR2, dim = c(1000L, 8L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.rod)

Mean.fix.rodNB = apply(varr.rod*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.rodMB = apply(varr.rod*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.rodSB = apply(varr.rod*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB


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


load('VPcondwb_bf.RDATA')
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


varr.bf <- aperm(array(TabMF_bf$TjurR2, dim = c(1000L, 57L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.bf)

Mean.fix.bfNB = apply(varr.bf*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.bfMB = apply(varr.bf*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.bfSB = apply(varr.bf*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB

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


load('VPcondwb_moth.RDATA')
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


varr.moth <- aperm(array(TabMF_Moth$TjurR2, dim = c(1000L, 319L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.moth)

Mean.fix.mothNB = apply(varr.moth*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.mothMB = apply(varr.moth*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.mothSB = apply(varr.moth*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB

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



load('VPcondwb_wg.RDATA')
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


varr.wg <- aperm(array(TabMF_wg$TjurR2, dim = c(1000L, 15L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.wg)

Mean.fix.wgNB = apply(varr.wg*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.wgMB = apply(varr.wg*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.wgSB = apply(varr.wg*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB

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


load('VPcondwb_bd.RDATA')
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


varr.bd <- aperm(array(TabMF_bd$TjurR2, dim = c(1000L, 102L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.bd)

Mean.fix.bdNB = apply(varr.bd*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.bdMB = apply(varr.bd*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.bdSB = apply(varr.bd*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB


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

Mean.fix.rodNB = Mean.fix.rodNB*100
Mean.fix.bfNB = Mean.fix.bfNB*100
Mean.fix.mothNB = Mean.fix.mothNB*100
Mean.fix.bdNB = Mean.fix.bdNB*100
Mean.fix.wgNB = Mean.fix.wgNB*100


VP.rodNB = data.frame(Climate = Mean.fix.rodNB[1,1,], habprop = Mean.fix.rodNB[3,3,], 
                        habConf = Mean.fix.rodNB[2,2,])
VP.rodNB$taxa = 'rodent'

VP.bdNB = data.frame(Climate = Mean.fix.bdNB[1,1,], habprop = Mean.fix.bdNB[3,3,], 
                       habConf = Mean.fix.bdNB[2,2,])
VP.bdNB$taxa = 'bird'

VP.bfNB = data.frame(Climate = Mean.fix.bfNB[1,1,], habprop = Mean.fix.bfNB[3,3,], 
                       habConf = Mean.fix.bfNB[2,2,])


VP.bfNB$taxa = 'butterfly'

VP.mothNB = data.frame(Climate = Mean.fix.mothNB[1,1,], habprop = Mean.fix.mothNB[3,3,], 
                         habConf = Mean.fix.mothNB[2,2,])

VP.mothNB$taxa = 'moth'


VP.wgNB = data.frame(Climate = Mean.fix.wgNB[1,1,], habprop = Mean.fix.wgNB[3,3,], 
                       habConf = Mean.fix.wgNB[2,2,])

VP.wgNB$taxa = 'large mammals'


VP.NB = rbind(VP.rodNB, VP.bfNB, VP.mothNB, VP.bdNB, VP.wgNB)
VP.NB$BG = 'NB'



Mean.fix.rodMB = Mean.fix.rodMB*100
Mean.fix.bfMB = Mean.fix.bfMB*100
Mean.fix.mothMB = Mean.fix.mothMB*100
Mean.fix.bdMB = Mean.fix.bdMB*100
Mean.fix.wgMB = Mean.fix.wgMB*100

VP.rodMB = data.frame(Climate = Mean.fix.rodMB[1,1,], habprop = Mean.fix.rodMB[3,3,], 
                        habConf = Mean.fix.rodMB[2,2,])

VP.rodMB$taxa = 'rodent'

VP.bdMB = data.frame(Climate = Mean.fix.bdMB[1,1,], habprop = Mean.fix.bdMB[3,3,], 
                       habConf = Mean.fix.bdMB[2,2,])
VP.bdMB$taxa = 'bird'

VP.bfMB = data.frame(Climate = Mean.fix.bfMB[1,1,], habprop = Mean.fix.bfMB[3,3,], 
                       habConf = Mean.fix.bfMB[2,2,])

VP.bfMB$taxa = 'butterfly'

VP.mothMB = data.frame(Climate = Mean.fix.mothMB[1,1,], habprop = Mean.fix.mothMB[3,3,], 
                         habConf = Mean.fix.mothMB[2,2,])

VP.mothMB$taxa = 'moth'

VP.wgMB = data.frame(Climate = Mean.fix.wgMB[1,1,], habprop = Mean.fix.wgMB[3,3,], 
                       habConf = Mean.fix.wgMB[2,2,])

VP.wgMB$taxa = 'large mammals'


VP.MB = rbind(VP.rodMB, VP.bfMB, VP.mothMB, VP.bdMB, VP.wgMB)
VP.MB$BG = 'MB'



Mean.fix.rodSB = Mean.fix.rodSB*100
Mean.fix.bfSB = Mean.fix.bfSB*100
Mean.fix.mothSB = Mean.fix.mothSB*100
Mean.fix.bdSB = Mean.fix.bdSB*100
Mean.fix.wgSB = Mean.fix.wgSB*100

VP.rodSB = data.frame(Climate = Mean.fix.rodSB[1,1,], habprop = Mean.fix.rodSB[3,3,], 
                        habConf = Mean.fix.rodSB[2,2,])

VP.rodSB$taxa = 'rodent'

VP.bdSB = data.frame(Climate = Mean.fix.bdSB[1,1,], habprop = Mean.fix.bdSB[3,3,], 
                       habConf = Mean.fix.bdSB[2,2,])

VP.bdSB$taxa = 'bird'

VP.bfSB = data.frame(Climate = Mean.fix.bfSB[1,1,], habprop = Mean.fix.bfSB[3,3,], 
                       habConf = Mean.fix.bfSB[2,2,])

VP.bfSB$taxa = 'butterfly'

VP.mothSB = data.frame(Climate = Mean.fix.mothSB[1,1,], habprop = Mean.fix.mothSB[3,3,], 
                         habConf = Mean.fix.mothSB[2,2,])

VP.mothSB$taxa = 'moth'

VP.wgSB = data.frame(Climate = Mean.fix.wgSB[1,1,], habprop = Mean.fix.wgSB[3,3,], 
                       habConf = Mean.fix.wgSB[2,2,])

VP.wgSB$taxa = 'large mammals'


VP.SB = rbind(VP.rodSB, VP.bfSB, VP.mothSB, VP.bdSB, VP.wgSB)
VP.SB$BG = 'SB'




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



VP.sp = rbind(VP.btw, VP.wit, VP.SB, VP.MB, VP.NB)
VP.sp$BG = factor(VP.sp$BG, levels = c('NB', 'MB', 'SB','btw', 'wit'))
VP.sp$taxa = factor(VP.sp$taxa, levels = c('bird', 'butterfly', 'moth', 'rodent', 'large mammals'))


# colors
library(wesanderson)
names(wes_palettes)


library(grid)
head(VP.sp)
summary(VP.sp)

VP.sp2 = VP.sp[VP.sp$Climate<100 & VP.sp$Climate>0,]
VP.sp3 = VP.sp2[VP.sp2$habprop<100 & VP.sp2$habprop>0,]
VP.sp = VP.sp3[VP.sp3$habConf<100 & VP.sp3$habConf>0,]


VP.sp$sum = VP.sp$Climate + VP.sp$habprop + VP.sp$habConf

VP.sp$Climate = VP.sp$Climate*100/VP.sp$sum
VP.sp$habprop = VP.sp$habprop*100/VP.sp$sum
VP.sp$habConf = VP.sp$habConf*100/VP.sp$sum

summary(VP.sp)

VP.sp.NB = VP.sp[VP.sp$BG == 'NB',]

VP.spmeanNB = VP.sp.NB[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanNB)[2:4] = c('Climate', 'habprop', 'habConf')


Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")

ternNB = ggtern(VP.sp.NB, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanNB, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanNB, aes(color = taxa), show.legend = F) + 
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
  labs(title = 'North Boreal', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternNB, file=paste0("VTernNB_occ",".png", sep = ""))

climplot <- ggdensity(VP.sp.NB, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_NB.png")

habPplot <- ggdensity(VP.sp.NB, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_NB.png")

habCplot <- ggdensity(VP.sp.NB, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_NB.png")




VP.sp.MB = VP.sp[VP.sp$BG == 'MB',]

VP.spmeanMB = VP.sp.MB[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanMB)[2:4] = c('Climate', 'habprop', 'habConf')

ternMB = ggtern(VP.sp.MB, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanMB, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanMB, aes(color = taxa), show.legend = F) + 
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
  labs(title = 'Middle Boreal', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternMB, file=paste0("VTernMB_occ",".png", sep = ""))

climplot <- ggdensity(VP.sp.MB, 'Climate', fill = 'taxa') +
  scale_fill_manual(values =Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_MB.png")

habPplot <- ggdensity(VP.sp.MB, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_MB.png")

habCplot <- ggdensity(VP.sp.MB, "habConf",  fill = 'taxa') +
  scale_fill_manual(values =Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_MB.png")



VP.sp.SB = VP.sp[VP.sp$BG == 'SB',]

VP.spmeanSB = VP.sp.SB[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanSB)[2:4] = c('Climate', 'habprop', 'habConf')

ternSB = ggtern(VP.sp.SB, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanSB, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanSB, aes(color = taxa), show.legend = F) + 
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
  labs(title = 'South Boreal', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternSB, file=paste0("VTernSB_occ",".png", sep = ""))

climplot <- ggdensity(VP.sp.SB, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_SB.png")

habPplot <- ggdensity(VP.sp.SB, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_SB.png")

habCplot <- ggdensity(VP.sp.SB, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_SB.png")





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




ggtern::grid.arrange(ternNB, ternMB, ternSB, ternbtw, ternwit,
             ncol=3, nrow=2)







#
#   Abundance
#

## Load species specific variance partition summarises: focusing on the normalized VP
##..........................................................................

load('VPcondwb_ab_Rod.RDATA')
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


varr.rod <- aperm(array(TabMF_rod$SR2, dim = c(2000L, 8L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.rod)

Mean.fix.rodNB = apply(varr.rod*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.rodMB = apply(varr.rod*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.rodSB = apply(varr.rod*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB


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


load('VPcondwb_ab_bf.RDATA')
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


varr.bf <- aperm(array(TabMF_bf$SR2, dim = c(1000L, 57L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.bf)

Mean.fix.bfNB = apply(varr.bf*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.bfMB = apply(varr.bf*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.bfSB = apply(varr.bf*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB

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


load('VPcondwb_ab_moth.RDATA')
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


varr.moth <- aperm(array(TabMF_Moth$SR2, dim = c(1000L, 319L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.moth)

Mean.fix.mothNB = apply(varr.moth*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.mothMB = apply(varr.moth*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.mothSB = apply(varr.moth*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB

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



load('VPcondwb_ab_wg.RDATA')
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


varr.wg <- aperm(array(TabMF_wg$SR2, dim = c(1000L, 15L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.wg)

Mean.fix.wgNB = apply(varr.wg*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.wgMB = apply(varr.wg*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.wgSB = apply(varr.wg*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB

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


load('VPcondwb_ab_bd.RDATA')
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


varr.bd <- aperm(array(TabMF_bd$SR2, dim = c(1000L, 102L, 3L, 3L)), perm = c(4L, 3L, 2L, 1L))
dim(varr.bd)

Mean.fix.bdNB = apply(varr.bd*VP.test_cond[1:3, 1:3,,,'NB'], c(1, 2, 3), mean, na.rm=T) # NB
Mean.fix.bdMB = apply(varr.bd*VP.test_cond[1:3, 1:3,,,'MB'], c(1, 2, 3), mean, na.rm=T) # MB
Mean.fix.bdSB = apply(varr.bd*VP.test_cond[1:3, 1:3,,,'SB'], c(1, 2, 3), mean, na.rm=T) # SB


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

Mean.fix.rodNB = Mean.fix.rodNB*100
Mean.fix.bfNB = Mean.fix.bfNB*100
Mean.fix.mothNB = Mean.fix.mothNB*100
Mean.fix.bdNB = Mean.fix.bdNB*100
Mean.fix.wgNB = Mean.fix.wgNB*100


VP.rodNB = data.frame(Climate = Mean.fix.rodNB[1,1,], habprop = Mean.fix.rodNB[3,3,], 
                      habConf = Mean.fix.rodNB[2,2,])
VP.rodNB$taxa = 'rodent'

VP.bdNB = data.frame(Climate = Mean.fix.bdNB[1,1,], habprop = Mean.fix.bdNB[3,3,], 
                     habConf = Mean.fix.bdNB[2,2,])
VP.bdNB$taxa = 'bird'

VP.bfNB = data.frame(Climate = Mean.fix.bfNB[1,1,], habprop = Mean.fix.bfNB[3,3,], 
                     habConf = Mean.fix.bfNB[2,2,])


VP.bfNB$taxa = 'butterfly'

VP.mothNB = data.frame(Climate = Mean.fix.mothNB[1,1,], habprop = Mean.fix.mothNB[3,3,], 
                       habConf = Mean.fix.mothNB[2,2,])

VP.mothNB$taxa = 'moth'


VP.wgNB = data.frame(Climate = Mean.fix.wgNB[1,1,], habprop = Mean.fix.wgNB[3,3,], 
                     habConf = Mean.fix.wgNB[2,2,])

VP.wgNB$taxa = 'large mammals'


VP.NB = rbind(VP.rodNB, VP.bfNB, VP.mothNB, VP.bdNB, VP.wgNB)
VP.NB$BG = 'NB'



Mean.fix.rodMB = Mean.fix.rodMB*100
Mean.fix.bfMB = Mean.fix.bfMB*100
Mean.fix.mothMB = Mean.fix.mothMB*100
Mean.fix.bdMB = Mean.fix.bdMB*100
Mean.fix.wgMB = Mean.fix.wgMB*100

VP.rodMB = data.frame(Climate = Mean.fix.rodMB[1,1,], habprop = Mean.fix.rodMB[3,3,], 
                      habConf = Mean.fix.rodMB[2,2,])

VP.rodMB$taxa = 'rodent'

VP.bdMB = data.frame(Climate = Mean.fix.bdMB[1,1,], habprop = Mean.fix.bdMB[3,3,], 
                     habConf = Mean.fix.bdMB[2,2,])
VP.bdMB$taxa = 'bird'

VP.bfMB = data.frame(Climate = Mean.fix.bfMB[1,1,], habprop = Mean.fix.bfMB[3,3,], 
                     habConf = Mean.fix.bfMB[2,2,])

VP.bfMB$taxa = 'butterfly'

VP.mothMB = data.frame(Climate = Mean.fix.mothMB[1,1,], habprop = Mean.fix.mothMB[3,3,], 
                       habConf = Mean.fix.mothMB[2,2,])

VP.mothMB$taxa = 'moth'

VP.wgMB = data.frame(Climate = Mean.fix.wgMB[1,1,], habprop = Mean.fix.wgMB[3,3,], 
                     habConf = Mean.fix.wgMB[2,2,])

VP.wgMB$taxa = 'large mammals'


VP.MB = rbind(VP.rodMB, VP.bfMB, VP.mothMB, VP.bdMB, VP.wgMB)
VP.MB$BG = 'MB'



Mean.fix.rodSB = Mean.fix.rodSB*100
Mean.fix.bfSB = Mean.fix.bfSB*100
Mean.fix.mothSB = Mean.fix.mothSB*100
Mean.fix.bdSB = Mean.fix.bdSB*100
Mean.fix.wgSB = Mean.fix.wgSB*100

VP.rodSB = data.frame(Climate = Mean.fix.rodSB[1,1,], habprop = Mean.fix.rodSB[3,3,], 
                      habConf = Mean.fix.rodSB[2,2,])

VP.rodSB$taxa = 'rodent'

VP.bdSB = data.frame(Climate = Mean.fix.bdSB[1,1,], habprop = Mean.fix.bdSB[3,3,], 
                     habConf = Mean.fix.bdSB[2,2,])

VP.bdSB$taxa = 'bird'

VP.bfSB = data.frame(Climate = Mean.fix.bfSB[1,1,], habprop = Mean.fix.bfSB[3,3,], 
                     habConf = Mean.fix.bfSB[2,2,])

VP.bfSB$taxa = 'butterfly'

VP.mothSB = data.frame(Climate = Mean.fix.mothSB[1,1,], habprop = Mean.fix.mothSB[3,3,], 
                       habConf = Mean.fix.mothSB[2,2,])

VP.mothSB$taxa = 'moth'

VP.wgSB = data.frame(Climate = Mean.fix.wgSB[1,1,], habprop = Mean.fix.wgSB[3,3,], 
                     habConf = Mean.fix.wgSB[2,2,])

VP.wgSB$taxa = 'large mammals'


VP.SB = rbind(VP.rodSB, VP.bfSB, VP.mothSB, VP.bdSB, VP.wgSB)
VP.SB$BG = 'SB'




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



VP.sp = rbind(VP.btw, VP.wit, VP.SB, VP.MB, VP.NB)
VP.sp$BG = factor(VP.sp$BG, levels = c('NB', 'MB', 'SB','btw', 'wit'))
VP.sp$taxa = factor(VP.sp$taxa, levels = c('bird', 'butterfly', 'moth', 'rodent', 'large mammals'))


# colors
library(wesanderson)
names(wes_palettes)


library(grid)
head(VP.sp)
summary(VP.sp)

VP.sp2 = VP.sp[VP.sp$Climate<100 & VP.sp$Climate>0,]
VP.sp3 = VP.sp2[VP.sp2$habprop<100 & VP.sp2$habprop>0,]
VP.sp = VP.sp3[VP.sp3$habConf<100 & VP.sp3$habConf>0,]


VP.sp$sum = VP.sp$Climate + VP.sp$habprop + VP.sp$habConf

VP.sp$Climate = VP.sp$Climate*100/VP.sp$sum
VP.sp$habprop = VP.sp$habprop*100/VP.sp$sum
VP.sp$habConf = VP.sp$habConf*100/VP.sp$sum

summary(VP.sp)

VP.sp.NB = VP.sp[VP.sp$BG == 'NB',]

VP.spmeanNB = VP.sp.NB[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanNB)[2:4] = c('Climate', 'habprop', 'habConf')


Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")

ternNB = ggtern(VP.sp.NB, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanNB, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanNB, aes(color = taxa), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"
        ),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 28)) + 
  guides(alpha='none', colour='none') + 
  labs(title = 'North Boreal', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternNB, file=paste0("VTernNB_ab",".png", sep = ""))

climplot <- ggdensity(VP.sp.NB, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_ab_NB.png")

habPplot <- ggdensity(VP.sp.NB, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_ab_NB.png")

habCplot <- ggdensity(VP.sp.NB, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_ab_NB.png")




VP.sp.MB = VP.sp[VP.sp$BG == 'MB',]

VP.spmeanMB = VP.sp.MB[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanMB)[2:4] = c('Climate', 'habprop', 'habConf')

ternMB = ggtern(VP.sp.MB, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanMB, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanMB, aes(color = taxa), show.legend = F) + 
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
  labs(title = 'Middle Boreal', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternMB, file=paste0("VTernMB_ab",".png", sep = ""))

climplot <- ggdensity(VP.sp.MB, 'Climate', fill = 'taxa') +
  scale_fill_manual(values =Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_ab_MB.png")

habPplot <- ggdensity(VP.sp.MB, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_ab_MB.png")

habCplot <- ggdensity(VP.sp.MB, "habConf",  fill = 'taxa') +
  scale_fill_manual(values =Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_ab_MB.png")



VP.sp.SB = VP.sp[VP.sp$BG == 'SB',]

VP.spmeanSB = VP.sp.SB[, c(1:4)] %>%
  group_by(taxa) %>%
  #summarise(across(everything(), list(mean)))
  summarise(across(everything(), list(mean)))
colnames(VP.spmeanSB)[2:4] = c('Climate', 'habprop', 'habConf')

ternSB = ggtern(VP.sp.SB, aes(habprop, Climate, habConf, fill = taxa, color=taxa)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = taxa), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = Newcol) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = Newcol) +
  geom_crosshair_tern(data = VP.spmeanSB, lty = 2, size = 3.2) +  # add mean information and lines
  geom_point(data = VP.spmeanSB, aes(color = taxa), show.legend = F) + 
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
  labs(title = 'South Boreal', x='', y='', z='', 
       xarrow = "Habitat proportions",
       yarrow = "Climate",
       zarrow = "Habitat configuration")
ggsave(ternSB, file=paste0("VTernSB_ab",".png", sep = ""))

climplot <- ggdensity(VP.sp.SB, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_ab_SB.png")

habPplot <- ggdensity(VP.sp.SB, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_ab_SB.png")

habCplot <- ggdensity(VP.sp.SB, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_ab_SB.png")





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
ggsave(ternbtw, file=paste0("VTernBTW_ab",".png", sep = ""))

climplot <- ggdensity(VP.sp.btw, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_ab_btw.png")

habPplot <- ggdensity(VP.sp.btw, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_ab_btw.png")

habCplot <- ggdensity(VP.sp.btw, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_ab_btw.png")



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
ggsave(ternwit, file=paste0("VTernWIT_ab",".png", sep = ""))


climplot <- ggdensity(VP.sp.wit, 'Climate', fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
ggsave("climglobal_ab_wit.png")

habPplot <- ggdensity(VP.sp.wit, "habprop",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("habPglobal_ab_wit.png")

habCplot <- ggdensity(VP.sp.wit, "habConf",  fill = 'taxa') +
  scale_fill_manual(values = Newcol) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_ab_wit.png")




ggtern::grid.arrange(ternNB, ternMB, ternSB, ternbtw, ternwit,
                     ncol=3, nrow=2)


,100) + guides(fill='none') + clean_theme() 
ggsave("landCglobal_ab_wit.png")




ggtern::grid.arrange(ternNB, ternMB, ternSB, ternbtw, ternwit,
                     ncol=3, nrow=2)

)

