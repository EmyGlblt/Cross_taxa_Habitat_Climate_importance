#
#             Summary plots to compare taxa
#
#
# --------------------------------------------------------------------


library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)
library(grid)
library(ggpubr)
library(ggeasy)
library(tidyr)

##...............................................  occurrence data

setwd("C:/Users/guilbaul/OneDrive - University of Helsinki/RECcII/ScriptData/Data/4.VPanalyses/VP_cond_trait")

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

VPfix.rodNB = varr.rod*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.rodMB = varr.rod*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.rodSB = varr.rod*VP.test_cond[1:3, 1:3,,,'SB'] # SB


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



VPfix.rodbtw = varr.rod*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.rodwit = varr.rod*VP.test_bw[1:3, 1:3,,,'wit'] # 


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

VPfix.bfNB = varr.bf*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.bfMB = varr.bf*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.bfSB = varr.bf*VP.test_cond[1:3, 1:3,,,'SB'] # SB

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



VPfix.bfbtw = varr.bf*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.bfwit = varr.bf*VP.test_bw[1:3, 1:3,,,'wit'] # 


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

VPfix.mothNB = varr.moth*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.mothMB = varr.moth*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.mothsB = varr.moth*VP.test_cond[1:3, 1:3,,,'SB'] # SB

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



VPfix.mothbtw = varr.moth*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.mothwit = varr.moth*VP.test_bw[1:3, 1:3,,,'wit'] # 



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

VPfix.wgNB = varr.wg*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.wgMB = varr.wg*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.wgSB = varr.wg*VP.test_cond[1:3, 1:3,,,'SB'] # SB

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



VPfix.wgbtw = varr.wg*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.wgwit = varr.wg*VP.test_bw[1:3, 1:3,,,'wit'] # 


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

VPfix.bdNB = varr.bd*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.bdMB = varr.bd*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.bdSB = varr.bd*VP.test_cond[1:3, 1:3,,,'SB'] # SB


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



VPfix.bdbtw = varr.bd*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.bdwit = varr.bd*VP.test_bw[1:3, 1:3,,,'wit'] # 




##.............................  Summarize


spvp_summaries = function(vpdat, Threshold){
  
  # Initialize a matrix to store the count of being the maximum
  # rows are the fixed effect only but columns are fixed effect and random effects we need to focus only on fixed egffect
  max_count <- array(0, dim = c(dim(vpdat)[1], dim(vpdat)[3]))
  #dim(max_count)  # fix effect # species # posterior sample
  
  # For each species and posterior sample, find the covariate with the maximum value
  for (s in 1:dim(vpdat)[3]) {
    for (p in 1:dim(vpdat)[4]) {
      # Find the index of the maximum value for this species and sample
      max_val <- max(diag(vpdat[,,s,p]))
      max_indices <- which(diag(vpdat[,,s,p]) == max_val, arr = TRUE)
      
      max_count[max_indices, s] <- max_count[max_indices, s] + 1
    }
  }
  
  # Calculate the proportion of posterior samples where each covariate is the maximum for each species
  # calculate among the species that are not NA (R2) so using ColSums instead of dim(Vpdat[4])
  prop_max <- max_count / colSums(max_count)
  
  
  prop_sp <- rowMeans(prop_max >= Threshold)
  
  # prop_species now contains the proportion for each covariate
  print(prop_sp)
  
  return(list(prop_sp = prop_sp))
}

# different posterior probabilities

## NB
# 90%
max90.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.90)
max90.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.90)
max90.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.90)
max90.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.90)
max90.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.90)


# 95%
max95.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.95)
max95.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.95)
max95.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.95)
max95.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.95)
max95.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.95)


# 80%
max80.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.80)
max80.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.80)
max80.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.80)
max80.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.80)
max80.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.80)

# 75%
max75.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.75)
max75.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.75)
max75.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.75)
max75.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.75)
max75.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.75)


maxvpNB = data.frame(rbind(rbind(max90.rodNB$prop_sp, max90.bfNB$prop_sp, max90.mothNB$prop_sp, max90.bdNB$prop_sp, max90.wgNB$prop_sp),
                         rbind(max95.rodNB$prop_sp, max95.bfNB$prop_sp, max95.mothNB$prop_sp, max95.bdNB$prop_sp, max95.wgNB$prop_sp),
                         rbind(max80.rodNB$prop_sp, max80.bfNB$prop_sp, max80.mothNB$prop_sp, max80.bdNB$prop_sp, max80.wgNB$prop_sp),
                         rbind(max75.rodNB$prop_sp, max75.bfNB$prop_sp, max75.mothNB$prop_sp, max75.bdNB$prop_sp, max75.wgNB$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                   cond = 'NB')
colnames(maxvpNB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## MB
# 90%
max90.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.90)
max90.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.90)
max90.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.90)
max90.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.90)
max90.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.90)


# 95%
max95.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.95)
max95.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.95)
max95.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.95)
max95.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.95)
max95.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.95)


# 80%
max80.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.80)
max80.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.80)
max80.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.80)
max80.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.80)
max80.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.80)

# 75%
max75.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.75)
max75.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.75)
max75.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.75)
max75.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.75)
max75.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.75)


maxvpMB = data.frame(rbind(rbind(max90.rodMB$prop_sp, max90.bfMB$prop_sp, max90.mothMB$prop_sp, max90.bdMB$prop_sp, max90.wgMB$prop_sp),
                         rbind(max95.rodMB$prop_sp, max95.bfMB$prop_sp, max95.mothMB$prop_sp, max95.bdMB$prop_sp, max95.wgMB$prop_sp),
                         rbind(max80.rodMB$prop_sp, max80.bfMB$prop_sp, max80.mothMB$prop_sp, max80.bdMB$prop_sp, max80.wgMB$prop_sp),
                         rbind(max75.rodMB$prop_sp, max75.bfMB$prop_sp, max75.mothMB$prop_sp, max75.bdMB$prop_sp, max75.wgMB$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                   cond = 'MB')
colnames(maxvpMB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## SB
# 90%
max90.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.90)
max90.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.90)
max90.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.90)
max90.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.90)
max90.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.90)


# 95%
max95.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.95)
max95.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.95)
max95.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.95)
max95.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.95)
max95.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.95)


# 80%
max80.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.80)
max80.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.80)
max80.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.80)
max80.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.80)
max80.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.80)

# 75%
max75.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.75)
max75.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.75)
max75.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.75)
max75.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.75)
max75.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.75)


maxvpSB = data.frame(rbind(rbind(max90.rodSB$prop_sp, max90.bfSB$prop_sp, max90.mothsB$prop_sp, max90.bdSB$prop_sp, max90.wgSB$prop_sp),
                         rbind(max95.rodSB$prop_sp, max95.bfSB$prop_sp, max95.mothsB$prop_sp, max95.bdSB$prop_sp, max95.wgSB$prop_sp),
                         rbind(max80.rodSB$prop_sp, max80.bfSB$prop_sp, max80.mothsB$prop_sp, max80.bdSB$prop_sp, max80.wgSB$prop_sp),
                         rbind(max75.rodSB$prop_sp, max75.bfSB$prop_sp, max75.mothsB$prop_sp, max75.bdSB$prop_sp, max75.wgSB$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                   cond = 'SB')
colnames(maxvpSB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## wit
# 90%
max90.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.90)
max90.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.90)
max90.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.90)
max90.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.90)
max90.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.90)


# 95%
max95.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.95)
max95.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.95)
max95.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.95)
max95.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.95)
max95.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.95)


# 80%
max80.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.80)
max80.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.80)
max80.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.80)
max80.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.80)
max80.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.80)

# 75%
max75.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.75)
max75.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.75)
max75.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.75)
max75.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.75)
max75.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.75)


maxvpwit = data.frame(rbind(rbind(max90.rodwit$prop_sp, max90.bfwit$prop_sp, max90.mothwit$prop_sp, max90.bdwit$prop_sp, max90.wgwit$prop_sp),
                         rbind(max95.rodwit$prop_sp, max95.bfwit$prop_sp, max95.mothwit$prop_sp, max95.bdwit$prop_sp, max95.wgwit$prop_sp),
                         rbind(max80.rodwit$prop_sp, max80.bfwit$prop_sp, max80.mothwit$prop_sp, max80.bdwit$prop_sp, max80.wgwit$prop_sp),
                         rbind(max75.rodwit$prop_sp, max75.bfwit$prop_sp, max75.mothwit$prop_sp, max75.bdwit$prop_sp, max75.wgwit$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                   cond = 'wit')
colnames(maxvpwit)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## btw
# 90%
max90.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.90)
max90.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.90)
max90.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.90)
max90.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.90)
max90.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.90)


# 95%
max95.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.95)
max95.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.95)
max95.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.95)
max95.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.95)
max95.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.95)


# 80%
max80.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.80)
max80.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.80)
max80.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.80)
max80.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.80)
max80.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.80)

# 75%
max75.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.75)
max75.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.75)
max75.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.75)
max75.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.75)
max75.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.75)


maxvpbtw = data.frame(rbind(rbind(max90.rodbtw$prop_sp, max90.bfbtw$prop_sp, max90.mothbtw$prop_sp, max90.bdbtw$prop_sp, max90.wgbtw$prop_sp),
                         rbind(max95.rodbtw$prop_sp, max95.bfbtw$prop_sp, max95.mothbtw$prop_sp, max95.bdbtw$prop_sp, max95.wgbtw$prop_sp),
                         rbind(max80.rodbtw$prop_sp, max80.bfbtw$prop_sp, max80.mothbtw$prop_sp, max80.bdbtw$prop_sp, max80.wgbtw$prop_sp),
                         rbind(max75.rodbtw$prop_sp, max75.bfbtw$prop_sp, max75.mothbtw$prop_sp, max75.bdbtw$prop_sp, max75.wgbtw$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                   cond = 'btw')
colnames(maxvpbtw)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')


maxvp = rbind(maxvpNB, maxvpMB, maxvpSB, maxvpbtw, maxvpwit)


maxvplg <- maxvp %>% 
  pivot_longer(
    cols = 'Climate':'HabComp', 
    names_to = "Driver",
    values_to = "prop"
  )


maxvplg$Driver = factor(maxvplg$Driver, levels = c('Climate', 'LandscapeConf', 'HabComp'))
maxvplg$taxa = factor(maxvplg$taxa,levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
maxvplg$cond = factor(maxvplg$cond, levels = c('NB', 'MB', 'SB', 'btw', 'wit'))
levels(maxvplg$cond) = c('North Boreal', 'Middle Boreal', 'South Boreal', 'Between group', 'Within group')

minors <- c(0.25, 0.5, 0.75, 1)


ggocc = ggplot(data = maxvplg[maxvplg$Prob_post == '90%',], aes(x=prop, y = Driver, group = Driver, fill = Driver)) +
  #geom_bar(stat = "identity") +
  geom_vline(mapping=NULL, xintercept=minors,colour='grey90') +
  geom_pointrange(aes(y = Driver, 
                      xmin = maxvplg[maxvplg$Prob_post == '95%',]$prop, 
                      xmax = maxvplg[maxvplg$Prob_post == '75%',]$prop,
                      color = Driver), 
                  linewidth = 2, size = 2) +
  scale_fill_manual(values = alpha(c("#66C2A5", "yellow3", 'orange3'), 0.5)) +
  scale_colour_manual(values = c("#66C2A5", "yellow3", 'orange3')) +
  facet_grid(cond~ taxa)  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab('Proportion of species') +
  coord_flip() +  theme_classic()  + 
  theme(axis.title.x = element_text(vjust=-1),
        axis.text.x = element_text(angle=45, vjust=1),
  ) + easy_remove_x_axis(c("ticks", "title", "line"))







#............................................... Abundance


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

VPfix.rodNB = varr.rod*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.rodMB = varr.rod*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.rodSB = varr.rod*VP.test_cond[1:3, 1:3,,,'SB'] # SB


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



VPfix.rodbtw = varr.rod*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.rodwit = varr.rod*VP.test_bw[1:3, 1:3,,,'wit'] # 


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

VPfix.bfNB = varr.bf*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.bfMB = varr.bf*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.bfSB = varr.bf*VP.test_cond[1:3, 1:3,,,'SB'] # SB

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



VPfix.bfbtw = varr.bf*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.bfwit = varr.bf*VP.test_bw[1:3, 1:3,,,'wit'] # 


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

VPfix.mothNB = varr.moth*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.mothMB = varr.moth*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.mothsB = varr.moth*VP.test_cond[1:3, 1:3,,,'SB'] # SB

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



VPfix.mothbtw = varr.moth*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.mothwit = varr.moth*VP.test_bw[1:3, 1:3,,,'wit'] # 



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

VPfix.wgNB = varr.wg*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.wgMB = varr.wg*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.wgSB = varr.wg*VP.test_cond[1:3, 1:3,,,'SB'] # SB

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



VPfix.wgbtw = varr.wg*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.wgwit = varr.wg*VP.test_bw[1:3, 1:3,,,'wit'] # 


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

VPfix.bdNB = varr.bd*VP.test_cond[1:3, 1:3,,,'NB'] # NB
VPfix.bdMB = varr.bd*VP.test_cond[1:3, 1:3,,,'MB'] # MB
VPfix.bdSB = varr.bd*VP.test_cond[1:3, 1:3,,,'SB'] # SB


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



VPfix.bdbtw = varr.bd*VP.test_bw[1:3, 1:3,,,'btw'] # 
VPfix.bdwit = varr.bd*VP.test_bw[1:3, 1:3,,,'wit'] # 



##.............................  Summarize


spvp_summaries = function(vpdat, Threshold){
  
  # Initialize a matrix to store the count of being the maximum
  # rows are the fixed effect only but columns are fixed effect and random effects we need to focus only on fixed egffect
  max_count <- array(0, dim = c(dim(vpdat)[1], dim(vpdat)[3]))
  #dim(max_count)  # fix effect # species # posterior sample
  
  # For each species and posterior sample, find the covariate with the maximum value
  for (s in 1:dim(vpdat)[3]) {
    for (p in 1:dim(vpdat)[4]) {
      # Find the index of the maximum value for this species and sample
      max_val <- max(diag(vpdat[,,s,p]))
      max_indices <- which(diag(vpdat[,,s,p]) == max_val, arr = TRUE)
      
      max_count[max_indices, s] <- max_count[max_indices, s] + 1
    }
  }
  
  # Calculate the proportion of posterior samples where each covariate is the maximum for each species
  # calculate among the species that are not NA (R2) so using ColSums instead of dim(Vpdat[4])
  prop_max <- max_count / colSums(max_count)
  
  
  prop_sp <- rowMeans(prop_max >= Threshold)
  
  # prop_species now contains the proportion for each covariate
  print(prop_sp)
  
  return(list(prop_sp = prop_sp))
}

# different posterior probabilities

## NB
# 90%
max90.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.90)
max90.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.90)
max90.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.90)
max90.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.90)
max90.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.90)


# 95%
max95.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.95)
max95.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.95)
max95.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.95)
max95.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.95)
max95.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.95)


# 80%
max80.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.80)
max80.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.80)
max80.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.80)
max80.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.80)
max80.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.80)

# 75%
max75.rodNB = spvp_summaries(VPfix.rodNB, Threshold = 0.75)
max75.bfNB = spvp_summaries(VPfix.bfNB, Threshold = 0.75)
max75.mothNB = spvp_summaries(VPfix.mothNB, Threshold = 0.75)
max75.wgNB = spvp_summaries(VPfix.wgNB, Threshold = 0.75)
max75.bdNB = spvp_summaries(VPfix.bdNB, Threshold = 0.75)


maxvpNB = data.frame(rbind(rbind(max90.rodNB$prop_sp, max90.bfNB$prop_sp, max90.mothNB$prop_sp, max90.bdNB$prop_sp, max90.wgNB$prop_sp),
                           rbind(max95.rodNB$prop_sp, max95.bfNB$prop_sp, max95.mothNB$prop_sp, max95.bdNB$prop_sp, max95.wgNB$prop_sp),
                           rbind(max80.rodNB$prop_sp, max80.bfNB$prop_sp, max80.mothNB$prop_sp, max80.bdNB$prop_sp, max80.wgNB$prop_sp),
                           rbind(max75.rodNB$prop_sp, max75.bfNB$prop_sp, max75.mothNB$prop_sp, max75.bdNB$prop_sp, max75.wgNB$prop_sp)),
                     Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                     taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                     cond = 'NB')
colnames(maxvpNB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## MB
# 90%
max90.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.90)
max90.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.90)
max90.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.90)
max90.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.90)
max90.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.90)


# 95%
max95.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.95)
max95.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.95)
max95.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.95)
max95.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.95)
max95.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.95)


# 80%
max80.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.80)
max80.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.80)
max80.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.80)
max80.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.80)
max80.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.80)

# 75%
max75.rodMB = spvp_summaries(VPfix.rodMB, Threshold = 0.75)
max75.bfMB = spvp_summaries(VPfix.bfMB, Threshold = 0.75)
max75.mothMB = spvp_summaries(VPfix.mothMB, Threshold = 0.75)
max75.wgMB = spvp_summaries(VPfix.wgMB, Threshold = 0.75)
max75.bdMB = spvp_summaries(VPfix.bdMB, Threshold = 0.75)


maxvpMB = data.frame(rbind(rbind(max90.rodMB$prop_sp, max90.bfMB$prop_sp, max90.mothMB$prop_sp, max90.bdMB$prop_sp, max90.wgMB$prop_sp),
                           rbind(max95.rodMB$prop_sp, max95.bfMB$prop_sp, max95.mothMB$prop_sp, max95.bdMB$prop_sp, max95.wgMB$prop_sp),
                           rbind(max80.rodMB$prop_sp, max80.bfMB$prop_sp, max80.mothMB$prop_sp, max80.bdMB$prop_sp, max80.wgMB$prop_sp),
                           rbind(max75.rodMB$prop_sp, max75.bfMB$prop_sp, max75.mothMB$prop_sp, max75.bdMB$prop_sp, max75.wgMB$prop_sp)),
                     Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                     taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                     cond = 'MB')
colnames(maxvpMB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## SB
# 90%
max90.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.90)
max90.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.90)
max90.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.90)
max90.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.90)
max90.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.90)


# 95%
max95.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.95)
max95.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.95)
max95.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.95)
max95.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.95)
max95.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.95)


# 80%
max80.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.80)
max80.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.80)
max80.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.80)
max80.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.80)
max80.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.80)

# 75%
max75.rodSB = spvp_summaries(VPfix.rodSB, Threshold = 0.75)
max75.bfSB = spvp_summaries(VPfix.bfSB, Threshold = 0.75)
max75.mothsB = spvp_summaries(VPfix.mothsB, Threshold = 0.75)
max75.wgSB = spvp_summaries(VPfix.wgSB, Threshold = 0.75)
max75.bdSB = spvp_summaries(VPfix.bdSB, Threshold = 0.75)


maxvpSB = data.frame(rbind(rbind(max90.rodSB$prop_sp, max90.bfSB$prop_sp, max90.mothsB$prop_sp, max90.bdSB$prop_sp, max90.wgSB$prop_sp),
                           rbind(max95.rodSB$prop_sp, max95.bfSB$prop_sp, max95.mothsB$prop_sp, max95.bdSB$prop_sp, max95.wgSB$prop_sp),
                           rbind(max80.rodSB$prop_sp, max80.bfSB$prop_sp, max80.mothsB$prop_sp, max80.bdSB$prop_sp, max80.wgSB$prop_sp),
                           rbind(max75.rodSB$prop_sp, max75.bfSB$prop_sp, max75.mothsB$prop_sp, max75.bdSB$prop_sp, max75.wgSB$prop_sp)),
                     Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                     taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                     cond = 'SB')
colnames(maxvpSB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## wit
# 90%
max90.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.90)
max90.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.90)
max90.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.90)
max90.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.90)
max90.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.90)


# 95%
max95.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.95)
max95.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.95)
max95.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.95)
max95.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.95)
max95.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.95)


# 80%
max80.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.80)
max80.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.80)
max80.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.80)
max80.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.80)
max80.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.80)

# 75%
max75.rodwit = spvp_summaries(VPfix.rodwit, Threshold = 0.75)
max75.bfwit = spvp_summaries(VPfix.bfwit, Threshold = 0.75)
max75.mothwit = spvp_summaries(VPfix.mothwit, Threshold = 0.75)
max75.wgwit = spvp_summaries(VPfix.wgwit, Threshold = 0.75)
max75.bdwit = spvp_summaries(VPfix.bdwit, Threshold = 0.75)


maxvpwit = data.frame(rbind(rbind(max90.rodwit$prop_sp, max90.bfwit$prop_sp, max90.mothwit$prop_sp, max90.bdwit$prop_sp, max90.wgwit$prop_sp),
                            rbind(max95.rodwit$prop_sp, max95.bfwit$prop_sp, max95.mothwit$prop_sp, max95.bdwit$prop_sp, max95.wgwit$prop_sp),
                            rbind(max80.rodwit$prop_sp, max80.bfwit$prop_sp, max80.mothwit$prop_sp, max80.bdwit$prop_sp, max80.wgwit$prop_sp),
                            rbind(max75.rodwit$prop_sp, max75.bfwit$prop_sp, max75.mothwit$prop_sp, max75.bdwit$prop_sp, max75.wgwit$prop_sp)),
                      Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                      taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                      cond = 'wit')
colnames(maxvpwit)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')



## btw
# 90%
max90.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.90)
max90.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.90)
max90.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.90)
max90.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.90)
max90.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.90)


# 95%
max95.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.95)
max95.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.95)
max95.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.95)
max95.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.95)
max95.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.95)


# 80%
max80.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.80)
max80.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.80)
max80.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.80)
max80.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.80)
max80.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.80)

# 75%
max75.rodbtw = spvp_summaries(VPfix.rodbtw, Threshold = 0.75)
max75.bfbtw = spvp_summaries(VPfix.bfbtw, Threshold = 0.75)
max75.mothbtw = spvp_summaries(VPfix.mothbtw, Threshold = 0.75)
max75.wgbtw = spvp_summaries(VPfix.wgbtw, Threshold = 0.75)
max75.bdbtw = spvp_summaries(VPfix.bdbtw, Threshold = 0.75)


maxvpbtw = data.frame(rbind(rbind(max90.rodbtw$prop_sp, max90.bfbtw$prop_sp, max90.mothbtw$prop_sp, max90.bdbtw$prop_sp, max90.wgbtw$prop_sp),
                            rbind(max95.rodbtw$prop_sp, max95.bfbtw$prop_sp, max95.mothbtw$prop_sp, max95.bdbtw$prop_sp, max95.wgbtw$prop_sp),
                            rbind(max80.rodbtw$prop_sp, max80.bfbtw$prop_sp, max80.mothbtw$prop_sp, max80.bdbtw$prop_sp, max80.wgbtw$prop_sp),
                            rbind(max75.rodbtw$prop_sp, max75.bfbtw$prop_sp, max75.mothbtw$prop_sp, max75.bdbtw$prop_sp, max75.wgbtw$prop_sp)),
                      Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                      taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4),
                      cond = 'btw')
colnames(maxvpbtw)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')


maxvp = rbind(maxvpNB, maxvpMB, maxvpSB, maxvpbtw, maxvpwit)


maxvplg <- maxvp %>% 
  pivot_longer(
    cols = 'Climate':'HabComp', 
    names_to = "Driver",
    values_to = "prop"
  )


maxvplg$Driver = factor(maxvplg$Driver, levels = c('Climate', 'LandscapeConf', 'HabComp'))
maxvplg$taxa = factor(maxvplg$taxa,levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
maxvplg$cond = factor(maxvplg$cond, levels = c('NB', 'MB', 'SB', 'btw', 'wit'))
levels(maxvplg$cond) = c('North Boreal', 'Middle Boreal', 'South Boreal', 'Between group', 'Within group')

minors <- c(0.25, 0.5, 0.75, 1)


ggAB = ggplot(data = maxvplg[maxvplg$Prob_post == '90%',], aes(x=prop, y = Driver, group = Driver, fill = Driver)) +
  #geom_bar(stat = "identity") +
  geom_vline(mapping=NULL, xintercept=minors,colour='grey90') +
  geom_pointrange(aes(y = Driver, 
                      xmin = maxvplg[maxvplg$Prob_post == '95%',]$prop, 
                      xmax = maxvplg[maxvplg$Prob_post == '75%',]$prop,
                      color = Driver), 
                  linewidth = 2, size = 2) +
  scale_fill_manual(values = alpha(c("#66C2A5", "yellow3", 'orange3'), 0.5)) +
  scale_colour_manual(values = c("#66C2A5", "yellow3", 'orange3')) +
  facet_grid(cond~ taxa)  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab('Proportion of species') +
  coord_flip() +  theme_classic()  + 
  theme(axis.title.x = element_text(vjust=-1),
        axis.text.x = element_text(angle=45, vjust=1),
  ) + easy_remove_x_axis(c("ticks", "title", "line"))

