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

setwd("~/Data/4.VPanalyses/VPCP") # to adapt depending on where it is saved

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


varr <- aperm(array(TabMF_rod$TjurR2, dim = c(8L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

VPfixR2.rod = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_bf$TjurR2, dim = c(57L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)


VPfixR2.bf = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_wg$TjurR2, dim = c(15L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)


VPfixR2.wg = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_bd$TjurR2, dim = c(102L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))

VPfixR2.bd = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_Moth$TjurR2, dim = c(319L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)


VPfixR2.moth = varr*VP.test3[1:3, 1:3, ,]




##.............................  Summarize


dim(VPfixR2.bd)  #  fixed eff # fixed eff  # species # posterior sample
dim(VPfixR2.moth)
dim(VPfixR2.bf)
dim(VPfixR2.rod)
dim(VPfixR2.wg)


# test
vpdat = VPfixR2.rod

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
# 90%
max90.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.90)
max90.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.90)
max90.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.90)
max90.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.90)
max90.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.90)


# 95%
max95.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.95)
max95.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.95)
max95.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.95)
max95.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.95)
max95.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.95)


# 80%
max80.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.80)
max80.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.80)
max80.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.80)
max80.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.80)
max80.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.80)

# 75%
max75.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.75)
max75.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.75)
max75.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.75)
max75.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.75)
max75.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.75)


maxvp = data.frame(rbind(rbind(max90.rod$prop_sp, max90.bf$prop_sp, max90.moth$prop_sp, max90.bd$prop_sp, max90.wg$prop_sp),
                         rbind(max95.rod$prop_sp, max95.bf$prop_sp, max95.moth$prop_sp, max95.bd$prop_sp, max95.wg$prop_sp),
                         rbind(max80.rod$prop_sp, max80.bf$prop_sp, max80.moth$prop_sp, max80.bd$prop_sp, max80.wg$prop_sp),
                         rbind(max75.rod$prop_sp, max75.bf$prop_sp, max75.moth$prop_sp, max75.bd$prop_sp, max75.wg$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4))
colnames(maxvp)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')


maxvplg <- maxvp %>% 
  pivot_longer(
    cols = 'Climate':'HabComp', 
    names_to = "Driver",
    values_to = "prop"
  )


maxvplg$Driver = factor(maxvplg$Driver, levels = c('Climate', 'LandscapeConf', 'HabComp'))
maxvplg$taxa = factor(maxvplg$taxa,levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
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
  facet_wrap(.~ taxa, nrow=1)  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab('Proportion of species') +
  coord_flip() +  theme_classic()  + 
  theme(axis.title.x = element_text(vjust=-1),
        axis.text.x = element_text(angle=45, vjust=1),
  ) + easy_remove_x_axis(c("ticks", "title", "line"))


##...............................................  Abundance data


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


varr <- aperm(array(TabMF_rod$SR2, dim = c(8L, 3L, 3L, 2000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

VPfixR2.rod = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_bf$SR2, dim = c(57L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)


VPfixR2.bf = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_wg$SR2, dim = c(15L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)


VPfixR2.wg = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_bd$SR2, dim = c(102L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))

VPfixR2.bd = varr*VP.test3[1:3, 1:3, ,]



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


varr <- aperm(array(TabMF_Moth$SR2, dim = c(319L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)


VPfixR2.moth = varr*VP.test3[1:3, 1:3, ,]




##.............................  Summarize


dim(VPfixR2.bd)  #  fixed eff # fixed eff + rd  # species # posterior sample
dim(VPfixR2.moth)
dim(VPfixR2.bf)
dim(VPfixR2.rod)
dim(VPfixR2.wg)


# test
# vpdat = VPfixR2.rod
# vpdat = VPfixR2.moth


spvp_summaries = function(vpdat, Threshold){

  # Initialize a matrix to store the count of being the maximum
  # rows are the fixed effect only but columns are fixed effect and random effects we need to focus only on fixed egffect
  max_count <- array(0, dim = c(dim(vpdat)[1], dim(vpdat)[3]))
  #dim(max_count)  # fix effect # species # posterior sample
  
  # For each species and posterior sample, find the covariate with the maximum value
  for (s in 1:dim(vpdat)[3]) {
    for (p in 1:dim(vpdat)[4]) {
      # Find the index of the maximum value for this species and sample
      # 
      is.na(vpdat[,,s,p]) = 0 # prevent numerical problens
      max_val <- max(diag(vpdat[,,s,p]), na.rm=T)
      max_indices <- which(diag(vpdat[,,s,p]) == max_val, arr = TRUE)
      
      max_count[max_indices, s] <- max_count[max_indices, s] + 1
    }
  }
  
  # Calculate the proportion of posterior samples where each covariate is the maximum for each species
  # calculate among the species that are not NA (R2) so using ColSums instead of dim(Vpdat[4])
  prop_max <- max_count / colSums(max_count)
  
  prop_sp <- rowMeans(prop_max >= Threshold, na.rm=T)
  
  # prop_species now contains the proportion for each covariate
  print(prop_sp)
  
  return(list(prop_sp = prop_sp))
}

# different posterior probabilities
# 90%
max90.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.90)
max90.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.90)
max90.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.90)
max90.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.90)
max90.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.90)


# 95%
max95.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.95)
max95.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.95)
max95.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.95)
max95.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.95)
max95.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.95)


# 80%
max80.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.80)
max80.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.80)
max80.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.80)
max80.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.80)
max80.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.80)

# 75%
max75.rod = spvp_summaries(VPfixR2.rod, Threshold = 0.75)
max75.bf = spvp_summaries(VPfixR2.bf, Threshold = 0.75)
max75.moth = spvp_summaries(VPfixR2.moth, Threshold = 0.75)
max75.wg = spvp_summaries(VPfixR2.wg, Threshold = 0.75)
max75.bd = spvp_summaries(VPfixR2.bd, Threshold = 0.75)


maxvpAB = data.frame(rbind(rbind(max90.rod$prop_sp, max90.bf$prop_sp, max90.moth$prop_sp, max90.bd$prop_sp, max90.wg$prop_sp),
                         rbind(max95.rod$prop_sp, max95.bf$prop_sp, max95.moth$prop_sp, max95.bd$prop_sp, max95.wg$prop_sp),
                         rbind(max80.rod$prop_sp, max80.bf$prop_sp, max80.moth$prop_sp, max80.bd$prop_sp, max80.wg$prop_sp),
                         rbind(max75.rod$prop_sp, max75.bf$prop_sp, max75.moth$prop_sp, max75.bd$prop_sp, max75.wg$prop_sp)),
                   Prob_post = c(rep('90%', 5), rep('95%', 5), rep('80%', 5), rep('75%', 5)),
                   taxa = rep(c('small mammals', 'butterflies', 'moths', 'birds', 'large mammals'), 4))
colnames(maxvpAB)[1:3] =  c('Climate', 'LandscapeConf', 'HabComp')


maxvplgAB <- maxvpAB %>% 
  pivot_longer(
    cols = 'Climate':'HabComp', 
    names_to = "Driver",
    values_to = "prop"
  )


maxvplgAB$Driver = factor(maxvplgAB$Driver, levels = c('Climate', 'LandscapeConf', 'HabComp'))
maxvplgAB$taxa = factor(maxvplgAB$taxa,levels = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals'))
minors <- c(0.25, 0.5, 0.75, 1)

ggAB = ggplot(data = maxvplgAB[maxvplgAB$Prob_post == '90%',], aes(x=prop, y = Driver, group = Driver, fill = Driver)) +
  #geom_bar(stat = "identity") +
  geom_vline(mapping=NULL, xintercept=minors,colour='grey90') +
  geom_pointrange(aes(y = Driver, 
                      xmin = maxvplgAB[maxvplgAB$Prob_post == '95%',]$prop, 
                      xmax = maxvplgAB[maxvplgAB$Prob_post == '75%',]$prop,
                      color = Driver), 
                  linewidth = 2, size = 2) +
  scale_fill_manual(values = alpha(c("#66C2A5", "yellow3", 'orange3'), 0.5)) +
  scale_colour_manual(values = c("#66C2A5", "yellow3", 'orange3')) +
  facet_wrap(.~ taxa, nrow=1)  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab('Proportion of species') +
  coord_flip() +  theme_classic() +
  theme(axis.title.x = element_text(vjust=-1),
        axis.text.x = element_text(angle=45, vjust=1),
  ) + easy_remove_x_axis(c("ticks", "title", "line")) 


## combine

ggarrange(ggocc, ggAB, nrow=2)
