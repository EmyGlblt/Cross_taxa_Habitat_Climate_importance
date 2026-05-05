#
#
#            VP trait full model & all taxa extract
#


## libraries
library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)

# #
##.......................................................................................................
##  Rodents
##.......................................................................................................
# #
# #    Occurrence


setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_rod.RDATA') # to get variance and covariance idf interested in. (VP.test3)

dim(VP.test3_rod$Vnorm)
dim(VP.test3_rod$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_rod$Cnorm))

for (i in 1:dim(VP.test3_rod$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_rod$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_rod$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_rod$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_rod$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_rod$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_rod$TjurR2, dim = c(8L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keepin the Normalized variance
# Posterior mean
Mean.fix.sprod = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)

# full distribution
VP.fix.ROD = varr*VP.test3[1:3, 1:3, ,]

rod.sp = TabMF_rod$species

setwd("~/Data/2.Prep_data/Traits")
load('Traits_rodmamal.RDATA')

dim(Rodent_trait)
colnames(Rodent_trait)
Rodent_trait_sub = Rodent_trait[, c('iucn2020_binomial', 'adult_mass_g',
                                    'adult_body_length_mm', 'litter_size_n',
                                    'litters_per_year_n', 'habitat_breadth_n',
                                    'det_diet_breadth_n', "dispersal_km")]

Rodent_trait_sub = Rodent_trait_sub[Rodent_trait_sub$iucn2020_binomial %in% rod.sp,]
dim(Rodent_trait_sub)

# combine data
VPclim = 100*Mean.fix.sprod[1,1,]
VPhabC = 100*Mean.fix.sprod[2,2,]
VPhabP = 100*Mean.fix.sprod[3,3,]

rodVP = data.frame(VPClim = VPclim,
                   VPhabC = VPhabC,
                   VPhabP = VPhabP)

dim(Rodent_trait_sub)
dim(rodVP)
rodVP$iucn2020_binomial = rod.sp ## order was kept when calculating the VP

rodDat_trait = merge(rodVP, Rodent_trait_sub, by= 'iucn2020_binomial')
dim(rodDat_trait)

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(rodDat_trait, file='VPtrait_rod.RDATA')


# #
# #    Abundance

setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_ab_rod.RDATA')
dim(VP.test3_rod$Vnorm)
dim(VP.test3_rod$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_rod$Cnorm))

for (i in 1:dim(VP.test3_rod$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_rod$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_rod$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_rod$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_rod$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_rod$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_rod$SR2, dim = c(8L, 3L, 3L, 2000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keeping the Normalized variance
# Posterior mean
Mean.fix.sprod = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)
VPclim = 100*Mean.fix.sprod[1,1,]
VPhabC = 100*Mean.fix.sprod[2,2,]
VPhabP = 100*Mean.fix.sprod[3,3,]

rodVP_ab = data.frame(VPClim = VPclim,
                   VPhabC = VPhabC,
                   VPhabP = VPhabP)

dim(Rodent_trait_sub)
dim(rodVP_ab)
rodVP_ab$iucn2020_binomial = rod.sp


rodDat_trait_ab = merge(rodVP_ab, Rodent_trait_sub, by= 'iucn2020_binomial')
dim(rodDat_trait_ab)

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(rodDat_trait_ab, file='VPtrait_ab_rod.RDATA')



##.......................................................................................................
##  Butterflies
##.......................................................................................................

# #
# #    Occurrence

setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_bf.RDATA') # to get variance and covariance idf interested in. (VP.test3)

dim(VP.test3_bf$Vnorm)
dim(VP.test3_bf$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bf$Cnorm))

for (i in 1:dim(VP.test3_bf$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bf$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bf$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bf$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bf$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bf$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_bf$TjurR2, dim = c(57L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keepin the Normalized variance
# Posterior mean
Mean.fix.spbf = apply(varr*VP.test3[1:3, 1:3, ,], c(1,2,3), mean, na.rm=T)

# full distribution
VP.fix.bf = varr*VP.test3[1:3, 1:3, ,]

bf.sp = TabMF_bf$species


setwd("~/Data/2.Prep_data/Traits")
bf_trait = read.csv('Lep_traits.csv', header=T, sep=',')

dim(bf_trait)
colnames(bf_trait)
bf_trait_sub = bf_trait[, c('names_rec', 'Body.size', 
                            'Specificity.of.larval.host.plant.use', 
                            'Breadth.of.habitat.use.adult.butterflies',
                            'Voltinism')]

bf_trait_sub = bf_trait_sub[bf_trait_sub$Specificity.of.larval.host.plant.use %in% c('monophagous', 'oligophagous', 'polyphagous'),]

VPclim = 100*Mean.fix.spbf[1,1,]
VPhabC = 100*Mean.fix.spbf[2,2,]
VPhabP = 100*Mean.fix.spbf[3,3,]

bfVP = data.frame(VPClim = VPclim,
                  VPhabC = VPhabC,
                  VPhabP = VPhabP)

dim(bf_trait_sub)
dim(bfVP)
bfVP$names_rec = bf.sp


BfDat_trait = merge(bfVP, bf_trait_sub, by= 'names_rec')

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(BfDat_trait, file='VPtrait_bf.RDATA')

# #
# #    Abundance

setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_ab_bf.RDATA')
dim(VP.test3_bf$Vnorm)
dim(VP.test3_bf$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bf$Cnorm))

for (i in 1:dim(VP.test3_bf$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bf$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bf$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bf$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bf$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bf$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_bf$SR2, dim = c(57L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keeping the Normalized variance
# Posterior mean
Mean.fix.spbf = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)
VPclim = 100*Mean.fix.spbf[1,1,]
VPhabC = 100*Mean.fix.spbf[2,2,]
VPhabP = 100*Mean.fix.spbf[3,3,]

bfVP_ab = data.frame(VPClim = VPclim,
                  VPhabC = VPhabC,
                  VPhabP = VPhabP)

dim(bf_trait_sub)
dim(bfVP_ab)
bfVP_ab$names_rec = bf.sp


BfDat_trait_ab = merge(bfVP_ab, bf_trait_sub, by= 'names_rec')

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(BfDat_trait_ab, file='VPtrait_ab_bf.RDATA')


##.......................................................................................................
##  Moths
##.......................................................................................................
# #
# #    Occurrence


setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_moth.RDATA') # to get variance and covariance idf interested in. (VP.test3)

dim(VP.test3_moth$Vnorm)
dim(VP.test3_moth$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_moth$Cnorm))

for (i in 1:dim(VP.test3_moth$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_moth$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_moth$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_moth$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_moth$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_moth$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_Moth$TjurR2, dim = c(319L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keepin the Normalized variance
# Posterior mean
Mean.fix.spmoth = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)

# full distribution
VP.fix.moth = varr*VP.test3[1:3, 1:3, ,]

moth.sp = TabMF_Moth$species

setwd("~/Data/2.Prep_data/Traits")
moth_trait = read.csv('Traits_Moths.csv', header=T, sep=',')

dim(moth_trait)
colnames(moth_trait)
head(moth_trait)
moth_trait_sub = moth_trait[, c('Genus.and.species', 'female.size..average.wingspan..mm.', 
                                'Diet', 
                                'HostPlantForm',
                                'Voltinism')]

colnames(moth_trait_sub)[2] = 'Body.Size'
moth_trait_sub$Body.Size = as.numeric(moth_trait_sub$Body.Size)
moth.sp[!moth.sp %in% moth_trait_sub$Genus.and.species]


moth_trait_sub = moth_trait_sub[moth_trait_sub$Genus.and.species %in% moth.sp,]
summary(moth_trait_sub)

## combine
VPclim = 100*Mean.fix.spmoth[1,1,]
VPhabC = 100*Mean.fix.spmoth[2,2,]
VPhabP = 100*Mean.fix.spmoth[3,3,]

mothVP = data.frame(VPClim = VPclim,
                    VPhabC = VPhabC,
                    VPhabP = VPhabP)

dim(moth_trait_sub)
dim(mothVP)
mothVP$Genus.and.species = moth.sp


mothDat_trait = merge(mothVP, moth_trait_sub, by= 'Genus.and.species')

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(mothDat_trait, file='VPtrait_moth.RDATA')


# #
# #    Abundance
setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_ab_moth.RDATA')
dim(VP.test3_moth$Vnorm)
dim(VP.test3_moth$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_moth$Cnorm))

for (i in 1:dim(VP.test3_moth$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_moth$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_moth$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_moth$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_moth$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_moth$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_Moth$SR2, dim = c(319L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keeping the Normalized variance
# Posterior mean
Mean.fix.spmoth = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)
VPclim = 100*Mean.fix.spmoth[1,1,]
VPhabC = 100*Mean.fix.spmoth[2,2,]
VPhabP = 100*Mean.fix.spmoth[3,3,]

mothVP_ab = data.frame(VPClim = VPclim,
                    VPhabC = VPhabC,
                    VPhabP = VPhabP)

dim(moth_trait_sub)
dim(mothVP_ab)
mothVP_ab$Genus.and.species = moth.sp


mothDat_trait_ab = merge(mothVP_ab, moth_trait_sub, by= 'Genus.and.species')

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(mothDat_trait_ab, file='VPtrait_ab_moth.RDATA')


##.......................................................................................................
##  Birds
##.......................................................................................................
# #
# #    Occurrence


setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_bd.RDATA') # to get variance and covariance idf interested in. (VP.test3)

dim(VP.test3_bd$Vnorm)
dim(VP.test3_bd$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bd$Cnorm))

for (i in 1:dim(VP.test3_bd$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bd$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bd$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bd$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bd$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bd$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_bd$TjurR2, dim = c(102L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keepin the Normalized variance
# Posterior mean
Mean.fix.spbd = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)

# full distribution
VP.fix.bd = varr*VP.test3[1:3, 1:3, ,]

bd.sp = TabMF_bd$species

setwd("~/Data/2.Prep_data/Traits")
bird_trait = read.csv('BirdTraits21112018.csv', header=T, sep=',')

dim(bird_trait)
colnames(bird_trait)

## STILL MISSING SOME NAMES COMPUTE THE MATCHING
Bname_match = read.csv('Birds_names_match.csv', header=T, sep=',')
head(Bname_match)

# direction old to new
colnames(Bname_match)[2] = 'Sp_Lat'

bird_trait2 = merge(bird_trait, Bname_match, by = 'Sp_Lat', all.x = T)

# direction new to old
colnames(Bname_match)[1:2] = c('Sp_Lat', 'old.names')
bird_trait3 = merge(bird_trait2, Bname_match, by = 'Sp_Lat', all.x = T)

# Other birds traits:
bird_spe = read.csv('Bird_specialism.csv', header=T, sep=',')
head(bird_spe)

length(bird_trait3$Sp_Lat %in% bird_spe$Species) # 245

bird_trait_sub = bird_trait3[(bird_trait3$Sp_Lat %in% bd.sp | 
                                bird_trait3$New.Species.names %in% bd.sp |
                                bird_trait3$old.names %in% bd.sp) ,]
dim(bird_trait_sub)
table(bird_trait_sub$Sp_Lat %in% bird_spe$Species) # 102


which.arenot = bird_trait_sub$Sp_Lat %in% bird_spe$Species
bird_trait_sub$Sp_Lat[!which.arenot]


colnames(bird_spe)[3] = 'Sp_Lat'
bird_trait_sub2 = merge(bird_trait_sub, bird_spe, by = 'Sp_Lat')
dim(bird_trait_sub2)  # 

# Combine
VPclim = 100*Mean.fix.spbd[1,1,]
VPhabC = 100*Mean.fix.spbd[2,2,]
VPhabP = 100*Mean.fix.spbd[3,3,]

bdVP = data.frame(VPClim = VPclim,
                  VPhabC = VPhabC,
                  VPhabP = VPhabP)

dim(bird_trait_sub2)

bird_trait_sub2$sp_name = bird_trait_sub2$Sp_Lat
bird_trait_sub2$sp_name[!is.na(bird_trait_sub2$old.names)] = bird_trait_sub2$old.names[!is.na(bird_trait_sub2$old.names)]
bird_trait_sub2$sp_name[!is.na(bird_trait_sub2$New.Species.names)] = bird_trait_sub2$New.Species.names[!is.na(bird_trait_sub2$New.Species.names)]

dim(bdVP)
bdVP$sp_name = bd.sp

BdDat_trait = merge(bdVP, bird_trait_sub2, by= 'sp_name')
dim(BdDat_trait)

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(BdDat_trait, file='VPtrait_bird2.RDATA')


# #
# #    Abundance
setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_ab_bd.RDATA')
dim(VP.test3_bd$Vnorm)
dim(VP.test3_bd$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_bd$Cnorm))

for (i in 1:dim(VP.test3_bd$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_bd$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_bd$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_bd$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_bd$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_bd$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_bd$SR2, dim = c(102L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keeping the Normalized variance
# Posterior mean
Mean.fix.spbd = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)
VPclim = 100*Mean.fix.spbd[1,1,]
VPhabC = 100*Mean.fix.spbd[2,2,]
VPhabP = 100*Mean.fix.spbd[3,3,]

bdVP_ab = data.frame(VPClim = VPclim,
                  VPhabC = VPhabC,
                  VPhabP = VPhabP)

dim(bird_trait_sub2)

bird_trait_sub2$sp_name = bird_trait_sub2$Sp_Lat
bird_trait_sub2$sp_name[!is.na(bird_trait_sub2$old.names)] = bird_trait_sub2$old.names[!is.na(bird_trait_sub2$old.names)]
bird_trait_sub2$sp_name[!is.na(bird_trait_sub2$New.Species.names)] = bird_trait_sub2$New.Species.names[!is.na(bird_trait_sub2$New.Species.names)]

dim(bdVP_ab)
bdVP_ab$sp_name = bd.sp

BdDat_trait_ab = merge(bdVP_ab, bird_trait_sub2, by= 'sp_name')
dim(BdDat_trait_ab)


setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(BdDat_trait_ab, file='VPtrait_ab_bird.RDATA')



##.......................................................................................................
##  Game triangle
##.......................................................................................................
# #
# #    Occurrence
setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_wg.RDATA') # to get variance and covariance idf interested in. (VP.test3)

dim(VP.test3_wg$Vnorm)
dim(VP.test3_wg$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_wg$Cnorm))

for (i in 1:dim(VP.test3_wg$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_wg$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_wg$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_wg$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_wg$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_wg$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_wg$TjurR2, dim = c(15L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keepin the Normalized variance
# Posterior mean
Mean.fix.spwg = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)

# full distribution
VP.fix.wg = varr*VP.test3[1:3, 1:3, ,]

wg.sp = TabMF_wg$species


setwd("~/Data/2.Prep_data/Traits")
load('Traits_rodmamal.RDATA')

dim(Mammals_trait)
colnames(Mammals_trait)
wintG_trait_sub = Mammals_trait[, c('iucn2020_binomial', 'adult_mass_g',
                                    'adult_body_length_mm', 'litter_size_n',
                                    'litters_per_year_n', 'habitat_breadth_n',
                                    'det_diet_breadth_n', "dispersal_km")]

wintG_trait_sub = wintG_trait_sub[wintG_trait_sub$iucn2020_binomial %in% wg.sp,]
dim(wintG_trait_sub)

# combine
VPclim = 100*Mean.fix.spwg[1,1,]
VPhabC = 100*Mean.fix.spwg[2,2,]
VPhabP = 100*Mean.fix.spwg[3,3,]

wgVP = data.frame(VPClim = VPclim,
                  VPhabC = VPhabC,
                  VPhabP = VPhabP)

dim(wintG_trait_sub)
dim(wgVP)
wgVP$iucn2020_binomial = wg.sp

wgDat_trait = merge(wgVP, wintG_trait_sub, by= 'iucn2020_binomial')
setwd("~/Data/4.VPanalyses/VP_cond_trait")

save(wgDat_trait, file='VPtrait_wg.RDATA')


# #
# #    Abundance


setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_ab_wg.RDATA')
dim(VP.test3_wg$Vnorm)
dim(VP.test3_wg$Cnorm)

VP.test3 = array(NA, dim = dim(VP.test3_wg$Cnorm))

for (i in 1:dim(VP.test3_wg$Vnorm)[4]){  ## mcmc samples
  for (j in 1:dim(VP.test3_wg$Vnorm)[3]){  ## species
    
    VP.test3[,,j,i] = VP.test3_wg$Cnorm[,,j,i]
    diag(VP.test3[,,j,i]) = VP.test3_wg$Vnorm[,,j,i]
    
    colnames( VP.test3[,,j,i]) = colnames(VP.test3_wg$Cnorm[,,j,i])
    rownames( VP.test3[,,j,i]) = rownames(VP.test3_wg$Cnorm[,,j,i])
    
  }
}

# prepare scaling based on tjurR2
varr <- aperm(array(TabMF_wg$SR2, dim = c(15L, 3L, 3L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)

# only keeping the Normalized variance
# Posterior mean
Mean.fix.spwg = apply(varr*VP.test3[1:3, 1:3, ,], c(1, 2,3), mean, na.rm=T)
VPclim = 100*Mean.fix.spwg[1,1,]
VPhabC = 100*Mean.fix.spwg[2,2,]
VPhabP = 100*Mean.fix.spwg[3,3,]

wgVP_ab = data.frame(VPClim = VPclim,
                  VPhabC = VPhabC,
                  VPhabP = VPhabP)

dim(wintG_trait_sub)
dim(wgVP_ab)
wgVP_ab$iucn2020_binomial = wg.sp


wgDat_trait_ab = merge(wgVP_ab, wintG_trait_sub, by= 'iucn2020_binomial')

setwd("~/Data/4.VPanalyses/VP_cond_trait")
save(wgDat_trait_ab, file='VPtrait_ab_wg.RDATA')
