#
#           Link cov and corr
#           and look at the VCMN variance to investigate 
#           covariance effects
# ...............................................................

library(dplyr)
library(tidyr)
library(abind)
library(quantmod)
library(tibble)
library(reshape2)

## load all VCNMP variance partitions
setwd("D:/Helsinki/RECcII/CSC_HMSC/VPres")
load('VPCP_bd.RDATA')
load('VPCP_bf.RDATA')
load('VPCP_rod.RDATA')
load('VPCP_wg.RDATA')
load('VPCP_moth.RDATA')

## load all linear predictor correlation
load('FRcorr_taxa1010.RDATA')

# example dimensions
dim(VP.test3_bf$Vnorm) #    1  /   predictors group  /  species /  mcmc samples
dim(FRcor_bf) #    predictors group (+ effort in our case)  /  predictors group (+ effort in our case)  /  species /  mcmc samples
groupnames1=c("Effort", "Climate", "HoccitatConf", 'Hoccitatprop')      # + 3 random effect so 7 in total


# Scaling by R2 all VP summaries
# 
# rod
varr <- aperm(array(TabMF_rod$TjurR2, dim = c(1000L, 8L, 6L, 1L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)
VP.test3_rod$Vnorm = varr * VP.test3_rod$Vnorm
VP.test3_rod$Vdiag = varr * VP.test3_rod$Vdiag


## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_rod$TjurR2, dim = c(1000L, 8L, 6L, 6L)), perm = c(4L, 3L, 2L, 1L))
VP.test3_rod$Cnorm = Carr * VP.test3_rod$Cnorm
VP.test3_rod$Vmarg = Carr * VP.test3_rod$Vmarg
VP.test3_rod$Vpart = Carr * VP.test3_rod$Vpart

# bf
varr <- aperm(array(TabMF_bf$TjurR2, dim = c(1000L, 57L, 6L, 1L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)
VP.test3_bf$Vnorm = varr * VP.test3_bf$Vnorm
VP.test3_bf$Vdiag = varr * VP.test3_bf$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_bf$TjurR2, dim = c(1000L, 57L, 6L, 6L)), perm = c(4L, 3L, 2L, 1L))
VP.test3_bf$Cnorm = Carr * VP.test3_bf$Cnorm
VP.test3_bf$Vmarg = Carr * VP.test3_bf$Vmarg
VP.test3_bf$Vpart = Carr * VP.test3_bf$Vpart

# moth
varr <- aperm(array(TabMF_Moth$TjurR2, dim = c(1000L, 319L, 6L, 1L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)
VP.test3_moth$Vnorm = varr * VP.test3_moth$Vnorm
VP.test3_moth$Vdiag = varr * VP.test3_moth$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_Moth$TjurR2, dim = c(1000L, 319L, 6L, 6L)), perm = c(4L, 3L, 2L, 1L))
VP.test3_moth$Cnorm = Carr * VP.test3_moth$Cnorm
VP.test3_moth$Vmarg = Carr * VP.test3_moth$Vmarg
VP.test3_moth$Vpart = Carr * VP.test3_moth$Vpart

# wg
varr <- aperm(array(TabMF_wg$TjurR2, dim = c(1000L, 15L, 6L, 1L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)
VP.test3_wg$Vnorm = varr * VP.test3_wg$Vnorm
VP.test3_wg$Vdiag = varr * VP.test3_wg$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_wg$TjurR2, dim = c(1000L, 15L, 6L, 6L)), perm = c(4L, 3L, 2L, 1L))
VP.test3_wg$Cnorm = Carr * VP.test3_wg$Cnorm
VP.test3_wg$Vmarg = Carr * VP.test3_wg$Vmarg
VP.test3_wg$Vpart = Carr * VP.test3_wg$Vpart


# bird
varr <- aperm(array(TabMF_bd$TjurR2, dim = c(1000L, 102L, 6L, 1L)), perm = c(4L, 3L, 2L, 1L))
dim(varr)
VP.test3_bd$Vnorm = varr * VP.test3_bd$Vnorm
VP.test3_bd$Vdiag = varr * VP.test3_bd$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_bd$TjurR2, dim = c(1000L, 102L, 6L, 6L)), perm = c(4L, 3L, 2L, 1L))
VP.test3_bd$Cnorm = Carr * VP.test3_bd$Cnorm
VP.test3_bd$Vmarg = Carr * VP.test3_bd$Vmarg
VP.test3_bd$Vpart = Carr * VP.test3_bd$Vpart


## combine all taxa
# ..........................................

## combined covariance
VP.list = list(VP.test3_bf, VP.test3_bd, VP.test3_wg, VP.test3_rod, VP.test3_moth)

## combined correlation
## removing the row and line occout effort
Matcorr = list(FRcor_bf[-1, -1,,], FRcor_bd[-1, -1,,], FRcor_wg[-1, -1,,], FRcor_rod[-1, -1,,], FRcor_moth[-1, -1,,])

# Link to all the vp values  
# ------------------------------------

##.......................
## define cov and corr thresholds
cvt = 0.01 # 1%
taxa_list = c('Butterfly', 'Bird', 'Large mammal', 'Rodent', 'Moth')

VP_summary = vector("list", 5) # 5 taxa and 10 probs

for (j in 1:5) {
  print(paste('*********    Taxa', j, sep=''))
  species_j = names(VP.list[[j]]$Vnorm[1,1,,1])
  
  Mat_post = list()
  
  for (k in 1:dim(VP.list[[j]]$Vnorm[,,,])[3]) { # run longer for 1000
    print(paste('It.  ', k, sep=''))
    
    Mat_taxa = list()
    
    for (i in 1:length(species_j)) {
      
      ## first we prepare the matrix that will have each quantity
      ## we fill the covariance from the lower triangular matrix (same as the upper triangular)
      mat = VP.list[[j]]$Cnorm[,,i,k]
      mat[lower.tri(mat)] <- t(mat)[lower.tri(t(mat))]
      
      Mat_taxa_cov = melt(mat)  # 
      colnames(Mat_taxa_cov) = c("row_name", 'col_name', 'Cnorm')
      
      ## remove the predictor to predictor (i.e. the VP because it is NA in Cnorm)
      Mat_taxa_cov = Mat_taxa_cov[!(Mat_taxa_cov$row_name == Mat_taxa_cov$col_name) ,]
      
      # Species i & posterior sample k
      Vnorm_cc = VP.list[[j]]$Vnorm[,,i, k]
      Vmarg_cc = VP.list[[j]]$Vmarg[,,i, k]
      Vpart_cc = VP.list[[j]]$Vpart[,,i, k]
      
      for (l in 1:dim(Mat_taxa_cov)[1]) {
        Mat_taxa_cov$V1[l] = Vnorm_cc[names(Vnorm_cc) == Mat_taxa_cov$row_name[l] ]
        Mat_taxa_cov$M1[l]  = Vmarg_cc[rownames(Vmarg_cc) == Mat_taxa_cov$row_name[l], colnames(Vmarg_cc) == Mat_taxa_cov$col_name[l]]
        Mat_taxa_cov$P1[l]  = Vpart_cc[rownames(Vpart_cc) == Mat_taxa_cov$row_name[l], colnames(Vpart_cc) == Mat_taxa_cov$col_name[l]]
      }
      
      Mat_taxa_cov$VP1 = Mat_taxa_cov$V1/Mat_taxa_cov$P1
      Mat_taxa_cov$VM1 = Mat_taxa_cov$M1/Mat_taxa_cov$V1
      
      # add number of species for later proportion calculation
      Mat_taxa_cov$NBsp = dim(VP.list[[j]]$Vnorm)[3]
      
      ## add correlation between predictor to the matrix
      Matcorr.temp = Matcorr[[j]][,,i,k]
      rownames(Matcorr.temp) = colnames(VP.list[[j]]$Vnorm)
      colnames(Matcorr.temp) = colnames(VP.list[[j]]$Vnorm)
      Mat_taxa_corr = melt(Matcorr.temp)  # 
      colnames(Mat_taxa_corr) = c("row_name", 'col_name', 'corr')
      
      ## remove the predictor to predictor (i.e. the VP because it is NA in Cnorm)
      Mat_taxa_corr = Mat_taxa_corr[!(Mat_taxa_corr$corr == 1),]
      
      dim(Mat_taxa_corr)## same order and dim
      dim(Mat_taxa_cov)
      
      Mat_taxa_cov$corr = Mat_taxa_corr$corr
      
      ## since we have set up the case1 and case 2 we don't need to do it twice to have for example climat hocc and hocc climate covariation understaning
      ## thus we filter out
      
      ## Order by cases
      #----------------------------------------
      # ------------------------------------------------------
      
      Mat_taxa_cov$case = NA
      
      # define some thresholds for comparison
      ## M ≈ 0
      M0 = 0.01  # 1%
      
      #................................................................
      
      Mat_taxa_cov = Mat_taxa_cov %>% mutate(case = case_when(
        #abs(Cnorm) < cvt ~ 'NA', # at least 0.01 rescaled covariance
        abs(corr) < cvt ~ 'Independent', # same but for correlation
        M1 > V1  & abs(V1) >= M0 & V1 > P1 ~ 'Amplify',
        abs(M1) < V1 & P1 < V1 ~ 'Suppress',
        abs(M1) <= M0 ~ 'Confounded',
        abs(VM1) > (1-M0) & abs(VM1) > (1+M0) & abs(VP1) > (1-M0) & abs(VP1) > (1 + M0)  ~ 'Independent',
        TRUE ~ "Other"
      )) 
      Mat_taxa[[i]] = Mat_taxa_cov
    }
    
    ## calculate proportion of species for each case
    Mat_taxa_count = do.call(rbind, Mat_taxa)
    
    Mat_taxa_prop = Mat_taxa_count %>%
      select(row_name, col_name, NBsp, case) %>%
      pivot_longer(case) %>%
      count(col_name, row_name, NBsp, name, value)
    
    Mat_taxa_prop$prop = Mat_taxa_prop$n / Mat_taxa_prop$NBsp
    
    Mat_taxa_prop$iteration = k
    Mat_post[[k]] = Mat_taxa_prop
  }
  
  # get posterior distribution of proportion of species for each case:
  VP_summary_prop = do.call(rbind, Mat_post)
  
  #.......................................................................................................
  
  ## save information that does not vary across iteration
  VP_summary_prop$taxa = taxa_list[[j]]
  VP_summary[[j]] = VP_summary_prop
  
}


length(VP_summary) ## 

VP_taxa = do.call(rbind, VP_summary)
dim(VP_taxa)
table(VP_taxa$value)

#  Amplify  Confounded Independent       Other    Suppress 
#  132732      114421      129528       25716      143354 


#.......................................................................................................

## combine all random
table(VP_taxa[VP_taxa$col_name == 'year',]$taxa)  # 
table(VP_taxa[VP_taxa$col_name == 'site',]$taxa)  # 
table(VP_taxa[VP_taxa$col_name == 'bg',]$taxa)  # 
## and remove between random
VP_taxa = VP_taxa[!(VP_taxa$col_name == VP_taxa$row_name) ,]


VP_C_all = VP_taxa[VP_taxa$value %in% c("Amplify", "Confounded", "Independent", "Suppress"),]

VP_C_all$row_name = factor(VP_C_all$row_name)
VP_C_all$col_name = factor(VP_C_all$col_name)


table(VP_C_all$value)

VP_C_all$row_name = factor(VP_C_all$row_name, 
                           levels = c("Climate", "HabitatConf", "Habitatprop", 'year', 'site', 'bg'),
                           labels = c("Climate", "LandscapeConf", "Habitatcomp",  "Random",  "Random",  "Random"))
VP_C_all$col_name = factor(VP_C_all$col_name, 
                           levels = c("Climate", "HabitatConf", "Habitatprop", 'year', 'site', 'bg'),
                           labels = c("Climate", "LandscapeConf", "Habitatcomp",  "Random",  "Random",  "Random"))
# remove the remaining random to random evaluations
VP_C_all = VP_C_all[!(VP_C_all$col_name == VP_C_all$row_name) ,]



VP_C_all$taxa = factor(VP_C_all$taxa, 
                       levels = c('Bird', 'Butterfly', 'Moth', 'Rodent','Large mammal'))


# prepare the plotting
rc_un = data.frame(col_name= c("Climate", "LandscapeConf", "Habitatcomp", "Random"),
                   row_name=c("Climate", "LandscapeConf", "Habitatcomp", "Random"))
rc_un$row_name = factor(rc_un$row_name, 
                        levels = c("Climate", "LandscapeConf","Habitatcomp",  "Random"))
rc_un$col_name = factor(rc_un$col_name, 
                        levels = c("Climate", "LandscapeConf", "Habitatcomp", "Random" ))



VP_C_all_occ = VP_C_all


setwd("D:/Helsinki/RECcII/CSC_HMSC/VPres")
save(VP_C_all_occ, file = 'VP_C1resc_all_occ1111.RDATA')


_C_all_occ = VP_C_all


setwd("D:/Helsinki/RECcII/Joint_interpretation_part")
save(VP_C_all_occ, file = 'VP_C1resc_all_occ1010.RDATA')


