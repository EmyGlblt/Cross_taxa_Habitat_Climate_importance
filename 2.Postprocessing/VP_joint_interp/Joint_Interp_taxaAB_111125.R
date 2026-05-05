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
setwd("~/Data/4.VPanalyses/VPCP")
load('VPCP_ab_bd.RDATA')
load('VPCP_ab_bf.RDATA')
load('VPCP_ab_rod.RDATA')
load('VPCP_ab_wg.RDATA')
load('VPCP_ab_moth.RDATA')

## load all linear predictor correlation
setwd("~/Data/4.VPanalyses/LF_corr")

load('FRcorr_taxa1010.RDATA')

# example dimensions
dim(VP.test3_bf$Vnorm) #    1  /   predictors group  /  species /  mcmc samples
dim(FRcorAB_bf) #    predictors group (+ effort in our case)  /  predictors group (+ effort in our case)  /  species /  mcmc samples
groupnames1=c("Effort", "Climate", "HabitatConf", 'Habitatprop')      # + 3 random effect so 7 in total


# Scaling by R2 all VP summaries
# 
# rod
varr <- aperm(array(TabMF_rod$SR2, dim = c(8L, 6L, 1L, 2000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_rod$Vnorm = varr * VP.test3_rod$Vnorm
VP.test3_rod$Vdiag = varr * VP.test3_rod$Vdiag


## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_rod$SR2, dim = c(8L, 6L, 6L, 2000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_rod$Cnorm = Carr * VP.test3_rod$Cnorm
VP.test3_rod$Vmarg = Carr * VP.test3_rod$Vmarg
VP.test3_rod$Vpart = Carr * VP.test3_rod$Vpart

# bf
varr <- aperm(array(TabMF_bf$SR2, dim = c(57L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_bf$Vnorm = varr * VP.test3_bf$Vnorm
VP.test3_bf$Vdiag = varr * VP.test3_bf$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_bf$SR2, dim = c(57L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_bf$Cnorm = Carr * VP.test3_bf$Cnorm
VP.test3_bf$Vmarg = Carr * VP.test3_bf$Vmarg
VP.test3_bf$Vpart = Carr * VP.test3_bf$Vpart

# moth
varr <- aperm(array(TabMF_Moth$SR2, dim = c(319L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_moth$Vnorm = varr * VP.test3_moth$Vnorm
VP.test3_moth$Vdiag = varr * VP.test3_moth$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_Moth$SR2, dim = c(319L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_moth$Cnorm = Carr * VP.test3_moth$Cnorm
VP.test3_moth$Vmarg = Carr * VP.test3_moth$Vmarg
VP.test3_moth$Vpart = Carr * VP.test3_moth$Vpart

# wg
varr <- aperm(array(TabMF_wg$SR2, dim = c(15L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_wg$Vnorm = varr * VP.test3_wg$Vnorm
VP.test3_wg$Vdiag = varr * VP.test3_wg$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_wg$SR2, dim = c(15L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_wg$Cnorm = Carr * VP.test3_wg$Cnorm
VP.test3_wg$Vmarg = Carr * VP.test3_wg$Vmarg
VP.test3_wg$Vpart = Carr * VP.test3_wg$Vpart


# bird
varr <- aperm(array(TabMF_bd$SR2, dim = c(102L, 6L, 1L, 1000L)), perm = c(3L, 2L, 1L, 4L))
dim(varr)
VP.test3_bd$Vnorm = varr * VP.test3_bd$Vnorm
VP.test3_bd$Vdiag = varr * VP.test3_bd$Vdiag

## Which species have large covariance and which covariance are they?
Carr <- aperm(array(TabMF_bd$SR2, dim = c(102L, 6L, 6L, 1000L)), perm = c(3L, 2L, 1L, 4L))
VP.test3_bd$Cnorm = Carr * VP.test3_bd$Cnorm
VP.test3_bd$Vmarg = Carr * VP.test3_bd$Vmarg
VP.test3_bd$Vpart = Carr * VP.test3_bd$Vpart



## combine all taxa
# ..........................................

## combined covariance
VP.list = list(VP.test3_bf, VP.test3_bd, VP.test3_wg, VP.test3_rod, VP.test3_moth)

## combined correlation
## removing the row and line about effort
Matcorr = list(FRcorAB_bf[-1, -1,,], FRcorAB_bd[-1, -1,,], FRcorAB_wg[-1, -1,,], FRcorAB_rod[-1, -1,,], FRcorAB_moth[-1, -1,,])

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
      
      ## since we have set up the case1 and case 2 we don't need to do it twice to have for example climat hab and hab climate covariation understaning
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
# 131182      124587      143683      115454      142886 



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



VP_C_all_ab = VP_C_all


setwd("D:/Helsinki/RECcII/CSC_HMSC/VPres")
save(VP_C_all_ab, file = 'VP_Cresc_all_ab1111.RDATA')





# Combine plotting from the different posterior probability

Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")

colnames(VP_C_all_ab)
head(VP_C_all_ab)


VP_C_all_ab$prop = as.numeric(VP_C_all_ab$prop)
table(VP_C_all_ab$value)

## 1 confounded 
# 3 independent

pal = c('#66C2A5', 'orange3', 'yellow3', 'gray60')

library(grid)
library(gtable)


gg_old = {ggplot() +
    geom_rect(data = rc_un, aes(fill = row_name), fill='darkgrey', 
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5)+ 
    #geom_point(data = VP_C_all_ab, aes(y=value, x=prop, color=taxa)) +
    geom_density_ridges(data = VP_C_all_ab, aes(y=value, x=prop, fill=taxa), 
                        alpha=0.6, quantile_lines = TRUE, trim = TRUE) +
    xlim(0, 1) +
    #scale_fill_manual(values = Newcol) +
    scale_fill_manual(values = Newcol) +
    facet_grid(row_name ~ col_name, drop = FALSE, switch = "y")+
    theme(panel.grid.major.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())}%>%
  modify_facet_appearance(strip.background.x.fill = pal,
                          strip.background.x.col = rep('black', 4),
                          strip.text.x.col = rep('black', 4),
                          strip.background.y.fill = pal,
                          strip.background.y.col =rep('black', 4),
                          strip.text.y.col = rep('black', 4))

library(ggridges)

VP_C_all_ab.bd = VP_C_all_ab[VP_C_all_ab$taxa == 'Bird',]

ggplot(data = VP_C_all_ab.bd, aes(y=prop, x=value, fill=taxa)) +
  #geom_point(data = VP_C_all_ab, aes(y=prop, x=value, color=taxa, size = Prob)) +
  geom_violin(position = "dodge", trim = T, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = Newcol) +
  facet_grid(row_name ~ col_name, drop = FALSE, switch = "y") + 
  coord_polar(start = -pi/12) +
  ylim(0, 0.75) +
  geom_hline(yintercept = seq(0, 1, by = 0.25), colour = "darkgrey", size = 0.35) +
  #geom_vline(xintercept = seq(1, 4, by = 1), colour = "darkgrey", size = 0.35) +
  labs(y= "Joint effect", x = "% of species") +
  theme(panel.grid.major.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = c(-48, 48, -48, 48)),
        panel.ontop =FALSE)+
  guides(color=guide_legend(title="Taxa", override.aes=list(size=4)))


gp = {ggplot() +
    #geom_point(data = VP_C_all_ab, aes(y=prop, x=value, color=taxa, size = Prob)) +
    geom_violin(data = VP_C_all_ab, aes(y=prop, x=value, fill=taxa)) +
    # geom_jitter(height = 0, width = 0.1) +
    # coord_cartesian(ylim=c(0, 30)) +
    scale_fill_manual(values = Newcol) +
    facet_grid(row_name ~ col_name, drop = FALSE, switch = "y") + 
    coord_flip() +
    coord_polar(start = -pi/12) +
    geom_hline(yintercept = seq(0, 30, by = 10), colour = "darkgrey", size = 0.35) +
    #geom_vline(xintercept = seq(1, 4, by = 1), colour = "darkgrey", size = 0.35) +
    labs(y= "Joint effect", x = "% of species") +
    theme(panel.grid.major.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.ontop =FALSE)+
    guides(color=guide_legend(title="Taxa", override.aes=list(size=4)))}%>%
  modify_facet_appearance(strip.background.x.fill = pal,
                          strip.background.x.col = rep('black', 4),
                          strip.text.x.col = rep('black', 4),
                          strip.background.y.fill = pal,
                          strip.background.y.col =rep('black', 4),
                          strip.text.y.col = rep('black', 4))

plot(gp)

# Remove facets
idx <- which(gp$layout$name %in% c("panel-1-1", "panel-2-2", "panel-3-3", "panel-4-4"))
for (i in idx) gp$grobs[[i]] <- gg_old$grobs[[i]]

## plot again
grid.newpage()
grid.draw(gp)+
  geom_rect(data = rc_un, aes(fill = row_name), fill='black', 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5)

trip.text.y.col = rep('black', 4))

plot(gp)

# Remove facets
idx <- which(gp$layout$name %in% c("panel-1-1", "panel-2-2", "panel-3-3", "panel-4-4"))
for (i in idx) gp$grobs[[i]] <- gg_old$grobs[[i]]

## plot again
grid.newpage()
grid.draw(gp)+
  geom_rect(data = rc_un, aes(fill = row_name), fill='black', 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5)

0.5)

trip.text.y.col = rep('black', 4))

plot(gp)

# Remove facets
idx <- which(gp$layout$name %in% c("panel-1-1", "panel-2-2", "panel-3-3", "panel-4-4"))
for (i in idx) gp$grobs[[i]] <- gg_old$grobs[[i]]

## plot again
grid.newpage()
grid.draw(gp)+
  geom_rect(data = rc_un, aes(fill = row_name), fill='black', 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5)

