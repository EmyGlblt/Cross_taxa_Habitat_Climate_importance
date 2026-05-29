library(dplyr)
library(tidyr)
library(abind)
#library(quantmod)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggridges)
library(tidybayes)
library(scales)

library(gt)

# Occurrence ...................................................................................................................
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
    
    mat.bd =   sweep(
      VP.test3_bd$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_bd$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )

    # to prevent unrealistic values
    # if M is too small
    mat.bd[abs(VP.test3_bd$Vmarg[,,j,i])<0.01] = -1#*sign(mat.bd[abs(VP.test3_bd$Vmarg[,,j,i])<0.01])
    
    mat.diff.bd[,,j,i] = mat.bd
  }
}
dim(mat.diff.bd)
dimnames(mat.diff.bd) = dimnames(VP.test3_bd$Vmarg)


mat.diff.bf = array(data=NA, dim = dim(VP.test3_bf$Vmarg))
for (i in 1:dim(VP.test3_bf$Vmarg)[4]){
  for (j in 1:dim(VP.test3_bf$Vmarg)[3]) {
    
    mat.bf =   sweep(
      VP.test3_bf$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_bf$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.bf[abs(VP.test3_bf$Vmarg[,,j,i])<0.01] = -1#*sign(mat.bf[abs(VP.test3_bf$Vmarg[,,j,i])<0.01])
    
    mat.diff.bf[,,j,i] = mat.bf    
  }
}
dimnames(mat.diff.bf) = dimnames(VP.test3_bf$Vmarg)



mat.diff.rod = array(data=NA, dim = dim(VP.test3_rod$Vmarg))
for (i in 1:dim(VP.test3_rod$Vmarg)[4]){
  for (j in 1:dim(VP.test3_rod$Vmarg)[3]) {
    
    mat.rod =   sweep(
      VP.test3_rod$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_rod$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.rod[abs(VP.test3_rod$Vmarg[,,j,i])<0.01] = -1#*sign(mat.rod[abs(VP.test3_rod$Vmarg[,,j,i])<0.01])
    
    mat.diff.rod[,,j,i] = mat.rod    
  }
}
dimnames(mat.diff.rod) = dimnames(VP.test3_rod$Vmarg)


mat.diff.moth = array(data=NA, dim = dim(VP.test3_moth$Vmarg))
for (i in 1:dim(VP.test3_moth$Vmarg)[4]){
  for (j in 1:dim(VP.test3_moth$Vmarg)[3]) {
    mat.moth =   sweep(
      VP.test3_moth$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_moth$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.moth[abs(VP.test3_moth$Vmarg[,,j,i])<0.01] = -1#*sign(mat.moth[abs(VP.test3_moth$Vmarg[,,j,i])<0.01])
    
    mat.diff.moth[,,j,i] = mat.moth    
  }
}
dimnames(mat.diff.moth) = dimnames(VP.test3_moth$Vmarg)


mat.diff.wg = array(data=NA, dim = dim(VP.test3_wg$Vmarg))
for (i in 1:dim(VP.test3_wg$Vmarg)[4]){
  for (j in 1:dim(VP.test3_wg$Vmarg)[3]) {
    
    mat.wg =   sweep(
      VP.test3_wg$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_wg$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.wg[abs(VP.test3_wg$Vmarg[,,j,i])<0.01] = -1#*sign(mat.wg[abs(VP.test3_wg$Vmarg[,,j,i])<0.01])
    
    mat.diff.wg[,,j,i] = mat.wg    
  }
}
dimnames(mat.diff.wg) = dimnames(VP.test3_wg$Vmarg)


## posterior mean

# re do mean after excluding extreme values?
mat.diff.bd[mat.diff.bd>1] = NA
DiffMV_bd = melt(apply(mat.diff.bd, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_bd)[1] = 'var'
DiffMV_bd$Taxa = 'birds'

mat.diff.bf[mat.diff.bf>1] = NA
DiffMV_bf = melt(apply(mat.diff.bf, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_bf)[1] = 'var'
DiffMV_bf$Taxa = 'butterflies'

mat.diff.moth[mat.diff.moth>1] = NA
DiffMV_moth = melt(apply(mat.diff.moth, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_moth)[1] = 'var'
DiffMV_moth$Taxa = 'moths'

mat.diff.rod[mat.diff.rod>1] = NA
DiffMV_rod = melt(apply(mat.diff.rod, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_rod)[1] = 'var'
DiffMV_rod$Taxa = 'small mammals'

mat.diff.wg[mat.diff.wg>1] = NA
DiffMV_wg = melt(apply(mat.diff.wg, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_wg)[1] = 'var'
DiffMV_wg$Taxa = 'large mammals'


## combine all and keep only the shared importance of interests

DiffMV = rbind(DiffMV_bd, DiffMV_bf, DiffMV_moth, DiffMV_rod, DiffMV_wg)
DiffMV = DiffMV[!(DiffMV$var == DiffMV$element) ,]
DiffMV = DiffMV[!(DiffMV$var %in% c('site', 'year', 'bg') | DiffMV$element %in% c('site', 'year', 'bg')) ,]

DiffMV$Taxa = factor(DiffMV$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
DiffMV$var = factor(DiffMV$var, levels = c("Climate", "HabitatConf","Habitatprop"),
                      labels = c("Climate", "LandscapeConf","Habitatcomp"))
DiffMV$element = factor(DiffMV$element, levels = c("Climate", "HabitatConf","Habitatprop"),
                        labels = c("Climate", "LandscapeConf","Habitatcomp"))



summary(DiffMV)
dim(DiffMV) #3006  5

DiffMV2 = DiffMV[DiffMV$value < 1,] # now only very few excluded
summary(DiffMV2)
dim(DiffMV2) # 2897



# Get the percentage
# 
confoun = na.omit(DiffMV[DiffMV$value < -0.99,])
confoun %>%
  group_by(var, element) %>%
  summarise(count = 100*length(value)/(102+57+319+8+15))

confoun = confoun %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)
confoun$tot = c(rep(102, 2), rep(319, 2), rep(15, 2)) # check for species if a driver1-driver2 confounding exist and adapt the numbers of rep (ex none for small mammals)
confoun$pct = 100*confoun$count/confoun$tot
confoun %>% 
  print(n = 100)

confoun$shared = 'confounded'


#.......................................................
indepdt = na.omit(DiffMV[DiffMV$value > -0.01 & DiffMV$value < 0.01,])
indepdt %>%
  group_by(var, element) %>%
  summarise(count = 100*length(value)/(102+57+319+8+15))

indepdt = indepdt %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)
indepdt$tot = c(rep(102, 5), rep(57, 3), rep(319, 5), rep(15, 1), rep(15, 2))  # check for species if a driver1-driver2 independence exist and adapt the numbers of rep (ex none for small mammals)
indepdt$pct = 100*indepdt$count/indepdt$tot
indepdt %>% 
  print(n = 100)

indepdt$shared = 'Independent'


#.......................................................

Suppress = na.omit(DiffMV[DiffMV$value > -0.99 & DiffMV$value < -0.01,])
Suppress = Suppress %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)

Suppress$tot = c(rep(102, 6), rep(57, 6), rep(319, 6), rep(8, 6), rep(15, 6))
Suppress$pct = 100*Suppress$count/Suppress$tot
Suppress %>% 
  print(n = 100)

Suppress$shared = 'Suppressing'


#.......................................................

reinf = na.omit(DiffMV[DiffMV$value > 0.01 & DiffMV$value < 1,])
reinf = reinf %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)

reinf$tot = c(rep(102, 6), rep(57, 4), rep(319, 6), rep(8, 1), 
              rep(15, 6))
reinf$pct = 100*reinf$count/reinf$tot
reinf %>% 
  print(n = 100)

reinf$shared = 'Reinforcing'


# altogether
Effects = rbind(reinf, indepdt, Suppress, confoun)
Effects$shared = factor(Effects$shared,
                        levels = c('Reinforcing', 'Independent', 'Suppressing', 'confounded'))


gt_table <- Effects %>%
  gt(groupname_col = "shared", rowname_col = 'Taxa', row_group_as_column = FALSE)%>%
  tab_spanner_delim(delim = "_") |>
  tab_header(
    title = "Occurrence data",
  ) %>%
  fmt_number(
    columns = c(pct),
    decimals = 2
  ) %>%
  cols_label(
    var = "Variable 1",
    element = "Variable 2",
    count = "Species #",
    tot = "Total #",
    pct = "Pct (%)"
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(everything())
  ) %>% tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3), # A bit more padding for summaries
    row_group.padding = px(6)    # And even more for our groups
  ) |> 
  opt_stylize(style = 3, color = 'gray')

# Print the table
gt_table |> gtsave("Occ_shared.png", expand = 10)



#.......................................................
#.......................................................
# row-wise

# prepare data properly
nb.sp = c('birds'=102, 'butterflies'=57, 'moths'=319,
          'small mammals'= 8, 'large mammals' = 15)
totals_df <- tibble(
  var = c(rep("Climate", 5), rep("LandscapeConf", 5), rep("Habitatcomp",5)),
  Taxa = c(rep(levels(DiffMV$Taxa), 3)),
  tot = c(rep(nb.sp, 3))
)
group_taxa = levels(DiffMV$Taxa)


# Prepare different bining to highlight independent and confounding effects
DiffMV = na.omit(DiffMV)

df_plot_indptconf <- DiffMV[(DiffMV$value >= -0.01) & (DiffMV$value <= 0.01) | (DiffMV$value <= -0.99),] %>%
  mutate(
    vbin = floor(value / 0.01) * 0.01
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

# öarger bins for others
df_plot_shared <- DiffMV[!((DiffMV$value> -0.01) & (DiffMV$value< 0.01)) & !(DiffMV$value < -0.99) ,] %>%
  mutate(
    vbin = floor(value / 0.05) * 0.05
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

df_plot = rbind(df_plot_indptconf, df_plot_shared)

df_plot$Taxa = factor(df_plot$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
df_plot$var = factor(df_plot$var, levels = c("Climate", "LandscapeConf","Habitatcomp"))
df_plot$element = factor(df_plot$element, levels = c("Climate", "LandscapeConf","Habitatcomp"))


plot = ggplot(df_plot, aes(vbin, pct, fill = element)) + 
  geom_col(width = 0.05, position = 'identity') +  
  scale_fill_manual(values = c(alpha("#66C2A5",0.4), alpha("yellow3",0.4), alpha('orange3',0.4))) +
  ylab('Species proportion, %') + xlab('Relative difference of the marginal variance \n partition (M) and the variance partition (V)') +
  geom_vline(xintercept = c(0), colour='gray50', linetype = "dashed") + ## Independent
  geom_vline(xintercept = c(-1, -0.99, 1), colour='gray50', linetype = "dashed") + ## confounded
  scale_x_continuous(
    limits = c(-1.025, 1.025),
    breaks = seq(-1, 1, by = 0.5),
    expand = c(0.05, 0)
  ) +  
  facet_grid(var ~ Taxa) +
  coord_flip(ylim = c(0, 50), clip = "on") +
  theme_classic() + theme(text = element_text(size = 18),
                          axis.text.x = element_text(angle = 45, hjust = 1),
                          axis.title.y = element_text(hjust = 0),
                          legend.position = "bottom",
                          panel.background = element_rect(fill = NA, color = "black")) +
  guides(fill=guide_legend(title=" "))



g <- ggplot_gtable(ggplot_build(plot))
strips <- which(grepl('strip-', g$layout$name))

#col_strip = c(rep("white",5), alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5))
#col_text = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", rep("black",3))

col_strip = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5))
col_text = c(rep("white",5), rep("black",3))

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_strip[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$col <- col_text[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- col_text[i]
  
}

plot(g)



# with random effects
#

## combine all and keep only the shared importance of interests
DiffMV = rbind(DiffMV_bd, DiffMV_bf, DiffMV_moth, DiffMV_rod, DiffMV_wg)
DiffMV = DiffMV[!(DiffMV$var == DiffMV$element) ,]

DiffMV$Taxa = factor(DiffMV$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))

DiffMV <- na.omit(DiffMV) %>%
  mutate(element = ifelse(element %in% c('year', 'site', 'bg'), "random", element)) %>%
  group_by(Taxa, var, element, Species) %>%
  summarise(value = mean(value, na.rm=T), .groups = "drop")

DiffMV <- na.omit(DiffMV) %>%
  mutate(var = ifelse(var %in% c('year', 'site', 'bg'), "random", var)) %>%
  group_by(Taxa, var, element, Species) %>%
  summarise(value = mean(value, na.rm=T), .groups = "drop")

DiffMV$var = factor(DiffMV$var, 
                    levels = c("1", "2", "3", 'random'),
                    labels = c("Climate", "LandscapeConf", "Habitatcomp",  "Random"))

DiffMV$element = factor(DiffMV$element, levels = c("1", "2", "3", 'random'),
                        labels = c("Climate", "LandscapeConf", "Habitatcomp",  "Random"))

# Since we have multiple random effects here let's summarize it for simplicity


totals_df <- tibble(
  var = c(rep("Climate", 5), rep("LandscapeConf", 5), rep("Habitatcomp",5), rep("Random",5)),
  Taxa = c(rep(levels(DiffMV$Taxa), 4)),
  tot = c(rep(nb.sp, 4))
)
group_taxa = levels(DiffMV$Taxa)

# Prepare different bining to highlight independent and confounding effects
DiffMV = na.omit(DiffMV)

df_plot_indptconf <- DiffMV[(DiffMV$value >= -0.01) & (DiffMV$value <= 0.01) | (DiffMV$value <= -0.99),] %>%
  mutate(
    vbin = floor(value / 0.01) * 0.01
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

# öarger bins for others
df_plot_shared <- DiffMV[!((DiffMV$value> -0.01) & (DiffMV$value< 0.01)) & !(DiffMV$value < -0.99) ,] %>%
  mutate(
    vbin = floor(value / 0.05) * 0.05
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

df_plot = rbind(df_plot_indptconf, df_plot_shared)

df_plot$Taxa = factor(df_plot$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
df_plot$var = factor(df_plot$var, levels = c("Climate", "LandscapeConf","Habitatcomp", "Random"))
df_plot$element = factor(df_plot$element, levels = c("Climate", "LandscapeConf","Habitatcomp", "Random"))



# row-wise


plot =  ggplot(df_plot, aes(vbin, pct, fill = element)) + 
  geom_col(width = 0.05, position = 'identity') +  
  scale_fill_manual(values = c(alpha("#66C2A5",0.4), alpha("yellow3",0.4), 
                               alpha('orange3',0.4), alpha('gray60',0.4))) +
  ylab('Species proportion, %') + xlab('Relative difference of the marginal variance \n partition (M) and the variance partition (V)') +
  geom_vline(xintercept = c(0), colour='gray50', linetype = "dashed") + ## Independent
  geom_vline(xintercept = c(-1, -0.99, 1), colour='gray50', linetype = "dashed") + ## confounded
  scale_x_continuous(
    limits = c(-1.025, 1.025),
    breaks = seq(-1, 1, by = 0.5),
    expand = c(0.05, 0)
  ) +  
  facet_grid(var ~ Taxa) +
  coord_flip(ylim = c(0, 80), clip = "on") +
  theme_classic() + theme(text = element_text(size = 18),
                          axis.text.x = element_text(angle = 45, hjust = 1),
                          legend.position = "bottom",
                          axis.title.y = element_text(hjust = 0),
                          panel.background = element_rect(fill = NA, color = "black")) +
  guides(fill=guide_legend(title=" "))




g <- ggplot_gtable(ggplot_build(plot))
strips <- which(grepl('strip-', g$layout$name))

#col_strip = c(rep("white",5), alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5))
#col_text = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", rep("black",3))

col_strip = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", 
              alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5), 'gray70')
col_text = c(rep("white",5), rep("black",4))

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_strip[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$col <- col_text[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- col_text[i]
  
}

plot(g)


#..........................................................................
# Abundance data

load('VPCP_ab_bd.RDATA')
load('VPCP_ab_bf.RDATA')
load('VPCP_ab_rod.RDATA')
load('VPCP_ab_wg.RDATA')
load('VPCP_ab_moth.RDATA')

## load all linear predictor correlation
#load('FRcorr_taxa1010.RDATA')

# example dimensions
dim(VP.test3_bf$Vnorm) #    1  /   predictors group  /  species /  mcmc samples
#dim(FRcor_bf) #    predictors group (+ effort in our case)  /  predictors group (+ effort in our case)  /  species /  mcmc samples
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

## Posterior differences between M and V
VP.test3_bd$Vmarg[,,1,1]
VP.test3_bd$Vnorm[,,1,1]

mat.diff.bd = array(data=NA, dim = dim(VP.test3_bd$Vmarg))
for (i in 1:dim(VP.test3_bd$Vmarg)[4]){
  for (j in 1:dim(VP.test3_bd$Vmarg)[3]) {
    
    mat.bd =   sweep(
      VP.test3_bd$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_bd$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.bd[abs(VP.test3_bd$Vmarg[,,j,i])<0.01] = -1#*sign(mat.bd[abs(VP.test3_bd$Vmarg[,,j,i])<0.01])
    
    mat.diff.bd[,,j,i] = mat.bd
  }
}
dim(mat.diff.bd)
dimnames(mat.diff.bd) = dimnames(VP.test3_bd$Vmarg)


mat.diff.bf = array(data=NA, dim = dim(VP.test3_bf$Vmarg))
for (i in 1:dim(VP.test3_bf$Vmarg)[4]){
  for (j in 1:dim(VP.test3_bf$Vmarg)[3]) {
    
    mat.bf =   sweep(
      VP.test3_bf$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_bf$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.bf[abs(VP.test3_bf$Vmarg[,,j,i])<0.01] = -1#*sign(mat.bf[abs(VP.test3_bf$Vmarg[,,j,i])<0.01])
    
    mat.diff.bf[,,j,i] = mat.bf    
  }
}
dimnames(mat.diff.bf) = dimnames(VP.test3_bf$Vmarg)



mat.diff.rod = array(data=NA, dim = dim(VP.test3_rod$Vmarg))
for (i in 1:dim(VP.test3_rod$Vmarg)[4]){
  for (j in 1:dim(VP.test3_rod$Vmarg)[3]) {
    
    mat.rod =   sweep(
      VP.test3_rod$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_rod$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.rod[abs(VP.test3_rod$Vmarg[,,j,i])<0.01] = -1#*sign(mat.rod[abs(VP.test3_rod$Vmarg[,,j,i])<0.01])
    
    mat.diff.rod[,,j,i] = mat.rod    
  }
}
dimnames(mat.diff.rod) = dimnames(VP.test3_rod$Vmarg)


mat.diff.moth = array(data=NA, dim = dim(VP.test3_moth$Vmarg))
for (i in 1:dim(VP.test3_moth$Vmarg)[4]){
  for (j in 1:dim(VP.test3_moth$Vmarg)[3]) {
    mat.moth =   sweep(
      VP.test3_moth$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_moth$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.moth[abs(VP.test3_moth$Vmarg[,,j,i])<0.01] = -1#*sign(mat.moth[abs(VP.test3_moth$Vmarg[,,j,i])<0.01])
    
    mat.diff.moth[,,j,i] = mat.moth    
  }
}
dimnames(mat.diff.moth) = dimnames(VP.test3_moth$Vmarg)


mat.diff.wg = array(data=NA, dim = dim(VP.test3_wg$Vmarg))
for (i in 1:dim(VP.test3_wg$Vmarg)[4]){
  for (j in 1:dim(VP.test3_wg$Vmarg)[3]) {
    
    mat.wg =   sweep(
      VP.test3_wg$Vmarg[,,j,i],
      MARGIN = 1,           # 2 = columns, 1 = rows
      STATS = VP.test3_wg$Vnorm[,,j,i],
      FUN = function(x, y) exp(sign(x)*log(abs(x)) - log(y))-1
    )
    
    # to prevent unrealistic values
    # if M is too small
    mat.wg[abs(VP.test3_wg$Vmarg[,,j,i])<0.01] = -1#*sign(mat.wg[abs(VP.test3_wg$Vmarg[,,j,i])<0.01])
    
    mat.diff.wg[,,j,i] = mat.wg    
  }
}
dimnames(mat.diff.wg) = dimnames(VP.test3_wg$Vmarg)


## posterior mean

DiffMV_bd = melt(apply(mat.diff.bd, c(1,2,3), mean)) 
colnames(DiffMV_bd)[1] = 'var'
DiffMV_bd$Taxa = 'birds'

DiffMV_bf = melt(apply(mat.diff.bf, c(1,2,3), mean)) 
colnames(DiffMV_bf)[1] = 'var'
DiffMV_bf$Taxa = 'butterflies'

DiffMV_moth = melt(apply(mat.diff.moth, c(1,2,3), mean)) 
colnames(DiffMV_moth)[1] = 'var'
DiffMV_moth$Taxa = 'moths'

DiffMV_rod = melt(apply(mat.diff.rod, c(1,2,3), mean)) 
colnames(DiffMV_rod)[1] = 'var'
DiffMV_rod$Taxa = 'small mammals'

DiffMV_wg = melt(apply(mat.diff.wg, c(1,2,3), mean)) 
colnames(DiffMV_wg)[1] = 'var'
DiffMV_wg$Taxa = 'large mammals'


## combine all and keep only the shared importance of interests

DiffMV = rbind(DiffMV_bd, DiffMV_bf, DiffMV_moth, DiffMV_rod, DiffMV_wg)
DiffMV = DiffMV[!(DiffMV$var == DiffMV$element) ,]
DiffMV = DiffMV[!(DiffMV$var %in% c('site', 'year', 'bg') | DiffMV$element %in% c('site', 'year', 'bg')) ,]

DiffMV$Taxa = factor(DiffMV$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
DiffMV$var = factor(DiffMV$var, levels = c("Climate", "HabitatConf","Habitatprop"),
                    labels = c("Climate", "LandscapeConf","Habitatcomp"))
DiffMV$element = factor(DiffMV$element, levels = c("Climate", "HabitatConf","Habitatprop"),
                        labels = c("Climate", "LandscapeConf","Habitatcomp"))


## plotting

summary(DiffMV)
dim(DiffMV) #3006  5

DiffMV2 = DiffMV[DiffMV$value < 1,] # now only very few excluded
summary(DiffMV2)
dim(DiffMV2) # 2937


# re do mean after excluding extreme values?
mat.diff.bd[mat.diff.bd>1] = NA
DiffMV_bd = melt(apply(mat.diff.bd, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_bd)[1] = 'var'
DiffMV_bd$Taxa = 'birds'

mat.diff.bf[mat.diff.bf>1] = NA
DiffMV_bf = melt(apply(mat.diff.bf, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_bf)[1] = 'var'
DiffMV_bf$Taxa = 'butterflies'

mat.diff.moth[mat.diff.moth>1] = NA
DiffMV_moth = melt(apply(mat.diff.moth, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_moth)[1] = 'var'
DiffMV_moth$Taxa = 'moths'

mat.diff.rod[mat.diff.rod>1] = NA
DiffMV_rod = melt(apply(mat.diff.rod, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_rod)[1] = 'var'
DiffMV_rod$Taxa = 'small mammals'

mat.diff.wg[mat.diff.wg>1] = NA
DiffMV_wg = melt(apply(mat.diff.wg, c(1,2,3), mean, na.rm=T)) 
colnames(DiffMV_wg)[1] = 'var'
DiffMV_wg$Taxa = 'large mammals'


## combine all and keep only the shared importance of interests
DiffMV = rbind(DiffMV_bd, DiffMV_bf, DiffMV_moth, DiffMV_rod, DiffMV_wg)
DiffMV = DiffMV[!(DiffMV$var == DiffMV$element) ,]
DiffMV = DiffMV[!(DiffMV$var %in% c('site', 'year', 'bg') | DiffMV$element %in% c('site', 'year', 'bg')) ,]

DiffMV$Taxa = factor(DiffMV$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
DiffMV$var = factor(DiffMV$var, levels = c("Climate", "HabitatConf","Habitatprop"),
                    labels = c("Climate", "LandscapeConf","Habitatcomp"))
DiffMV$element = factor(DiffMV$element, levels = c("Climate", "HabitatConf","Habitatprop"),
                        labels = c("Climate", "LandscapeConf","Habitatcomp"))


# Get the percentage
# 
confoun = na.omit(DiffMV[DiffMV$value < -0.99,])
confoun %>%
  group_by(var, element) %>%
  summarise(count = 100*length(value)/(102+57+319+8+15))

confoun = confoun %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)
confoun$tot = c(rep(102, 6), rep(57, 6), rep(319, 4), rep(15, 6)) # check for species if a driver1-driver2 confounding exist and adapt the numbers of rep (ex none for small mammals)
confoun$pct = 100*confoun$count/confoun$tot
confoun %>% 
  print(n = 100)

confoun$shared = 'confounded'

#.......................................................
indepdt = na.omit(DiffMV[DiffMV$value > -0.01 & DiffMV$value < 0.01,])
indepdt %>%
  group_by(var, element) %>%
  summarise(count = 100*length(value)/(102+57+319+8+15))

indepdt = indepdt %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)
indepdt$tot = c(rep(102, 5), rep(57, 4), rep(319, 4), rep(15, 2))  # check for species if a driver1-driver2 independence exist and adapt the numbers of rep (ex none for small mammals)
indepdt$pct = 100*indepdt$count/indepdt$tot
indepdt %>% 
  print(n = 100)
indepdt$shared = 'Independent'


#.......................................................

Suppress = na.omit(DiffMV[DiffMV$value > -0.99 & DiffMV$value < -0.01,])
Suppress = Suppress %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)
Suppress$tot = c(rep(102, 6), rep(57, 6), rep(319, 6), rep(8, 6), rep(15, 6))
Suppress$pct = 100*Suppress$count/Suppress$tot
Suppress %>% 
  print(n = 100)

Suppress$shared = 'Suppressing'

#.......................................................

reinf = na.omit(DiffMV[DiffMV$value > 0.01 & DiffMV$value < 1,])
reinf = reinf %>%
  group_by(Taxa, var, element) %>%
  summarise(count = length(value))%>% 
  print(n = 100)

reinf$tot = c(rep(102, 6), rep(57, 4), rep(319, 6), #rep(8, 6), 
              rep(15, 6))
reinf$pct = 100*reinf$count/reinf$tot
reinf %>% 
  print(n = 100)


reinf$shared = 'Reinforcing'


# altogether
Effects = rbind(reinf, indepdt, Suppress, confoun)
Effects$shared = factor(Effects$shared,
                        levels = c('Reinforcing', 'Independent', 'Suppressing', 'confounded'))


gt_table2 <- Effects %>%
  gt(groupname_col = "shared", rowname_col = 'Taxa', row_group_as_column = FALSE)%>%
  tab_spanner_delim(delim = "_") |>
  tab_header(
    title = "Abundance data",
  ) %>%
  fmt_number(
    columns = c(pct),
    decimals = 2
  ) %>%
  cols_label(
    var = "Variable 1",
    element = "Variable 2",
    count = "Species #",
    tot = "Total #",
    pct = "Pct (%)"
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(everything())
  ) %>% tab_options(
    data_row.padding = px(2),
    summary_row.padding = px(3), # A bit more padding for summaries
    row_group.padding = px(6)    # And even more for our groups
  ) |> 
  opt_stylize(style = 3, color = 'gray')

# Print the table
gt_table2 |> gtsave("AB_shared.png", expand = 10)



#.......................................................
#.......................................................

# Plots
# row-wise

# prepare data properly
nb.sp = c('birds'=102, 'butterflies'=57, 'moths'=319,
                   'small mammals'=8, 'large mammals' = 15)
totals_df <- tibble(
  var = c(rep("Climate", 5), rep("LandscapeConf", 5), rep("Habitatcomp",5)),
  Taxa = c(rep(levels(DiffMV$Taxa), 3)),
  tot = c(rep(nb.sp, 3))
)
group_taxa = levels(DiffMV$Taxa)


# Prepare different bining to highlight independent and confounding effects
DiffMV = na.omit(DiffMV)

df_plot_indptconf <- DiffMV[(DiffMV$value >= -0.01) & (DiffMV$value <= 0.01) | (DiffMV$value <= -0.99),] %>%
  mutate(
    vbin = floor(value / 0.01) * 0.01
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

# öarger bins for others
df_plot_shared <- DiffMV[!((DiffMV$value> -0.01) & (DiffMV$value< 0.01)) & !(DiffMV$value < -0.99) ,] %>%
  mutate(
    vbin = floor(value / 0.05) * 0.05
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

df_plot = rbind(df_plot_indptconf, df_plot_shared)

df_plot$Taxa = factor(df_plot$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
df_plot$var = factor(df_plot$var, levels = c("Climate", "LandscapeConf","Habitatcomp"))
df_plot$element = factor(df_plot$element, levels = c("Climate", "LandscapeConf","Habitatcomp"))


plot = ggplot(df_plot, aes(vbin, pct, fill = element)) + 
  geom_col(width = 0.05, position = 'identity') +  # Not completely correct around 1
  scale_fill_manual(values = c(alpha("#66C2A5",0.4), alpha("yellow3",0.4), alpha('orange3',0.4))) +
  ylab('Species proportion, %') + xlab('Relative difference of the marginal variance \n partition (M) and the variance partition (V)') +
  geom_vline(xintercept = c(0), colour='gray50', linetype = "dashed") + ## Independent
  geom_vline(xintercept = c(-1, -0.99, 1), colour='gray50', linetype = "dashed") + ## confounded
  scale_x_continuous(
    limits = c(-1.025, 1.025),
    breaks = seq(-1, 1, by = 0.5),
    expand = c(0.05, 0)
  ) + 
  facet_grid(var ~ Taxa) +
  coord_flip(ylim = c(0, 40), clip = "on") +
  theme_classic() + theme(text = element_text(size = 18),
                          axis.text.x = element_text(angle = 45, hjust = 1),
                          axis.title.y = element_text(hjust = 0),
                          legend.position = "bottom",
                          panel.background = element_rect(fill = NA, color = "black")) +
  guides(fill=guide_legend(title=" "))



g <- ggplot_gtable(ggplot_build(plot))
strips <- which(grepl('strip-', g$layout$name))

#col_strip = c(rep("white",5), alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5))
#col_text = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", rep("black",3))

col_strip = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5))
col_text = c(rep("white",5), rep("black",3))

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_strip[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$col <- col_text[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- col_text[i]
  
}

plot(g)



# with random effects
#

## combine all and keep only the shared importance of interests
DiffMV = rbind(DiffMV_bd, DiffMV_bf, DiffMV_moth, DiffMV_rod, DiffMV_wg)
DiffMV = DiffMV[!(DiffMV$var == DiffMV$element) ,]

DiffMV$Taxa = factor(DiffMV$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))

DiffMV <- na.omit(DiffMV) %>%
  mutate(element = ifelse(element %in% c('year', 'site', 'bg'), "random", element)) %>%
  group_by(Taxa, var, element, Species) %>%
  summarise(value = mean(value, na.rm=T), .groups = "drop")

DiffMV <- na.omit(DiffMV) %>%
  mutate(var = ifelse(var %in% c('year', 'site', 'bg'), "random", var)) %>%
  group_by(Taxa, var, element, Species) %>%
  summarise(value = mean(value, na.rm=T), .groups = "drop")

DiffMV$var = factor(DiffMV$var, 
                    levels = c("1", "2", "3", 'random'),
                    labels = c("Climate", "LandscapeConf", "Habitatcomp",  "Random"))

DiffMV$element = factor(DiffMV$element, levels = c("1", "2", "3", 'random'),
                        labels = c("Climate", "LandscapeConf", "Habitatcomp",  "Random"))

# Since we have multiple random effects here let's summarize it for simplicity


totals_df <- tibble(
  var = c(rep("Climate", 5), rep("LandscapeConf", 5), rep("Habitatcomp",5), rep("Random",5)),
  Taxa = c(rep(levels(DiffMV$Taxa), 4)),
  tot = c(rep(nb.sp, 4))
)
group_taxa = levels(DiffMV$Taxa)


# Prepare different bining to highlight independent and confounding effects
DiffMV = na.omit(DiffMV)

df_plot_indptconf <- DiffMV[(DiffMV$value >= -0.01) & (DiffMV$value <= 0.01) | (DiffMV$value <= -0.99),] %>%
  mutate(
    vbin = floor(value / 0.01) * 0.01
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

# öarger bins for others
df_plot_shared <- DiffMV[!((DiffMV$value> -0.01) & (DiffMV$value< 0.01)) & !(DiffMV$value < -0.99) ,] %>%
  mutate(
    vbin = floor(value / 0.05) * 0.05
  ) %>%
  count(Taxa, element, var, vbin) %>%
  left_join(totals_df, by = c("var", "Taxa")) %>%
  mutate(
    pct = 100*n / tot
  )

df_plot = rbind(df_plot_indptconf, df_plot_shared)


df_plot$Taxa = factor(df_plot$Taxa, levels = c('birds', 'butterflies', 'moths', 'small mammals','large mammals'))
df_plot$var = factor(df_plot$var, levels = c("Climate", "LandscapeConf","Habitatcomp", "Random"))
df_plot$element = factor(df_plot$element, levels = c("Climate", "LandscapeConf","Habitatcomp", "Random"))



# row-wise


plot =  ggplot(df_plot, aes(vbin, pct, fill = element)) + 
  geom_col(width = 0.05, position = 'identity') +  
  scale_fill_manual(values = c(alpha("#66C2A5",0.4), alpha("yellow3",0.4), 
                               alpha('orange3',0.4), alpha('gray60',0.4))) +
  ylab('Species proportion, %') + xlab('Relative difference of the marginal variance \n partition (M) and the variance partition (V)') +
  geom_vline(xintercept = c(0), colour='gray50', linetype = "dashed") + ## Independent
  geom_vline(xintercept = c(-1, -0.99, 1), colour='gray50', linetype = "dashed") + ## confounded
  scale_x_continuous(
    limits = c(-1.025, 1.025),
    breaks = seq(-1, 1, by = 0.5),
    expand = c(0.05, 0)
  ) +  
  facet_grid(var ~ Taxa) +
  coord_flip(ylim = c(0, 80), clip = "on") +
  theme_classic() + theme(text = element_text(size = 18),
                          axis.text.x = element_text(angle = 45, hjust = 1),
                          legend.position = "bottom",
                          axis.title.y = element_text(hjust = 0),
                          panel.background = element_rect(fill = NA, color = "black")) +
  guides(fill=guide_legend(title=" "))




g <- ggplot_gtable(ggplot_build(plot))
strips <- which(grepl('strip-', g$layout$name))

#col_strip = c(rep("white",5), alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5))
#col_text = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", rep("black",3))

col_strip = c("#C1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310", 
              alpha("#66C2A5",0.5), alpha("yellow3",0.5), alpha('orange3',0.5), 'gray70')
col_text = c(rep("white",5), rep("black",4))

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_strip[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$col <- col_text[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- col_text[i]
  
}

plot(g)
