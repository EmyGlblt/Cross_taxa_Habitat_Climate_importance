library(dplyr)
library(grid)
library(gtable)
library(ggplot2)
library(ggridges)
library(geomtextpath)

# for plotting later
setwd('C:/Users/guilbaul/OneDrive - University of Helsinki/RECcII/ScriptData/Code/x.New_functions')
source("function_ggplot_strip_color.R")

setwd('C:/Users/guilbaul/OneDrive - University of Helsinki/RECcII/ScriptData/Data/4.VPanalyses/Joint_interp')

## Occurrence data .................................................
load('VP_C1resc_all_occ1111.RDATA')

## Combine plotting from the different posterior probability

Newcol = c("#D1A10A", "#05662C", "#B5BF99", "#9D6058", "#762310")

colnames(VP_C_all_occ)
head(VP_C_all_occ)


VP_C_all_occ$prop = as.numeric(VP_C_all_occ$prop)
table(VP_C_all_occ$value)

## 1 confounded 
# 3 independent

pal = c('#66C2A5', 'yellow3', 'orange3', 'gray60')


# prepare the plotting
rc_un = data.frame(col_name= c("Climate", "LandscapeConf", "HabitatComp", 'Random'),
                   row_name=c("Climate", "LandscapeConf", "HabitatComp", 'Random'))
rc_un$row_name = factor(rc_un$row_name, 
                        levels = c("Climate", "LandscapeConf","HabitatComp", 'Random'))
rc_un$col_name = factor(rc_un$col_name, 
                        levels = c("Climate", "LandscapeConf", "HabitatComp", 'Random'))


VP_C_all_occ.sub = VP_C_all_occ[VP_C_all_occ$col_name %in% c("Climate", "LandscapeConf","Habitatcomp", 'Random') &
                                  VP_C_all_occ$row_name %in% c("Climate", "LandscapeConf","Habitatcomp", 'Random'),]

VP_C_all_occ.sub$col_name = droplevels(VP_C_all_occ.sub$col_name)
VP_C_all_occ.sub$row_name = droplevels(VP_C_all_occ.sub$row_name)

## cleaning up names
VP_C_all_occ.sub$value[VP_C_all_occ.sub$value == 'Amplify'] = 'Reinforce'

levels(VP_C_all_occ.sub$col_name)[levels(VP_C_all_occ.sub$col_name)== 'Habitatcomp'] = 'HabitatComp'
levels(VP_C_all_occ.sub$row_name)[levels(VP_C_all_occ.sub$row_name)== 'Habitatcomp'] = 'HabitatComp'


## Abundance data .................................................
load('VP_Cresc_all_ab1111.RDATA')

colnames(VP_C_all_ab)
head(VP_C_all_ab)


VP_C_all_ab$prop = as.numeric(VP_C_all_ab$prop)
table(VP_C_all_ab$value)


VP_C_all_ab.sub = VP_C_all_ab[VP_C_all_ab$col_name %in% c("Climate", "LandscapeConf","Habitatcomp", 'Random') &
                                VP_C_all_ab$row_name %in% c("Climate", "LandscapeConf","Habitatcomp", 'Random'),]

VP_C_all_ab.sub$col_name = droplevels(VP_C_all_ab.sub$col_name)
VP_C_all_ab.sub$row_name = droplevels(VP_C_all_ab.sub$row_name)


# cleaning up names
VP_C_all_ab.sub$value[VP_C_all_ab.sub$value == 'Amplify'] = 'Reinforce'

levels(VP_C_all_ab.sub$col_name)[levels(VP_C_all_ab.sub$col_name)== 'Habitatcomp'] = 'HabitatComp'
levels(VP_C_all_ab.sub$row_name)[levels(VP_C_all_ab.sub$row_name)== 'Habitatcomp'] = 'HabitatComp'

## ...................................................................................
## ...................................................................................
## ...................................................................................
## comparison


# Occurrence
VP_C_mean_occ = VP_C_all_occ.sub %>%
  group_by(col_name, row_name, NBsp, name, value, taxa) %>%
  dplyr::summarise(Nmean = mean(n), Propmean = mean(prop),
                   Propup = quantile(prop, probs= 0.975),
                   Proplow = quantile(prop, probs= 0.025))

VP_C_mean_occ$Propmean = as.numeric(VP_C_mean_occ$Propmean)
levels(VP_C_mean_occ$taxa) = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals')

gg_old = {ggplot() +
    geom_rect(data = rc_un, aes(fill = row_name), fill='black', 
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5)+ 
    geom_point(data = VP_C_mean_occ, aes(y=value, x=Propmean, 
                                                   colour=taxa)) + 
    scale_fill_manual(values = c('birds' = "#D1A10A", 'butterflies' = "#05662C", 'moths' = "#B5BF99", 
                                 'small mammals' = "#9D6058", 'large mammals' = "#762310")) +
    facet_grid(row_name ~ col_name, drop = FALSE, switch = "y")+
    theme(panel.grid.major.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.margin=unit(c(0,0,0,0), "cm"))}%>%
  modify_facet_appearance(strip.background.x.fill = pal,
                          strip.background.x.col = rep('black', 4),
                          strip.text.x.col = rep('black', 4),
                          strip.background.y.fill = pal,
                          strip.background.y.col =rep('black', 4),
                          strip.text.y.col = rep('black', 4))

gp = {ggplot(data = VP_C_mean_occ, aes(y=Propmean, x=value, fill=taxa)) +
    geom_bar(stat = 'identity', position = "dodge") +
    geom_errorbar(data = VP_C_mean_occ, aes(x=value, ymin=Proplow, ymax=Propup, colour=taxa), 
                  position = "dodge", alpha=0.9, size=0.8)+
    scale_fill_manual(values = c('birds' = "#D1A10A", 'butterflies' = "#05662C", 'moths' = "#B5BF99", 
                                 'small mammals' = "#9D6058", 'large mammals' = "#762310")) +
    scale_colour_manual(values = c('birds' = "#D1A10A", 'butterflies' = "#05662C", 'moths' = "#B5BF99", 
                                   'small mammals' = "#9D6058", 'large mammals' = "#762310")) +
    facet_grid(row_name ~ col_name, drop = FALSE, switch = "y") + 
    coord_curvedpolar(start = -pi/12) +
    geom_hline(yintercept = seq(0, 0.75, by = 0.25), colour = "darkgrey", size = 0.55) +
    geom_hline(yintercept = 0.25, colour = "gray40", size = 0.75) +
    geom_hline(yintercept = 0.75, colour = "gray40", size = 0.75) +
    geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5), colour = "gray40", size = 0.85) +
    scale_y_continuous(breaks = seq(0, 0.75, by = 0.25)) +
    labs(y= "% of species", x = "") +
    theme(panel.grid.major.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(colour = c("darkgrey", "gray40", "darkgrey", "gray40")),
          axis.text.x = element_text(vjust = 0.65),
          panel.ontop =FALSE,
          plot.margin=unit(c(0,0,0,0), "cm"))+
    guides(color='none',
           fill=guide_legend(title="Taxa", override.aes=list(size=4)))}%>%
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


par(mar = c(0,0,0,0))

## plot again
grid.newpage()
png("JointInterSupp_sc_occ.png") 
grid.draw(gp) 
dev.off()


 
 ## Abundance


VP_C_mean_ab = VP_C_all_ab.sub %>%
  group_by(col_name, row_name, NBsp, name, value, taxa) %>%
  dplyr::summarise(Nmean = mean(n), Propmean = mean(prop),
                   Propup = quantile(prop, probs= 0.975),
                   Proplow = quantile(prop, probs= 0.025))


VP_C_mean_ab$Propmean = as.numeric(VP_C_mean_ab$Propmean)
levels(VP_C_mean_ab$taxa) = c('birds', 'butterflies', 'moths', 'small mammals', 'large mammals')

gg_old = {ggplot() +
    geom_rect(data = rc_un, aes(fill = row_name), fill='black', 
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5)+ 
    geom_point(data = VP_C_mean_ab, aes(y=Propmean, x=value, colour=taxa)) +
    scale_colour_manual(values = c('birds' = "#D1A10A", 'butterflies' = "#05662C", 'moths' = "#B5BF99", 
                                   'small mammals' = "#9D6058", 'large mammals' = "#762310")) +
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

gp_ab = {ggplot(data = VP_C_mean_ab, aes(y=Propmean, x=value, fill=taxa)) +
    geom_bar(stat = 'identity', position = "dodge") +
    geom_errorbar(data = VP_C_mean_ab, aes(x=value, ymin=Proplow, ymax=Propup, colour=taxa), 
                  position = "dodge", alpha=0.9, size=0.8)+
    scale_fill_manual(values = c('birds' = "#D1A10A", 'butterflies' = "#05662C", 'moths' = "#B5BF99", 
                                 'small mammals' = "#9D6058", 'large mammals' = "#762310")) +
    scale_colour_manual(values = c('birds' = "#D1A10A", 'butterflies' = "#05662C", 'moths' = "#B5BF99", 
                                   'small mammals' = "#9D6058", 'large mammals' = "#762310")) +
    facet_grid(row_name ~ col_name, drop = FALSE, switch = "y") + 
    coord_curvedpolar(start = -pi/12) +
    geom_hline(yintercept = seq(0, 0.75, by = 0.25), colour = "darkgrey", size = 0.55) +
    geom_hline(yintercept = 0.25, colour = "gray40", size = 0.75) +
    geom_hline(yintercept = 0.75, colour = "gray40", size = 0.75) +
    geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5), colour = "gray40", size = 0.85) +
    scale_y_continuous(breaks = seq(0, 0.75, by = 0.25)) +
    labs(y= "% of species", x = "") +
    theme(panel.grid.major.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(colour = c("darkgrey", "gray40", "darkgrey", "gray40")),
          axis.text.x = element_text(vjust = 0.65),
          panel.ontop =FALSE)+
    guides(color='none',
           fill=guide_legend(title="Taxa", override.aes=list(size=4)))}%>%
  modify_facet_appearance(strip.background.x.fill = pal,
                          strip.background.x.col = rep('black', 4),
                          strip.text.x.col = rep('black', 4),
                          strip.background.y.fill = pal,
                          strip.background.y.col =rep('black', 4),
                          strip.text.y.col = rep('black', 4))

plot(gp_ab)

# Remove facets
idx <- which(gp_ab$layout$name %in% c("panel-1-1", "panel-2-2", "panel-3-3", "panel-4-4"))
for (i in idx) gp_ab$grobs[[i]] <- gg_old$grobs[[i]]

## plot again
grid.newpage()
grid.draw(gp_ab)

png("JointInterSupp_sc_ab.png") 
grid.draw(gp_ab) 
dev.off()
