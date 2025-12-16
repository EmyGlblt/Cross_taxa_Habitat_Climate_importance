#
#             Summary plots to compare taxa
#
#
# --------------------------------------------------------------------



library(ggtern)
library(Hmsc)
library(reshape2) # necessary for margina cal-> dcast function used
library(corrplot)
library(ggExtra)
library(cowplot)
library(patchwork)


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(grid)


setwd("D:/Helsinki/RECcII/Results_VP_Sept_2025")

load('VPtrait_ab_rod.RDATA')
load('VPtrait_ab_bf.RDATA')
load('VPtrait_ab_bird.RDATA')
load('VPtrait_ab_wg.RDATA')
load('VPtrait_ab_moth.RDATA')

BdDat_trait_ab$taxa = 'Bird'
BfDat_trait_ab$taxa = 'Butterfly'
mothDat_trait_ab$taxa = 'Moth'
rodDat_trait_ab$taxa = 'Rodent'
wgDat_trait_ab$taxa = 'Winter Game'

#.......................................................  Bird
BdDat_trait_ab2 = BdDat_trait_ab[,c(1:4, 10, 11, 29, 34)]
summary(BdDat_trait_ab2)


### cat - body length
hist(BdDat_trait_ab2$WingLengthMean)
BdDat_trait_ab2$BodySize = NA
BdDat_trait_ab2$BodySize[BdDat_trait_ab2$WingLengthMean <= 108] = 'small'
BdDat_trait_ab2$BodySize[BdDat_trait_ab2$WingLengthMean > 108] = 'Large'

table(BdDat_trait_ab2$BodySize)


### cat - generation 
BdDat_trait_ab2$Nr_Broods_S
BdDat_trait_ab2$Lpyear = NA
BdDat_trait_ab2$Lpyear[BdDat_trait_ab2$Nr_Broods_S == 1] = 'slow'
BdDat_trait_ab2$Lpyear[BdDat_trait_ab2$Nr_Broods_S == 2] = 'fast'

table(BdDat_trait_ab2$Lpyear)

### cat - Hab breadth
table(BdDat_trait_ab2$Habitat.specialism)

BdDat_trait_ab2$Habspe = NA
BdDat_trait_ab2$Habspe[BdDat_trait_ab2$Habitat.specialism > 0.5] = 'specialist'
BdDat_trait_ab2$Habspe[BdDat_trait_ab2$Habitat.specialism <= 0.5] = 'generalist'

table(BdDat_trait_ab2$Habspe)


# .......................................................  Moth
mothDat_trait_ab2 = mothDat_trait_ab[,c(1:4, 5, 8, 6, 9)]
summary(mothDat_trait_ab2)

### cat - body length
hist(mothDat_trait_ab2$Body.size)
mothDat_trait_ab2$BodySize = NA

mothDat_trait_ab2$BodySize[mothDat_trait_ab2$Body.size <= 32] = 'small'
mothDat_trait_ab2$BodySize[mothDat_trait_ab2$Body.size > 32] = 'Large'

table(mothDat_trait_ab2$BodySize)


### cat - generation per year
hist(mothDat_trait_ab2$Voltinism)
mothDat_trait_ab2$Lpyear = NA
mothDat_trait_ab2$Lpyear[mothDat_trait_ab2$Voltinism %in% c('0', '1')] = 'slow'
mothDat_trait_ab2$Lpyear[mothDat_trait_ab2$Voltinism %in% c('3', '2')] = 'fast'

table(mothDat_trait_ab2$Lpyear)


### cat - Hab breadth
table(mothDat_trait_ab2$Diet)

mothDat_trait_ab2$Habspe = NA
mothDat_trait_ab2$Habspe[mothDat_trait_ab2$Diet %in% c('monophag', 'oligophag', 'lichen', 'mushroom')] = 'specialist'
mothDat_trait_ab2$Habspe[mothDat_trait_ab2$Diet %in% 'polyphag'] = 'generalist'

table(mothDat_trait_ab2$Habspe)


# ....................................................... butterfly
BfDat_trait_ab2 = BfDat_trait_ab
summary(BfDat_trait_ab2)


### cat - body length
hist(BfDat_trait_ab2$Body.size)
BfDat_trait_ab2$BodySize = NA
BfDat_trait_ab2$BodySize[BfDat_trait_ab2$Body.size <= 37] = 'small'
BfDat_trait_ab2$BodySize[BfDat_trait_ab2$Body.size > 37] = 'Large'

table(BfDat_trait_ab2$BodySize)


### cat - generation per year
hist(BfDat_trait_ab2$Voltinism)
BfDat_trait_ab2$Lpyear = NA
BfDat_trait_ab2$Lpyear[BfDat_trait_ab2$Voltinism %in% c('0', '1')] = 'slow'
BfDat_trait_ab2$Lpyear[BfDat_trait_ab2$Voltinism %in% c('3', '2')] = 'fast'

table(BfDat_trait_ab2$Lpyear)


### cat - Hab breadth
table(BfDat_trait_ab2$Breadth.of.habitat.use.adult.butterflies)

BfDat_trait_ab2$Habspe = NA
BfDat_trait_ab2$Habspe[BfDat_trait_ab2$Specificity.of.larval.host.plant.use %in% c('monophagous', 'oligophagous')] = 'specialist'
BfDat_trait_ab2$Habspe[BfDat_trait_ab2$Specificity.of.larval.host.plant.use %in% 'polyphagous'] = 'generalist'
table(BfDat_trait_ab2$Habspe)

BfDat_trait_ab2 = BfDat_trait_ab2[, c(1:5, 8, 6, 9:12)]

#.......................................................  rodents
rodDat_trait_ab2 = rodDat_trait_ab[,c(1:4, 6, 8, 9, 12)]
summary(rodDat_trait_ab2)

## If group cat

### cat - body length
hist(rodDat_trait_ab2$adult_body_length_mm)
rodDat_trait_ab2$BodySize = NA
rodDat_trait_ab2$BodySize[rodDat_trait_ab2$adult_body_length_mm  <= 100] = 'small'
rodDat_trait_ab2$BodySize[rodDat_trait_ab2$adult_body_length_mm  > 100] = 'Large'

table(rodDat_trait_ab2$BodySize )


### cat - Litter per year
table(rodDat_trait_ab2$litters_per_year_n)
rodDat_trait_ab2$Lpyear = NA
rodDat_trait_ab2$Lpyear[rodDat_trait_ab2$litters_per_year_n <= 2.8] = 'slow'
rodDat_trait_ab2$Lpyear[rodDat_trait_ab2$litters_per_year_n > 2.8] = 'fast'

table(rodDat_trait_ab2$Lpyear)


### cat - Hab breadth
table(rodDat_trait_ab2$habitat_breadth_n)

rodDat_trait_ab2$Habspe = NA
rodDat_trait_ab2$Habspe[rodDat_trait_ab2$habitat_breadth_n %in% c('2')] = 'specialist'
rodDat_trait_ab2$Habspe[rodDat_trait_ab2$habitat_breadth_n %in% c('5','6','3')] = 'generalist'

table(rodDat_trait_ab2$Habspe)


#....................................................... Winter game triangle
wgDat_trait_ab2 = wgDat_trait_ab[,c(1:4, 6, 8, 9, 12)]
summary(wgDat_trait_ab2)

### cat - body length
hist(wgDat_trait_ab2$adult_body_length_mm)
wgDat_trait_ab2$BodySize = NA
wgDat_trait_ab2$BodySize[wgDat_trait_ab2$adult_body_length_mm  <= 610] = 'small'
wgDat_trait_ab2$BodySize[wgDat_trait_ab2$adult_body_length_mm  > 610] = 'Large'

table(wgDat_trait_ab2$BodySize )


### cat - Litter per year
hist(wgDat_trait_ab2$litters_per_year_n)
wgDat_trait_ab2$Lpyear = NA
wgDat_trait_ab2$Lpyear[wgDat_trait_ab2$litters_per_year_n <= 1] = 'slow'
wgDat_trait_ab2$Lpyear[wgDat_trait_ab2$litters_per_year_n > 1] = 'fast'

table(wgDat_trait_ab2$Lpyear)


### cat - Hab breadth
table(wgDat_trait_ab2$habitat_breadth_n)

wgDat_trait_ab2$Habspe = NA
wgDat_trait_ab2$Habspe[wgDat_trait_ab2$habitat_breadth_n <= 2] = 'specialist'
wgDat_trait_ab2$Habspe[wgDat_trait_ab2$habitat_breadth_n > 2] = 'generalist'

table(wgDat_trait_ab2$Habspe)


#.......................................................
colnames(BdDat_trait_ab2)[1:8] = c('species', 'Climate', 'habConf', 'habprop', 'BodySizeCont', 'FSlifehistCont', 'HabspeCont', 'Taxa')
colnames(BfDat_trait_ab2)[1:8] = c('species', 'Climate', 'habConf', 'habprop', 'BodySizeCont', 'FSlifehistCont', 'HabspeCont', 'Taxa')
colnames(mothDat_trait_ab2)[1:8] = c('species', 'Climate', 'habConf', 'habprop', 'BodySizeCont', 'FSlifehistCont', 'HabspeCont', 'Taxa')
colnames(rodDat_trait_ab2)[1:8] = c('species', 'Climate', 'habConf', 'habprop', 'BodySizeCont', 'FSlifehistCont', 'HabspeCont', 'Taxa')
colnames(wgDat_trait_ab2)[1:8] = c('species', 'Climate', 'habConf', 'habprop', 'BodySizeCont', 'FSlifehistCont', 'HabspeCont', 'Taxa')

Dat_trait_all = rbind(BdDat_trait_ab2, BfDat_trait_ab2, mothDat_trait_ab2, rodDat_trait_ab2, wgDat_trait_ab2)
head(Dat_trait_all)
colnames(Dat_trait_all)
summary(Dat_trait_all)


Dat_trait_all2 = Dat_trait_all[Dat_trait_all$Climate<100 & Dat_trait_all$Climate>0,]
Dat_trait_all3 = Dat_trait_all2[Dat_trait_all2$habprop<100 & Dat_trait_all2$habprop>0,]
Dat_trait_all = Dat_trait_all3[Dat_trait_all3$habConf<100 & Dat_trait_all3$habConf>0,]

Dat_trait_all$sum = Dat_trait_all$Climate + Dat_trait_all$habprop + Dat_trait_all$habConf

Dat_trait_all$Climate = Dat_trait_all$Climate*100/Dat_trait_all$sum
Dat_trait_all$habprop = Dat_trait_all$habprop*100/Dat_trait_all$sum
Dat_trait_all$habConf = Dat_trait_all$habConf*100/Dat_trait_all$sum



## tern plots ..................

summary(Dat_trait_all)

Dat_trait_all = Dat_trait_all[!is.na(Dat_trait_all$BodySizeCont) ,]
Dat_trait_all = Dat_trait_all[!is.na(Dat_trait_all$FSlifehistCont) ,]

VP.birds = Dat_trait_all[Dat_trait_all$Taxa == 'Bird',]

VP.birdsHLmean = VP.birds %>%
  group_by(Lpyear) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.birdsHLmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternBirdHL = ggtern(VP.birds, aes(habprop, Climate, habConf, fill = Lpyear, color=Lpyear)) +
  stat_density_tern(aes(alpha = ..level.., fill = Lpyear), geom = "polygon", 
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_color_manual(values = c('#99CC96','#009969')) +
  geom_crosshair_tern(data = VP.birdsHLmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.birdsHLmean, aes(color = Lpyear), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternBirdHL

ggdensity(VP.birds, 'Climate', fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
#ggsave("climglobal.png")

ggdensity(VP.birds, "habprop",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
#ggsave("habPglobal.png")

ggdensity(VP.birds, "habConf",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

#ggarrange(habPplot, habCplot, climplot,
 #         ncol = 1)



VP.birdsHSmean = VP.birds %>%
  group_by(Habspe) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.birdsHSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternBirdHS = ggtern(VP.birds, aes(habprop, Climate, habConf, fill = Habspe, color=Habspe)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = Habspe), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values = c('#FFCC83','#795653')) +
  geom_crosshair_tern(data = VP.birdsHSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.birdsHSmean, aes(color = Habspe), show.legend = F) + 
  #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
  #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
  #theme_showarrows() +
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternBirdHS

ggdensity(VP.birds, 'Climate', fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
#ggsave("climglobal.png")

ggdensity(VP.birds, "habprop",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
#ggsave("habPglobal.png")

ggdensity(VP.birds, "habConf",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

#ggarrange(habPplot, habCplot, climplot,
#          ncol = 1)



VP.birdsBSmean = VP.birds %>%
  group_by(BodySize) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.birdsBSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternBirdBS = ggtern(VP.birds, aes(habprop, Climate, habConf, fill = BodySize, color=BodySize)) +
  #geom_point()+
  stat_density_tern(aes(alpha = ..level.., fill = BodySize), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  #position = "identity", bins = 10, h=1.5, show.legend = F, base = 'ilr') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  #geom_point(aes(color = taxa), shape = 4, size = 0.4) +
  scale_color_manual(values =c('#FF35FF','#665099')) +
  geom_crosshair_tern(data = VP.birdsBSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.birdsBSmean, aes(color = BodySize), show.legend = F) + 
  #scale_color_manual(values = c("#D8B90D", "#01401D", "#A1B142", "#41A66D", "#672A25")) +
  #labs(title  = "VP taxa", Larrow = "% habitat prop", Tarrow = "% Climate", Rarrow = "% habitat configuration") +
  #theme_showarrows() +
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternBirdBS

ggdensity(VP.birds, 'Climate', fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 
#ggsave("climglobal.png")

ggdensity(VP.birds, "habprop",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 
#ggsave("habPglobal.png")

ggdensity(VP.birds, "habConf",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

#ggarrange(habPplot, habCplot, climplot,
#          ncol = 1)



### ......................................................
VP.buttfly = Dat_trait_all[Dat_trait_all$Taxa == 'Butterfly',]

VP.buttflyHLmean = VP.buttfly %>%
  group_by(Lpyear) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.buttflyHLmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternbuttflyHL = ggtern(VP.buttfly, aes(habprop, Climate, habConf, fill = Lpyear, color=Lpyear)) +
  stat_density_tern(aes(alpha = ..level.., fill = Lpyear), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_color_manual(values = c('#99CC96','#009969')) +
  geom_crosshair_tern(data = VP.buttflyHLmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.buttflyHLmean, aes(color = Lpyear), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternbuttflyHL

ggdensity(VP.buttfly, 'Climate', fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.buttfly, "habprop",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.buttfly, "habConf",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.buttflyHSmean = VP.buttfly %>%
  group_by(Habspe) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.buttflyHSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternbuttflyHS = ggtern(VP.buttfly, aes(habprop, Climate, habConf, fill = Habspe, color=Habspe)) +
  stat_density_tern(aes(alpha = ..level.., fill = Habspe), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_color_manual(values = c('#FFCC83','#795653')) +
  geom_crosshair_tern(data = VP.buttflyHSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.buttflyHSmean, aes(color = Habspe), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternbuttflyHS

ggdensity(VP.buttfly, 'Climate', fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.buttfly, "habprop",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.buttfly, "habConf",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.buttflyBSmean = VP.buttfly %>%
  group_by(BodySize) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.buttflyBSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternbuttflyBS = ggtern(VP.buttfly, aes(habprop, Climate, habConf, fill = BodySize, color=BodySize)) +
  stat_density_tern(aes(alpha = ..level.., fill = BodySize), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_color_manual(values =c('#FF35FF','#665099')) +
  geom_crosshair_tern(data = VP.buttflyBSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.buttflyBSmean, aes(color = BodySize), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternbuttflyBS

ggdensity(VP.buttfly, 'Climate', fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.buttfly, "habprop",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.buttfly, "habConf",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 




### ......................................................
VP.Moth = Dat_trait_all[Dat_trait_all$Taxa == 'Moth',]

VP.Moth[is.na(VP.Moth$Climate),]

VP.MothHLmean = VP.Moth %>%
  group_by(Lpyear) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.MothHLmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternMothHL = ggtern(VP.Moth, aes(habprop, Climate, habConf, fill = Lpyear, color=Lpyear)) +
  stat_density_tern(aes(alpha = ..level.., fill = Lpyear), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_color_manual(values = c('#99CC96','#009969')) +
  geom_crosshair_tern(data = VP.MothHLmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.MothHLmean, aes(color = Lpyear), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none',
         fill=guide_legend(title="Life pace")) + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternMothHL

ggdensity(VP.Moth, 'Climate', fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.Moth, "habprop",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.Moth, "habConf",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.MothHSmean = VP.Moth %>%
  group_by(Habspe) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.MothHSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternMothHS = ggtern(VP.Moth, aes(habprop, Climate, habConf, fill = Habspe, color=Habspe)) +
  stat_density_tern(aes(alpha = ..level.., fill = Habspe), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_color_manual(values = c('#FFCC83','#795653')) +
  geom_crosshair_tern(data = VP.MothHSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.MothHSmean, aes(color = Habspe), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none',
         fill=guide_legend(title="habitat/diet")) + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternMothHS

ggdensity(VP.Moth, 'Climate', fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.Moth, "habprop",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.Moth, "habConf",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.MothBSmean = VP.Moth %>%
  group_by(BodySize) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.MothBSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternMothBS = ggtern(VP.Moth, aes(habprop, Climate, habConf, fill = BodySize, color=BodySize)) +
  stat_density_tern(aes(alpha = ..level.., fill = BodySize), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_color_manual(values =c('#FF35FF','#665099')) +
  geom_crosshair_tern(data = VP.MothBSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.MothBSmean, aes(color = BodySize), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none',
         fill=guide_legend(title="Body size")) + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternMothBS

ggdensity(VP.Moth, 'Climate', fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.Moth, "habprop",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.Moth, "habConf",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 



ggternBBM = ggtern::grid.arrange(ternBirdHL, ternbuttflyHL, ternMothHL,
          ternBirdHS, ternbuttflyHS, ternMothHS,
          ternBirdBS, ternbuttflyBS, ternMothBS,
          ncol=3, nrow=3)

ggsave(ggternBBM, file=paste0("ggternBBM_ab",".png", sep = ""))




### ......................................................
VP.Rodent = Dat_trait_all[Dat_trait_all$Taxa == 'Rodent',]

VP.RodentHLmean = VP.Rodent %>%
  group_by(Lpyear) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.RodentHLmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternRodentHL = ggtern(VP.Rodent, aes(habprop, Climate, habConf, fill = Lpyear, color=Lpyear)) +
  stat_density_tern(aes(alpha = ..level.., fill = Lpyear), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_color_manual(values = c('#99CC96','#009969')) +
  geom_crosshair_tern(data = VP.RodentHLmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.RodentHLmean, aes(color = Lpyear), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternRodentHL

ggdensity(VP.Rodent, 'Climate', fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.Rodent, "habprop",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.Rodent, "habConf",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.RodentHSmean = VP.Rodent %>%
  group_by(Habspe) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.RodentHSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternRodentHS = ggtern(VP.Rodent, aes(habprop, Climate, habConf, fill = Habspe, color=Habspe)) +
  stat_density_tern(aes(alpha = ..level.., fill = Habspe), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_color_manual(values = c('#FFCC83','#795653')) +
  geom_crosshair_tern(data = VP.RodentHSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.RodentHSmean, aes(color = Habspe), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternRodentHS

ggdensity(VP.Rodent, 'Climate', fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.Rodent, "habprop",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.Rodent, "habConf",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.RodentBSmean = VP.Rodent %>%
  group_by(BodySize) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.RodentBSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternRodentBS = ggtern(VP.Rodent, aes(habprop, Climate, habConf, fill = BodySize, color=BodySize)) +
  stat_density_tern(aes(alpha = ..level.., fill = BodySize), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = F, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_color_manual(values =c('#FF35FF','#665099')) +
  geom_crosshair_tern(data = VP.RodentBSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.RodentBSmean, aes(color = BodySize), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none') + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternRodentBS

ggdensity(VP.Rodent, 'Climate', fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.Rodent, "habprop",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.Rodent, "habConf",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 



### ......................................................
VP.WG = Dat_trait_all[Dat_trait_all$Taxa == 'Winter Game',]

VP.WGHLmean = VP.WG %>%
  group_by(Lpyear) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.WGHLmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternWGHL = ggtern(VP.WG, aes(habprop, Climate, habConf, fill = Lpyear, color=Lpyear)) +
  stat_density_tern(aes(alpha = ..level.., fill = Lpyear), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_color_manual(values = c('#99CC96','#009969')) +
  geom_crosshair_tern(data = VP.WGHLmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.WGHLmean, aes(color = Lpyear), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none',
         fill=guide_legend(title="Life pace")) + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternWGHL

ggdensity(VP.WG, 'Climate', fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.WG, "habprop",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.WG, "habConf",  fill = 'Lpyear') +
  scale_fill_manual(values = c('#99CC66','#009999')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.WGHSmean = VP.WG %>%
  group_by(Habspe) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.WGHSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternWGHS = ggtern(VP.WG, aes(habprop, Climate, habConf, fill = Habspe, color=Habspe)) +
  stat_density_tern(aes(alpha = ..level.., fill = Habspe), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_color_manual(values = c('#FFCC83','#795653')) +
  geom_crosshair_tern(data = VP.WGHSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.WGHSmean, aes(color = Habspe), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none',
         fill=guide_legend(title="Habitat/Diet")) + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternWGHS

ggdensity(VP.WG, 'Climate', fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.WG, "habprop",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,120) + guides(fill='none') + clean_theme() 

ggdensity(VP.WG, "habConf",  fill = 'Habspe') +
  scale_fill_manual(values = c('#FFCC33','#996633')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 


VP.WGBSmean = VP.WG %>%
  group_by(BodySize) %>%
  summarise(across(everything(), list(mean)))
colnames(VP.WGBSmean)[3:5] = c('Climate', 'habConf', 'habprop')

ternWGBS = ggtern(VP.WG, aes(habprop, Climate, habConf, fill = BodySize, color=BodySize)) +
  stat_density_tern(aes(alpha = ..level.., fill = BodySize), geom = "polygon", #bdl = c(0, 100),
                    position = "identity", bins = 10, h=0.25, show.legend = T, base = 'identity') +
  scale_alpha_continuous(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_color_manual(values =c('#FF35FF','#665099')) +
  geom_crosshair_tern(data = VP.WGBSmean, lty = 2, size = 2.1) +  # add mean information and lines
  geom_point(data = VP.WGBSmean, aes(color = BodySize), show.legend = F) + 
  theme(tern.axis.arrow.L = element_line(size=0, color='white'), 
        tern.axis.arrow.T = element_line(size=0, color="white"),
        tern.axis.arrow.R = element_line(size=0, color="white"),
        tern.axis.line.L = element_line(color='orange3',size=2),
        tern.axis.line.T = element_line(color='#66C2A5',size=2),
        tern.axis.line.R = element_line(color='yellow3',size=2),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17)) + 
  guides(alpha='none', colour='none',
         fill=guide_legend(title="Body size")) + 
  labs(title = '', x='', y='', z='', 
       xarrow = "",
       yarrow = "",
       zarrow = "")

ternWGBS

ggdensity(VP.WG, 'Climate', fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  scale_y_reverse()+ scale_x_reverse() +
  xlim(100,0) + guides(fill='none') + clean_theme() 

ggdensity(VP.WG, "habprop",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

ggdensity(VP.WG, "habConf",  fill = 'BodySize') +
  scale_fill_manual(values = c('#FF85FF','#660099')) +
  xlim(0,100) + guides(fill='none') + clean_theme() 

# empty tern plot to keep proportion in group plotting
tmp <- data.frame(x=1/3,
                  y=1/3,
                  z=1/3,
                  Min=1/3-1/6,
                  Max=1/3+1/6)
gg_blank = ggtern(data=tmp,aes(x,y,z)) +
  geom_point() +
  geom_errorbarT(aes(Tmin=Min,Tmax=Max))+
  geom_errorbarL(aes(Lmin=Min,Lmax=Max))+
  geom_errorbarR(aes(Rmin=Min,Rmax=Max))



ggternRLM = ggtern::grid.arrange(ternRodentHL, ternWGHL, gg_blank,
             ternRodentHS, ternWGHS, gg_blank,
             ternRodentBS, ternWGBS, gg_blank,
             ncol=3, nrow=3)

ggsave(ggternRLM, file=paste0("ggternRLM_ab",".png", sep = ""))
g", sep = ""))
