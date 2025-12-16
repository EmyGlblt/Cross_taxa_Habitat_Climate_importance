#-.......................................................................................

## Useful libraries .................................................
library(stringr)  # data manip
library(dplyr)
library(tidyr)
library(class)
library(purrr)

library(ggplot2)  # for plotting
library(ggridges)

library(maps)   # map and spatial objects
library(mapdata)
library(sp)
library(mapsFinland)
library(maptools)
library(sf)  # opening spatial object


## Finland map 
par(mfrow=c(1,1), oma= c(2,2,2,2))
Finmap =  maps::map('worldHires','Finland', plot=FALSE, fill = TRUE)
class(Finmap)

Finmap$names[2] = "Finland1:?" 
Finmap$names[5] = "Finland2:?" 
Finmap$names[6] = "Finland3:?" 


Fin.sp = map2SpatialPolygons(map = Finmap, IDs = Finmap$names, 
                             proj4string = CRS("+init=epsg:4326"), 
                             checkHoles=FALSE)
plot(Fin.sp)
class(Fin.sp)

newMFin.sp <- spTransform(Fin.sp, CRS("+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs"))
class(newMFin.sp)

## boreal regions
setwd("P:/h570/rec_team/Cartographic_material/BioclimaticZonesFinland")
Fin_bioclim = st_read("BioclimaticZonesFinland.shp")
st_crs(Fin_bioclim) # should be the same crs as the grid

plot(Fin_bioclim)

mapDataALL.sf = st_as_sf(newMFin.sp, coords = c('X', 'Y'))
# dissolve_sf <- st_combine(mapDataALL.sf)
# cast_sf = st_cast(st_cast(dissolve_sf,"MULTIPOLYGON"),"POLYGON")
mapFin.sf = st_transform(mapDataALL.sf, st_crs(Fin_bioclim))

Fin_join <- st_intersection(mapFin.sf[1], Fin_bioclim)  # order is kept
# Fin_join2 = Fin_join
# st_geometry(Fin_join2) <- NULL ## makes a dataframe
# FIn_grid_BR = Fin_join2

Fin_join$BR = NA

Fin_join[Fin_join$Nimi == "Hemiboreaalinen vyöhyke",]$BR = 'Southern boreal' ## combine the two southern most regions
Fin_join[Fin_join$Nimi == "Eteläboreaalinen vyöhyke",]$BR = 'Southern boreal'
Fin_join[Fin_join$Nimi == "Keskiboreaalinen vyöhyke",]$BR = 'Middle boreal'
Fin_join[Fin_join$Nimi == "Pohjoisboreaalinen vyöhyke",]$BR = 'Northern boreal'
Fin_join$BR = factor(Fin_join$BR, levels = c('Northern boreal', 'Middle boreal', 'Southern boreal'))

ggplot(data = Fin_join, aes(fill = BR)) +
  geom_sf() + scale_fill_manual(values = rev(c("#FDE725FF","#35B779FF", "#31688EFF")))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_void()




ggplot(data = Fin_join, ) +
  geom_sf() + #scale_fill_manual(values = rep('white',3))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_minimal()


lower_bound <- 0.05
upper_bound <- 1 - lower_bound

ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1),
                color = "darkorange3", linewidth = 3) + ylab("") +
  scale_y_continuous(breaks = NULL) +
  stat_function(
    fun = function(x) ifelse(qnorm(lower_bound) < x  & x < qnorm(upper_bound), dnorm(x), NA), 
    geom = "area", 
    fill = "darkorange4", 
    alpha = 0.5
  ) +
  scale_x_continuous(
    breaks = qnorm(c(lower_bound, 0.5, upper_bound)),
    labels = round(qnorm(c(lower_bound, 0.5, upper_bound)), digits = 2)
  ) +
  labs(x = "z", y = "Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 


#................................ with taxa

setwd("D:/Helsinki/RECcII/CSC_HMSC/fullmodel/Results")

load('bftest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelbf = models
coord_bf <- modelbf[[1]]$ranLevels$site$s

coord_bfsf = st_as_sf(x = coord_bf, 
                      coords = c('Longitude', 'Latitude') ,
                      crs = st_crs(Fin_join))

ggplot() +
  geom_sf(data = Fin_join) + 
  geom_sf(data = Fin_join[,8], aes(fill = BR)) + 
  scale_fill_manual(values = rev(c("#FDE725FF","#35B779FF", "#31688EFF")))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_void() +
  geom_point(data = coord_bf, aes(x=Longitude, y=Latitude))


load('MothTest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelmoth = models
coord_moth <- modelmoth[[1]]$ranLevels$site$s

coord_mothsf = st_as_sf(x = coord_moth, 
                      coords = c('Longitude', 'Latitude') ,
                      crs = st_crs(Fin_join))

ggplot() +
  geom_sf(data = Fin_join) + 
  geom_sf(data = Fin_join[,8], aes(fill = BR)) + 
  scale_fill_manual(values = rev(c("#FDE725FF","#35B779FF", "#31688EFF")))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_void() +
  geom_point(data = coord_moth, aes(x=Longitude, y=Latitude))



load('Birdtest_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelBirds = models
coord_Birds <- modelBirds[[1]]$ranLevels$site$s

coord_Birdssf = st_as_sf(x = coord_Birds, 
                        coords = c('Longitude', 'Latitude') ,
                        crs = st_crs(Fin_join))

ggplot() +
  geom_sf(data = Fin_join) + 
  geom_sf(data = Fin_join[,8], aes(fill = BR)) + 
  scale_fill_manual(values = rev(c("#FDE725FF","#35B779FF", "#31688EFF")))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_void() +
  geom_point(data = coord_Birds, aes(x=Longitude, y=Latitude))



load('WintGtest_ch1_models_thin_1000_samples_250_chainsN1_1.Rdata')
modelWGs = models
coord_WGs <- modelWGs[[1]]$ranLevels$site$s

coord_WGssf = st_as_sf(x = coord_WGs, 
                         coords = c('Longitude', 'Latitude') ,
                         crs = st_crs(Fin_join))

ggplot() +
  geom_sf(data = Fin_join) + 
  geom_sf(data = Fin_join[,8], aes(fill = BR)) + 
  scale_fill_manual(values = rev(c("#FDE725FF","#35B779FF", "#31688EFF")))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_void() +
  geom_point(data = coord_WGs, aes(x=Longitude, y=Latitude))



load('Rodtest2_ch1_models_thin_1000_samples_250_chains_1.Rdata')
modelRods = models
coord_Rods <- modelRods[[1]]$ranLevels$site$s

coord_Rodssf = st_as_sf(x = coord_Rods, 
                         coords = c('Longitude', 'Latitude') ,
                         crs = st_crs(Fin_join))

ggplot() +
  geom_sf(data = Fin_join) + 
  geom_sf(data = Fin_join[,8], aes(fill = BR)) + 
  scale_fill_manual(values = rev(c("#FDE725FF","#35B779FF", "#31688EFF")))+ 
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  theme(plot.title = element_text(size=7.5), axis.text=element_text(size=6),
        axis.title = element_text(size = 6), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) + theme_void() +
  geom_point(data = coord_Rods, aes(x=Longitude, y=Latitude))

