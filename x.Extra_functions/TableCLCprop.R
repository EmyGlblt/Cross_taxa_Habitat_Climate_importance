library(tidyr)
library(dplyr)

setwd("E:/RECcII/Data_prep/RDATA")
load(file="Moth_subtest.RDATA")

head(Env.samp)
Mothmedian = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), median))
Mothmedian = round(as.matrix(Mothmedian), 3)
Mothmedian = as.data.frame(Mothmedian)
Mothmedian$taxa = 'Moth'


load(file="Bird_subtest.RDATA")

head(Env.samp)
Birdmedian = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), median))

Birdmedian = round(Birdmedian, 3)
Birdmedian = as.data.frame(Birdmedian)
Birdmedian$taxa = 'bird'

load(file="bf_subtest.RDATA")

head(Env.samp)
bfmedian = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), median))

bfmedian = round(bfmedian, 3)
bfmedian = as.data.frame(bfmedian)
bfmedian$taxa = 'Butterfly'


load(file="WG_subtest.RDATA")

head(Env.samp)
WGmedian = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), median))

WGmedian = round(WGmedian, 3)
WGmedian = as.data.frame(WGmedian)
WGmedian$taxa = 'Large mammal'


load(file="Rod2_subtest.RDATA")

head(Env.samp)
Rodmedian = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), median))

Rodmedian = round(Rodmedian, 3)
Rodmedian = as.data.frame(Rodmedian)
Rodmedian$taxa = 'Rodent'

All_median = rbind(Birdmedian, bfmedian, Mothmedian, WGmedian, Rodmedian)
head(All_median)

AM_rs = All_median %>%
  gather(var, Proportion, Crops_500m:Wetland_6km) %>%  ## Makes wide data long
  separate(var, c("habitat", "buffer"), sep = '_')    ## Splits up a column



AMbird = AM_rs[AM_rs$taxa =='bird', -1]
AMbird = AMbird[order(AMbird$habitat),]
colnames(AMbird)[3] = 'Bird'
AMbird2 = AMbird[!AMbird$habitat %in% c('Others', 'WaterBd'),]

AMbf = AM_rs[AM_rs$taxa =='Butterfly', -1]
AMbf = AMbf[order(AMbf$habitat),]
colnames(AMbf)[3] = 'Butterfly'
AMbf2 = AMbf[!AMbf$habitat %in% c('Others', 'WaterBd'),]

AMMoth = AM_rs[AM_rs$taxa =='Moth', -1]
AMMoth = AMMoth[order(AMMoth$habitat),]
colnames(AMMoth)[3] = 'Moth'
AMMoth2 = AMMoth[!AMMoth$habitat %in% c('Others', 'WaterBd'),]

AMWG = AM_rs[AM_rs$taxa =='Large mammal', -1]
AMWG = AMWG[order(AMWG$habitat),]
colnames(AMWG)[3] = 'Large mammal'
AMWG2 = AMWG[!AMWG$habitat %in% c('Others', 'WaterBd'),]

AMrod = AM_rs[AM_rs$taxa =='Rodent', -1]
AMrod = AMrod[order(AMrod$habitat),]
colnames(AMrod)[3] = 'Rodent'
AMrod2 = AMrod[!AMrod$habitat %in% c('Others', 'WaterBd'),]

Dataprop = data.frame(AMbird, 'Butterfly' = AMbf$Butterfly, 'Moth' = AMMoth$Moth, 
                      'Large mammal' = AMWG[,3], 'Rodent' = AMrod$Rodent)
head(Dataprop)
Dataprop2 = Dataprop[!Dataprop$habitat %in% c('Others', 'WaterBd'),]


library("rempsyc")

nice_table(Dataprop2)


library(formattable)
#rownames(Dataprop2) <- Dataprop2[,1]
formattable(Dataprop2)


library(kableExtra)
# kbl(Dataprop2)
# 
# maT = t(Dataprop2)
# 
# maT %>%
#   kbl() %>%
#   kable_classic_2(full_width = T)
# 
# 
# kbl(maT[-c(1:2),], caption = "Habitat proportion") %>%
#   kable_classic_2(full_width = F, html_font = "Cambria") %>%
#   add_header_above(c(" ", "500m" = 1, "1km" = 1, "2km" = 1, "4km" = 1, "6km" = 1, 
#                      "500m" = 1, "1km" = 1, "2km" = 1, "4km" = 1, "6km" = 1, 
#                      "500m" = 1, "1km" = 1, "2km" = 1, "4km" = 1, "6km" = 1, 
#                      "500m" = 1, "1km" = 1, "2km" = 1, "4km" = 1, "6km" = 1, 
#                      "500m" = 1, "1km" = 1, "2km" = 1, "4km" = 1, "6km" = 1)) %>%
#   add_header_above(c(" " = 1, "Crops" = 5, "Forest" = 5, 'HSP' = 5,
#                      'Urban' = 5, 'Wetland'= 5))



kbl(Dataprop2[, -c(1)], caption = "Habitat proportion median", row.names=FALSE) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  pack_rows("Crops", 1, 5) %>%
  pack_rows("Forest", 6, 10)%>%
  pack_rows("HSP", 11, 15) %>%
  pack_rows("Urban", 16, 20) %>%
  pack_rows("Wetland", 21, 25)



##  IQR

setwd("E:/RECcII/Data_prep/RDATA")
load(file="Moth_subtest.RDATA")

head(Env.samp)
MothIQR = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), IQR))
MothIQR = round(as.matrix(MothIQR), 3)
MothIQR = as.data.frame(MothIQR)
MothIQR$taxa = 'Moth'


load(file="Bird_subtest.RDATA")

head(Env.samp)
BirdIQR = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), IQR))

BirdIQR = round(BirdIQR, 3)
BirdIQR = as.data.frame(BirdIQR)
BirdIQR$taxa = 'bird'

load(file="bf_subtest.RDATA")

head(Env.samp)
bfIQR = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), IQR))

bfIQR = round(bfIQR, 3)
bfIQR = as.data.frame(bfIQR)
bfIQR$taxa = 'Butterfly'


load(file="WG_subtest.RDATA")

head(Env.samp)
WGIQR = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), IQR))

WGIQR = round(WGIQR, 3)
WGIQR = as.data.frame(WGIQR)
WGIQR$taxa = 'Large mammal'


load(file="Rod2_subtest.RDATA")

head(Env.samp)
RodIQR = Env.samp[, 32:66 ] %>%
  dplyr:: summarise(across(everything(), IQR))

RodIQR = round(RodIQR, 3)
RodIQR = as.data.frame(RodIQR)
RodIQR$taxa = 'Rodent'

All_IQR = rbind(BirdIQR, bfIQR, MothIQR, WGIQR, RodIQR)
head(All_IQR)

AM_rs = All_IQR %>%
  gather(var, Proportion, Crops_500m:Wetland_6km) %>%  ## Makes wide data long
  separate(var, c("habitat", "buffer"), sep = '_')    ## Splits up a column



AMbird = AM_rs[AM_rs$taxa =='bird', -1]
AMbird = AMbird[order(AMbird$habitat),]
colnames(AMbird)[3] = 'Bird'

AMbf = AM_rs[AM_rs$taxa =='Butterfly', -1]
AMbf = AMbf[order(AMbf$habitat),]
colnames(AMbf)[3] = 'Butterfly'

AMMoth = AM_rs[AM_rs$taxa =='Moth', -1]
AMMoth = AMMoth[order(AMMoth$habitat),]
colnames(AMMoth)[3] = 'Moth'

AMWG = AM_rs[AM_rs$taxa =='Large mammal', -1]
AMWG = AMWG[order(AMWG$habitat),]
colnames(AMWG)[3] = 'Large mammal'

AMrod = AM_rs[AM_rs$taxa =='Rodent', -1]
AMrod = AMrod[order(AMrod$habitat),]
colnames(AMrod)[3] = 'Rodent'


Dataprop = data.frame(AMbird, 'Butterfly' = AMbf$Butterfly, 'Moth' = AMMoth$Moth, 
                      'Large mammal' = AMWG[,3], 'Rodent' = AMrod$Rodent)
head(Dataprop)
Dataprop3 = Dataprop[!Dataprop$habitat %in% c('Others', 'WaterBd'),]


library("rempsyc")

nice_table(Dataprop3)


kbl(Dataprop3[, -c(1)], caption = "Habitat proportion IQR", row.names=FALSE) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  pack_rows("Crops", 1, 5) %>%
  pack_rows("Forest", 6, 10)%>%
  pack_rows("HSP", 11, 15) %>%
  pack_rows("Urban", 16, 20) %>%
  pack_rows("Wetland", 21, 25)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Both together?
knitr::kables(list(
  kbl(Dataprop2[, -c(1)], caption = "Habitat proportion median", row.names=FALSE) %>%
    kable_classic_2(full_width = F, html_font = "Cambria") %>%
    pack_rows("Crops", 1, 5) %>%
    pack_rows("Forest", 6, 10)%>%
    pack_rows("HSP", 11, 15) %>%
    pack_rows("Urban", 16, 20) %>%
    pack_rows("Wetland", 21, 25),
  kbl(Dataprop3[, -c(1)], caption = "Habitat proportion IQR", row.names=FALSE) %>%
    kable_classic_2(full_width = F, html_font = "Cambria") %>%
    pack_rows("Crops", 1, 5) %>%
    pack_rows("Forest", 6, 10)%>%
    pack_rows("HSP", 11, 15) %>%
    pack_rows("Urban", 16, 20) %>%
    pack_rows("Wetland", 21, 25)
  
)
) %>% kable_styling()
