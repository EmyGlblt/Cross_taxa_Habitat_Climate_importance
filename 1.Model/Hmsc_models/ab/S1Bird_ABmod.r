library(Hmsc)

load(file="prepBird_Jan24_Strat2_1km.RDATA")
# Loads pre-prepared Y, X, S and xy coordinate matrices.


sp_spread = sp40Y_spread
dim(sp_spread)
head(sp_spread)

Y = sp_spread[, -c(1,2)]
# Convert Y abundance into presence-absence data
Ypa = 1*(Y>0)
Ypa[is.na(Ypa)] = 0
#range(colMeans(Ypa))
#range(colSums(Ypa))
# Some species are very rare, others ubiquitous

#Here, I exclude species with < 40 observations and > total obs -40.
threshold = 40
rarespp = which(colSums(Ypa)<threshold)
commonspp = which(colSums(Ypa)>(nrow(Ypa)-threshold))
exclude = c(rarespp, commonspp)

Ypa = Ypa[, -exclude]

dim(Ypa)


## for abundance data - pure model

Yab = Y[, colnames(Y) %in% colnames(Ypa)]
Yab[is.na(Yab)] = 0

#Yab = log(Yab)


XFormula = ~Effort + poly(precip, 2, raw = T) + poly(temp, 2, raw = T) + poly(Snow_count, 2, raw = T) + FracCV + DivL + 
  poly(Crops, 2, raw = T) + poly(Forest, 2, raw = T) + poly(HSP, 2, raw = T) + poly(Urban, 2, raw = T) +  
  poly(Wetland, 2, raw = T)


### scales test


chs = c("_ch1", "_ch2", "_ch3", "_ch4")
#array job at CSC with 5 parallel runs = one for each of five measurement chs of the habitat data:
args = commandArgs(trailingOnly=TRUE)
jobin = as.numeric(args[1])

ch = chs[jobin]

# Env info
head(Env.br2)
XData = Env.br2[, c(3:4, 6:8, 12:14, 16, 18, 23)]
summary(XData)

XData$FracCV = as.numeric(XData$FracCV)
XData$DivL = as.numeric(XData$DivL)


Ypa = Ypa[!is.na(XData$Effort),]
Env.br2 = Env.br2[!is.na(XData$Effort),]
XData = XData[!is.na(XData$Effort),]

dim(XData)
dim(Ypa)

studyDesign = data.frame(site = as.factor(Env.br2$SiteID), 
                         year = as.factor(Env.br2$YearStart), 
                         bg = as.factor(Env.br2$BZ2))
dim(studyDesign)




# coord

birdCoords <- read.csv("BirdTransects_centroids.csv", header=T)
head(birdCoords)

#birdCoords$Numero = as.numeric(birdCoords$Numero)
birdCoords2 = birdCoords[birdCoords$Numero %in% unique(Env.br2$SiteID), ]

coord = birdCoords2[order(match(birdCoords2$Numero, unique(Env.br2$SiteID))),]

xy = data.frame(unique(coord[, 3:4]))
rownames(xy) = unique(coord[, 2])



xyKnots = constructKnots(sData = xy, nKnots = 20, knotDist = NULL, minKnotDist = NULL)
rL.site = HmscRandomLevel(sData = xy, sMethod = "GPP", sKnot = xyKnots)

Yr = studyDesign$year
rL.year = HmscRandomLevel(units = levels(Yr))

Biogeo = studyDesign$bg
rL.biogeo = HmscRandomLevel(units = levels(Biogeo))
m1 = Hmsc(Y=Yab, XData = XData,  XFormula = XFormula,     ## 
          distr="lognormal poisson", YScale=TRUE,  ## Yscale = TRUE only neede really for normal distrib does not influence here
          studyDesign=studyDesign,
          ranLevels={list("site"=rL.site, "year" = rL.year, "bg" = rL.biogeo)})

models = list(m1)
modelnames = c("abundance_condPA")

save(models,modelnames,file = paste0("unfitbirdAB_models", ch, ".Rdata"))

#Fit models
nChains = 1
samples = 250
thin = 1000

ptm = proc.time()
for(n in 1:1)
{
  m = models[[n]]
  m = sampleMcmc(m, samples = samples, thin=thin,
                 adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
                 transient = ceiling(0.5*samples*thin),
                 nChains = nChains, nParallel = nChains)
  models[[n]] = m
}

computational.time2 = proc.time() - ptm
computational.time2



filename_out = paste("BirdAB", ch, "_models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(models, modelnames, file=filename_out)

