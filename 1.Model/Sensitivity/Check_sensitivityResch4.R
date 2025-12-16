##
#
#             Check Results HMSC buffer sensitivity analysis
#
#
################################################################################

library(ggplot2)
library(RColorBrewer)
#library(purrr)

setwd("E:/RECcII/CSC_HMSC/Results/thin1000")

## Butterfly

# 500m
load('Birdtestperf_500m_MF_models_thin_1000_samples_250_chainsN1_1.Rdata')
MF_500_1 = MF
MFCV_500_1 = MFCV
WAIC_500_1 = WAIC

load('Birdtestperf_500m_MF_models_thin_1000_samples_250_chainsN2_1.Rdata')
MF_500_2 = MF
MFCV_500_2 = MFCV
WAIC_500_2 = WAIC

load('Birdtestperf_500m_MF_models_thin_1000_samples_250_chainsN3_1.Rdata')
MF_500_3 = MF
MFCV_500_3 = MFCV
WAIC_500_3 = WAIC

load('Birdtestperf_500m_MF_models_thin_1000_samples_250_chainsN4_1.Rdata')
MF_500_4 = MF
MFCV_500_4 = MFCV
WAIC_500_4 = WAIC

MFtjurR2 = list(MF_500_1[[1]]$TjurR2, MF_500_2[[1]]$TjurR2, MF_500_3[[1]]$TjurR2, MF_500_4[[1]]$TjurR2)
MFCVtjurR2 = list(MFCV_500_1[[1]]$TjurR2, MFCV_500_2[[1]]$TjurR2, MFCV_500_3[[1]]$TjurR2, MFCV_500_4[[1]]$TjurR2)
WAIC = list(WAIC_500_1[[1]], WAIC_500_2[[1]], WAIC_500_3[[1]], WAIC_500_4[[1]])

MF_500 = Reduce('+', MFtjurR2)/4
MFCV_500 = Reduce('+', MFCVtjurR2)/4
WAIC_500 = Reduce('+', WAIC)/4



# 1km
load('Birdtestperf_1km_MF_models_thin_1000_samples_250_chainsN1_1.Rdata')
MF_1k_1 = MF
MFCV_1k_1 = MFCV
WAIC_1k_1 = WAIC

load('Birdtestperf_1km_MF_models_thin_1000_samples_250_chainsN2_1.Rdata')
MF_1k_2 = MF
MFCV_1k_2 = MFCV
WAIC_1k_2 = WAIC

load('Birdtestperf_1km_MF_models_thin_1000_samples_250_chainsN3_1.Rdata')
MF_1k_3 = MF
MFCV_1k_3 = MFCV
WAIC_1k_3 = WAIC

load('Birdtestperf_1km_MF_models_thin_1000_samples_250_chainsN4_1.Rdata')
MF_1k_4 = MF
MFCV_1k_4 = MFCV
WAIC_1k_4 = WAIC

MFtjurR2 = list(MF_1k_1[[1]]$TjurR2, MF_1k_2[[1]]$TjurR2, MF_1k_3[[1]]$TjurR2, MF_1k_4[[1]]$TjurR2)
MFCVtjurR2 = list(MFCV_1k_1[[1]]$TjurR2, MFCV_1k_2[[1]]$TjurR2, MFCV_1k_3[[1]]$TjurR2, MFCV_1k_4[[1]]$TjurR2)
WAIC = list(WAIC_1k_1[[1]], WAIC_1k_2[[1]], WAIC_1k_3[[1]], WAIC_1k_4[[1]])

MF_1k = Reduce('+', MFtjurR2)/4
MFCV_1k = Reduce('+', MFCVtjurR2)/4
WAIC_1k = Reduce('+', WAIC)/4


# 2km
load('Birdtestperf_2km_MF_models_thin_1000_samples_250_chainsN1_1.Rdata')
MF_2k_1 = MF
MFCV_2k_1 = MFCV
WAIC_2k_1 = WAIC

load('Birdtestperf_2km_MF_models_thin_1000_samples_250_chainsN2_1.Rdata')
MF_2k_2 = MF
MFCV_2k_2 = MFCV
WAIC_2k_2 = WAIC

load('Birdtestperf_2km_MF_models_thin_1000_samples_250_chainsN3_1.Rdata')
MF_2k_3 = MF
MFCV_2k_3 = MFCV
WAIC_2k_3 = WAIC

load('Birdtestperf_2km_MF_models_thin_1000_samples_250_chainsN4_1.Rdata')
MF_2k_4 = MF
MFCV_2k_4 = MFCV
WAIC_2k_4 = WAIC

MFtjurR2 = list(MF_2k_1[[1]]$TjurR2, MF_2k_2[[1]]$TjurR2, MF_2k_3[[1]]$TjurR2, MF_2k_4[[1]]$TjurR2)
MFCVtjurR2 = list(MFCV_2k_1[[1]]$TjurR2, MFCV_2k_2[[1]]$TjurR2, MFCV_2k_3[[1]]$TjurR2, MFCV_2k_4[[1]]$TjurR2)
WAIC = list(WAIC_2k_1[[1]], WAIC_2k_2[[1]], WAIC_2k_3[[1]], WAIC_2k_4[[1]])

MF_2k = Reduce('+', MFtjurR2)/4
MFCV_2k = Reduce('+', MFCVtjurR2)/4
WAIC_2k = Reduce('+', WAIC)/4


# 4km
load('Birdtestperf_4km_MF_models_thin_1000_samples_250_chainsN1_1.Rdata')
MF_4k_1 = MF
MFCV_4k_1 = MFCV
WAIC_4k_1 = WAIC

load('Birdtestperf_4km_MF_models_thin_1000_samples_250_chainsN2_1.Rdata')
MF_4k_2 = MF
MFCV_4k_2 = MFCV
WAIC_4k_2 = WAIC

load('Birdtestperf_4km_MF_models_thin_1000_samples_250_chainsN3_1.Rdata')
MF_4k_3 = MF
MFCV_4k_3 = MFCV
WAIC_4k_3 = WAIC

load('Birdtestperf_4km_MF_models_thin_1000_samples_250_chainsN4_1.Rdata')
MF_4k_4 = MF
MFCV_4k_4 = MFCV
WAIC_4k_4 = WAIC

MFtjurR2 = list(MF_4k_1[[1]]$TjurR2, MF_4k_2[[1]]$TjurR2, MF_4k_3[[1]]$TjurR2, MF_4k_4[[1]]$TjurR2)
MFCVtjurR2 = list(MFCV_4k_1[[1]]$TjurR2, MFCV_4k_2[[1]]$TjurR2, MFCV_4k_3[[1]]$TjurR2, MFCV_4k_4[[1]]$TjurR2)
WAIC = list(WAIC_4k_1[[1]], WAIC_4k_2[[1]], WAIC_4k_3[[1]], WAIC_4k_4[[1]])

MF_4k = Reduce('+', MFtjurR2)/4
MFCV_4k = Reduce('+', MFCVtjurR2)/4
WAIC_4k = Reduce('+', WAIC)/4


# 6km
load('Birdtestperf_6km_MF_models_thin_1000_samples_250_chainsN1_1.Rdata')
MF_6k_1 = MF
MFCV_6k_1 = MFCV
WAIC_6k_1 = WAIC

load('Birdtestperf_6km_MF_models_thin_1000_samples_250_chainsN2_1.Rdata')
MF_6k_2 = MF
MFCV_6k_2 = MFCV
WAIC_6k_2 = WAIC

load('Birdtestperf_6km_MF_models_thin_1000_samples_250_chainsN3_1.Rdata')
MF_6k_3 = MF
MFCV_6k_3 = MFCV
WAIC_6k_3 = WAIC

load('Birdtestperf_6km_MF_models_thin_1000_samples_250_chainsN4_1.Rdata')
MF_6k_4 = MF
MFCV_6k_4 = MFCV
WAIC_6k_4 = WAIC

MFtjurR2 = list(MF_6k_1[[1]]$TjurR2, MF_6k_2[[1]]$TjurR2, MF_6k_3[[1]]$TjurR2, MF_6k_4[[1]]$TjurR2)
MFCVtjurR2 = list(MFCV_6k_1[[1]]$TjurR2, MFCV_6k_2[[1]]$TjurR2, MFCV_6k_3[[1]]$TjurR2, MFCV_6k_4[[1]]$TjurR2)
WAIC = list(WAIC_6k_1[[1]], WAIC_6k_2[[1]], WAIC_6k_3[[1]], WAIC_6k_4[[1]])

MF_6k = Reduce('+', MFtjurR2)/4
MFCV_6k = Reduce('+', MFCVtjurR2)/4
WAIC_6k = Reduce('+', WAIC)/4


### comparison

# WAIC
WAIC.df = data.frame(WAIC = c(WAIC_500, WAIC_1k, WAIC_2k,
                              WAIC_4k, WAIC_6k),
                     Scale = c('500m', '1km', '2km',
                               '4km', '6km'))

WAIC.df$Scale = factor(WAIC.df$Scale, levels = c('500m', '1km', '2km',
                                                 '4km', '6km'))


ggplot(WAIC.df, aes(x = Scale, y = WAIC)) + 
  geom_point(aes(size=2)) + theme_bw()


TRCV.df = data.frame(MFCV = c(MFCV_500, MFCV_1k,
                              MFCV_2k, MFCV_4k, MFCV_6k),
                     Scale = c(rep('500m', length(MFCV_500)), 
                               rep('1km', length(MFCV_1k)), 
                               rep('2km', length(MFCV_2k)), 
                               rep('4km', length(MFCV_4k)), 
                               rep('6km', length(MFCV_6k))),
                     WAIC = c(rep(WAIC_500, length(MFCV_500)), 
                              rep(WAIC_1k, length(MFCV_1k)), 
                              rep(WAIC_2k, length(MFCV_2k)), 
                              rep(WAIC_4k, length(MFCV_4k)), 
                              rep(WAIC_6k, length(MFCV_6k))))

TRCV.df$Scale = factor(TRCV.df$Scale, levels = c('500m', '1km','2km', 
                                                 '4km', '6km'))

colorchoice = brewer.pal(8, 'Dark2')
colorchoice = colorchoice[c(3, 5:8)]

ylim.prim <- c(summary(TRCV.df$MFCV)[1], summary(TRCV.df$MFCV)[6])   # min max
ylim.sec <- c(summary(TRCV.df$WAIC)[1], summary(TRCV.df$WAIC)[6])    # min max

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1] # there was a bug here



ggplot(TRCV.df, aes(x = Scale, fill=Scale)) + 
  geom_violin(aes(y = MFCV), draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(aes(y = a+WAIC*b), shape=18, size=3, col='red2') + 
  scale_fill_manual(values=colorchoice) +
  scale_y_continuous("TjurR2 CV", sec.axis = sec_axis(~(.-a)/b, name = "WAIC")) +
  theme_bw() + theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14,face="bold"))


ggplot(TRCV.df, aes(x = Scale, fill=Scale)) + 
  geom_violin(aes(y = MFCV), draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_point(aes(y = a+WAIC*b, shape=Scale, size=Scale)) + 
  scale_fill_manual(values=colorchoice) +
  scale_y_continuous("TjurR2 CV") +
  theme_bw() + theme(axis.text=element_text(size=12),
                     axis.title=element_text(size=14,face="bold"))
