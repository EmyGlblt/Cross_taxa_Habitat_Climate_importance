library(Hmsc)

#array job at CSC with 5 parallel runs = one for each of five measurement scales of the habitat data:
args = commandArgs(trailingOnly=TRUE)
jobin = as.numeric(args[1])

scales = c("_500m", "_1km", "_2km", "_4km", "_6km")
scale = scales[jobin]

nChains = 2
samples = 250
thin = 1000

inputfile = paste0("bftestAB", scale, "_models_thin_", thin, "_samples_", samples, "_chains_", nChains, ".Rdata") 
load(file = inputfile)

MF = list()
MFCV = list()
WAIC = list()

for(n in 1:1){
  m = models[[n]]
  preds = computePredictedValues(m)
  MF[[n]] = evaluateModelFit(hM=m, predY=preds)
  partition = createPartition(m, nfolds = 2)
  preds = computePredictedValues(m, partition=partition, nParallel = nChains)
  MFCV[[n]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[n]] = computeWAIC(m)       
}

filename_out = paste("bftestperfAB", scale, "_MF_models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains2_",as.character(nChains),
                     ".Rdata",sep = "")

save(MF, MFCV, WAIC, modelnames, file = filename_out)
