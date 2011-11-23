# chrisD_compKC.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: aCGH data Chris D - matched primary and metastasis
# Description: Run comparative KC smart on the data and save the
#              result. This can then be used for the indiviually
#              smoothed profiles
# -------------------------------------------------------------------

# Run on medoid

library(KCsmart)
data(mmMirrorLocs)
setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')

# Check if sample info and data are in the same order.
sampleInfo$fileName <- paste(sampleInfo$Slide, sampleInfo$Spot, 
  sep='')
all.equal(colnames(allKC[, 3:ncol(allKC)]), sampleInfo$fileName)

# Do comparative KC-smart and save the result
# Calculate the Spm collection, use Tumor and Metastasis as difference.

KCclassLab <- rep(0, times=nrow(sampleInfo))
KCclassLab[sampleInfo$Type == 'metastasis'] <- 1
TvsLNColl <- calcSpmCollection(data=allKC, mirrorLocs=mmMirrorLocs, 
  cl=KCclassLab)
save(file='chrisD_compKCresult.Rda', list=c('TvsLNColl'))
