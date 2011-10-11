# chrisD_DonRecDeltaLinear.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Code for the difference between the donors and 
#               recepients.
#               Use linear regression normalization to make the arrays
#               Comparable, then index the differences between the
#               paired samples.
# -------------------------------------------------------------------

# Run on medoid

# Working dir
setwd("/home/klijn/data/smallproj/chrisD/")

# Code
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_delta_functions.R')
library(KCsmart)
library(DNAcopy)
library(foreach)
library(doMC)

# Data
load('rawData_chrisD.Rda')
load('chrisD_segmentedKC.Rda')
data(mmMirrorLocs)
altMirrorLocs <- mmMirrorLocs[-21]
attributes(altMirrorLocs) <- attributes(mmMirrorLocs)

# Fix the segmented data, remove the appended X
colnames(allKCseg$data) <- gsub('X', '', colnames(allKCseg$data))

# Register multicores
registerDoMC(8)

# Check if sampleinfo and the data are ordered the same
all.equal(colnames(allKC[,3:ncol(allKC)]), paste(sampleInfo$Slide,
  sampleInfo$Spot, sep=''))

# Remove the negative control sample

negativeControl <- grep('NegativeControl', sampleInfo$SampleID)
if (length(negativeControl) > 0) {
  allKC <- allKC[,-(negativeControl+2)]
  sampleInfo <- sampleInfo[-negativeControl,] 
}

# Assign tumor numbers to sample, NA warning for the negative control
sampleInfo$tumNum <- as.numeric(gsub('[A-Z|a-z]','', sampleInfo$DRSet))
# Assign CGHID
sampleInfo$CGHID <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')

# Deltas between donor tumors and recepient tumors
# Aggregate per tumor

# First define the donor hybs and name them

tumNums <- unique(sampleInfo$tumNum)
diffList <- vector(mode='list', length=length(tumNums))
names(diffList) <- paste('T', tumNums, sep='')

for (t in tumNums) {
  tempSampInfo <- subset(sampleInfo, tumNum == t)
  donorSample <- tempSampInfo$CGHID[grep('D', tempSampInfo$DRSet)]
  recepientSamples <- tempSampInfo$CGHID[grepl('R', tempSampInfo$DRSet) &
   tempSampInfo$Site == 'Primary']
  names(recepientSamples) <- tempSampInfo$DRSet[grepl('R', 
    tempSampInfo$DRSet) & tempSampInfo$Site == 'Primary']
  
  resultList <- vector(mode='list', length=length(recepientSamples))
  names(resultList) <- names(recepientSamples)

  for (r in 1:length(recepientSamples)) {
    tempKC <- allKC[,c('chrom', 'maploc', recepientSamples[r],
      donorSample)]
    tempSeg <- subset(allKCseg, 
      samplelist=c(recepientSamples[r], donorSample))
    resultList[[names(recepientSamples)[r]]] <-
      deltaLinear(comb= c(recepientSamples[r], donorSample), 
      tempKC, tempSeg, thres=.2)
  }

  diffList[[paste('T', t, sep='')]] <- 
    cbind(allKC[, c('chrom', 'maploc')], resultList)

}

x11()
plotRawCghDotPlot(KCdataSet=allKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=24, chromosomes=6, setcex=2)
plotRawCghDotPlot(KCdataSet=allKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=25, chromosomes=6, setcex=2)
plotRawCghDotPlot(KCdataSet=diffList[['T1']], mirrorLocs=altMirrorLocs, doFilter=T, samples=1
