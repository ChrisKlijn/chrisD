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

# Register multicores
registerDoMC(8)

# Check if sampleinfo and the data are ordered the same
all.equal(colnames(allKC[,3:ncol(allKC)]), paste(sampleInfo$Slide,
  sampleInfo$Spot, sep=''))

# Deltas between donor tumors and recepient tumors
# Aggregate per tumor

# First define the donor hybs and name them

donorIndex <- grep('D[1-3]', sampleInfo$DRSet)
names(donorIndex) <- sampleInfo$DRSet[donorIndex]
donorList <- vector(mode='list', length=length(donorIndex))
names(donorList) <- names(donorIndex)

for (d in donorIndex) {
  recepientIndex <- with(sampleInfo, intersect(grep(gsub('D', 'R', DRSet[d]), DRSet), which(Site == 'Primary')))
  names(recepientIndex) <- sampleInfo$DRSet[recepientIndex]

  recepientList <- vector(mode='list', length=length(recepientIndex))
  names(recepientList) <- names(recepientIndex)

  foreach(r=recepientIndex) %dopar% {
    
    recepientList[names(r)] <- deltaLinear(comb=c(d+2, r+2), KC=allKC, KCseg=allKCseg)

  }

  donorList[names(d)] <- recepientList

}
