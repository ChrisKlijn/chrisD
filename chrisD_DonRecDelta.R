# chrisD_DonRecDelta.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Code for the difference between the donors and 
#               recepients.
#               Use quantile normalization to make the arrays
#               Comparable, then index the differences between the
#               paired samples.
# -------------------------------------------------------------------

# Run on medoid

library(KCsmart)
library(DNAcopy)
library(gplots)
library(preprocessCore)
data(mmMirrorLocs)
setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Quantile normalization
dataMatrix <- as.matrix(allKC[,3:ncol(allKC)])
dataMatrix <- normalize.quantiles(dataMatrix)
KCnorm <- allKC
KCnorm[,3:ncol(KCnorm)] <- dataMatrix

# Instantiate KC frame for the delta profile
diffKC <- KCnorm[,c(1,2)]

# Check if sample info is in the same order as the data
all.equal(paste(sampleInfo$Slide, sampleInfo$Spot, sep=''), colnames(KCnorm[,3:ncol(KCnorm)]))

# Make delta profiles for each of the pairs of donors with their primary tumors

# First define the donor hybs and name them
donorIndex <- grep('D[1-3]', sampleInfo$DRSet)
donorHybs <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[donorIndex]
names(donorHybs) <- sampleInfo$DRSet[donorIndex]

for (d in 1:length(donorHybs)) {
  
  # Match the donor 
  recepientIndex <- intersect(grep(gsub('D', 'R', names(donorNames))[d], sampleInfo$DRSet),
    which(sampleInfo$Site == 'Primary'))
  recepientHybs <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[recepientIndex]
  names(recepientHybs) <- sampleInfo$DRSet[recepientIndex]
  
  for (r in 1:length(recepientHybs)) {
    # Calculate the probewise difference between the donor and the recepient
    diffKC$temp <- KCnorm[, recepientHybs[r]] - KCnorm[, donorHybs[d]]
    # Name the column containing the difference
    comparisonTitle <- paste(names(recepientHybs)[r], '-', names(donorHybs)[d], sep='')
    colnames(diffKC) <- gsub('temp', comparisonTitle, colnames(diffKC))
  }
  
}  


#----------------------PLOTTING---------------------------

# Plot rainbowplots of the donor, transplant and the delta profile

plotCombs <- colnames(diffKC)[3:ncol(diffKC)]

altMirrorLocs <- mmMirrorLocs[-21]
attributes(altMirrorLocs) <- attributes(mmMirrorLocs)

for (i in plotCombs) {
  
  donIndex <- grep(gsub('.*[-]', '', i), sampleInfo$DRSet)
  recIndex <- intersect(grep(gsub('[-].*', '', i), sampleInfo$DRSet),
    which(sampleInfo$Site == 'Primary'))
  
  plotKC <- KCnorm[,c(1,2)]
  plotKC$donor <- KCnorm[,donIndex + 2]
  plotKC$recepient <- KCnorm[,recIndex + 2]
  plotKC$delta <- diffKC[,i]
  
  fileName <- paste('Figures/dotplot_Delta_', i, '.png', sep='')
  png(file=fileName, width=1200, height=1200)
    par(mfrow=c(3, 1))
    plotRawCghDotPlot(KCdataSet=plotKC, mirrorLocs=altMirrorLocs, doFilter=T, 
      samples=1, plotTitle=paste('Donor tumor', sampleInfo$DRSet[donIndex]))
    plotRawCghDotPlot(KCdataSet=plotKC, mirrorLocs=altMirrorLocs, doFilter=T, 
      samples=2, plotTitle=paste('Recepient tumor', sampleInfo$DRSet[recIndex]))
    plotRawCghDotPlot(KCdataSet=plotKC, mirrorLocs=altMirrorLocs, doFilter=T, 
      samples=3, plotTitle=paste('Delta', i))
  dev.off()
  
  # Save a CN file for IGV viewing
  
  fileNameCN <- paste('CNfiles/DonRec_', i, '.cn', sep='')
  
  colnames(plotKC) <- c(colnames(plotKC)[c(1,2)], paste(colnames(plotKC)[3:ncol(plotKC)], i))
  exportCNformat(KCdataSet=plotKC, fileName=fileNameCN)
  
}



