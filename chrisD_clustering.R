# chrisD_clustering.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Cluster the samples based on KCsmart profiles to
#               to see how similar they are
# -------------------------------------------------------------------

# Run on medoid

library(KCsmart)
library(gplots)
data(mmMirrorLocs)
source('~/gitCodeChris/chrisD/chrisD_functions.R')

setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')
load('chrisD_compKCresult.Rda')

# Check if sample info and data are in the same order.
sampleInfo$fileName <- paste(sampleInfo$Slide, sampleInfo$Spot, 
  sep='')
all.equal(colnames(allKC[, 3:ncol(allKC)]), sampleInfo$fileName)

# Generate readable labes for the cluster plots
classLabels <- rep(1, times=nrow(sampleInfo))
classLabels[grep('1', sampleInfo$DRSet)] <- 2
classLabels[grep('2', sampleInfo$DRSet)] <- 3
classLabels[grep('3', sampleInfo$DRSet)] <- 4
labs <- with(sampleInfo, paste(DRSet, Type, Site))

# All samples

png(file='Figures/corrMatAll.png', width=1000, height=1000)
  plotCorMap(TvsLNColl@data, classLabels, labs,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

postscript(file="Figures/corrMatAll.eps", width=8, height=8, paper='special', horizontal=F)
  plotCorMap(TvsLNColl@data, classLabels, labs,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

# All Recepients and their metastases
# Excluded the donors and the n.a. sample

subsetAllRM <- grepl('R', sampleInfo$DRSet)
TvsLNCollAllRM <- TvsLNColl@data[, subsetAllRM]
labsAllRM <- labs[subsetAllRM]
classLabelsAllRM <- classLabels[subsetAllRM]

png(file='Figures/corrMatAllRM.png', width=1000, height=1000)
  plotCorMap(TvsLNCollAllRM, classLabelsAllRM, labsAllRM,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

postscript(file="Figures/corrMatAllRM.eps", width=8, height=8, paper='special', horizontal=F)
  plotCorMap(TvsLNCollAllRM, classLabelsAllRM, labsAllRM,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

# Donors and receptients

subsetDR <- grepl('R|D', sampleInfo$DRSet) & 
  grepl('tumor', sampleInfo$Type)
TvsLNCollDR <- TvsLNColl@data[, subsetDR]
labsDR <- labs[subsetDR]
classLabelsDR <- classLabels[subsetDR]

png(file='Figures/corrMatDR.png', width=1000, height=1000)
  plotCorMap(TvsLNCollDR, classLabelsDR, labsDR,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

postscript(file="Figures/corrMatDR.eps", width=8, height=8, paper='special', horizontal=F)
  plotCorMap(TvsLNCollDR, classLabelsDR, labsDR,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

# Plot heatmaps for each recepient and its metastases

uniqDonors <- c('1','2','3')

for (i in uniqDonors) {
  
  subsetRec <- grepl(i, sampleInfo$DRSet)
  recData <- TvsLNColl@data[, subsetRec]
  labsRec <- labs[subsetRec]
  classLabelsRec <- sampleInfo$Type[subsetRec]
  tempTitle <- paste('Tumors donor', i)
  fileName <- paste('Figures/corrMatDonor', i, sep='')

  postscript(file=paste(fileName, '.eps', sep=''), 
    width=5, height=5, paper='special', horizontal=F)
  plotCorMap(recData, classLabelsRec, labsRec, plotTitle=tempTitle)
  dev.off()
  
  png(file=paste(fileName, '.png', sep=''), width=500, height=500)
  plotCorMap(recData, classLabelsRec, labsRec, plotTitle=tempTitle)
  dev.off()
}




