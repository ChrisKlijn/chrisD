
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
setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')

# Check if sample info and data are in the same order.
sampleInfo$fileName <- paste(sampleInfo$Slide, sampleInfo$Spot, 
  sep='')
all.equal(colnames(allKC[, 3:ncol(allKC)]), sampleInfo$fileName)

# Calculate the Spm collection, use Tumor and Metastasis as difference.
# This is not that important as we will probably not look at that difference with
# a full comparative KC smart

classLabels <- rep(1, times=nrow(sampleInfo))
# Mets are 1
classLabels[grep('1', sampleInfo$DRSet)] <- 2
classLabels[grep('2', sampleInfo$DRSet)] <- 3
classLabels[grep('3', sampleInfo$DRSet)] <- 4

plotCol <- colors()[c(24, 181, 200, 242)]

TvsLNColl <- calcSpmCollection(data=allKC, mirrorLocs=mmMirrorLocs, cl=classLabels)
sampCorMat <- cor(TvsLNColl@data, use='na.or.complete')
labs <- with(sampleInfo, paste(DRSet, Type, Site))

png(file='Figures/corrMatAll.png', width=1000, height=1000)
heatmap(sampCorMat, scale='none', labRow=labs, labCol=labs, 
  col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
  ColSideColors=plotCol[classLabels],
  RowSideColors=plotCol[classLabels],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

postscript(file="Figures/corrMatAll.eps", width=8, height=8, paper='special', horizontal=F)
heatmap(sampCorMat, scale='none', labRow=labs, labCol=labs, 
  col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
  ColSideColors=plotCol[classLabels],
  RowSideColors=plotCol[classLabels],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

# Donors and receptients

subsetDR <- grepl('R|D', sampleInfo$DRSet) & 
  grepl('tumor', sampleInfo$Type)
TvsLNCollDR <- TvsLNColl@data[, subsetDR]
sampCorMatDR <- cor(TvsLNCollDR, use='na.or.complete')
labsDR <- labs[subsetDR]
classLabelsDR <- classLabels[subsetDR]

png(file='Figures/corrMatDR.png', width=1000, height=1000)
heatmap(sampCorMatDR, scale='none', labRow=labsDR, labCol=labsDR, 
  col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
  ColSideColors=plotCol[classLabelsDR],
  RowSideColors=plotCol[classLabelsDR],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

postscript(file="Figures/corrMatDR.eps", width=8, height=8, paper='special', horizontal=F)
heatmap(sampCorMatDR, scale='none', labRow=labsDR, labCol=labsDR, 
  col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
  ColSideColors=plotCol[classLabelsDR],
  RowSideColors=plotCol[classLabelsDR],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()






