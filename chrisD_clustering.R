# chrisD_clustering.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Cluster the samples based on KCsmart profiles to
#               to see how similar they are
# -------------------------------------------------------------------

library(KCsmart)
library(gplots)
data(mmMirrorLocs)
setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')

# Calculate the Spm collection, use Tumor and Metastasis as difference.
# This is not that important as we will probably not look at that difference with
# a full comparative KC smart

classLabels <- rep(0, times=nrow(sampleInfo))
# Mets are 1
classLabels[grep('metastasis', sampleInfo$Type)] <- 1

TvsLNColl <- calcSpmCollection(data=allKC, mirrorLocs=mmMirrorLocs, cl=classLabels)
sampCorMat <- cor(TvsLNColl@data, use='na.or.complete')
labs <- with(sampleInfo, paste(DRSet, Type, Site))

png(file='Figures/corrMat.png', width=1000, height=1000)
heatmap(sampCorMat, scale='none', labRow=labs, labCol=labs, 
  col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
  ColSideColors=colors()[c(122, 148)][classLabels+1],
  RowSideColors=colors()[c(122, 148)][classLabels+1],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

