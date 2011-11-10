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
source('~/codeChris/smallProjects/chrisD/chrisD_functions.R')

setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')

# Check if sample info and data are in the same order.
sampleInfo$fileName <- paste(sampleInfo$Slide, sampleInfo$Spot, 
  sep='')
all.equal(colnames(allKC[, 3:ncol(allKC)]), sampleInfo$fileName)

# Functions

plotCorMap <- function (dataMat, classLabels, labs, plotTitle='Temp',
  plotCol=NULL) {
  
  require(gplots)
     
  corMat <- cor(dataMat, use='na.or.complete')

  if (is.null(plotCol)) {
    plotCol <- colorpanel(n=length(unique(classLabels)),
      low=colors()[24],
      high=colors()[349])
  }
  
  colIndex <- as.numeric(as.factor(classLabelsDR))

  heatmap(corMat, scale='none', labRow=labs, labCol=labs, 
    col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
    ColSideColors=plotCol[colIndex],
    RowSideColors=plotCol[colIndex],
    main=plotTitle)


}

# Do comparative KC-smart and save the result

KCclassLab <- rep(0, times=nrow(sampleInfo))
KCclassLab[sampleInfo$Type == 'metastasis'] <- 1
TvsLNColl <- calcSpmCollection(data=allKC, mirrorLocs=mmMirrorLocs, 
  cl=KCclassLab)
save(file='chrisD_compKCresult.Rda', list=c('TvsLNColl'))

# Calculate the Spm collection, use Tumor and Metastasis as difference.
# This is not that important as we will probably not look at that difference with
# a full comparative KC smart

classLabels <- rep(1, times=nrow(sampleInfo))
classLabels[grep('1', sampleInfo$DRSet)] <- 2
classLabels[grep('2', sampleInfo$DRSet)] <- 3
classLabels[grep('3', sampleInfo$DRSet)] <- 4
labs <- with(sampleInfo, paste(DRSet, Type, Site))

png(file='Figures/corrMatAll.png', width=1000, height=1000)
  plotCorMap(TvsLNColl@data, classLabels, labs,
    plotTitle='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

postscript(file="Figures/corrMatAll.eps", width=8, height=8, paper='special', horizontal=F)
  plotCorMap(TvsLNColl@data, classLabels, labs,
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




