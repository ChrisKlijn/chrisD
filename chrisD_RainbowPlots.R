# chrisD_RainbowPlots.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Plot rainbow plots to check the data
# -------------------------------------------------------------------

# Plot rainbowplots for checking

library(KCsmart)
data(mmMirrorLocs)
setwd("/home/klijn/data/smallproj/chrisD/")
source('~/gitCodeChris/generalFunctionsR/chris_cghdata_analysis.R')
load('rawData_chrisD.Rda')

# Remove the CMF control 
sampleInfo <- sampleInfo[-which(sampleInfo$DRSet == 'n.a.'),]

altMirrorLocs <- mmMirrorLocs[-21]
attributes(altMirrorLocs) <- attributes(mmMirrorLocs)

# All samples separately

for (i in 1:nrow(sampleInfo)) {
  sampleName <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[i]
  plotTitle <- with(sampleInfo, paste('Sample:', DRSet[i], Site[i], Type[i]))
  fileName <- paste('Figures/dotplot_', sampleName, '_screen.png', sep='')
  
  png(file=fileName, width=1200, height=600) 
    plotRawCghDotPlot(KCdataSet=allKC[,c('chrom', 'maploc', sampleName)], 
      mirrorLocs=altMirrorLocs, doFilter=T, samples=1, plotTitle=plotTitle)
  dev.off()
}

# Donor tumor and their primary transplants

uniqTransplNum <- unique(gsub('[a-z|A-Z|n.a.]', '', sampleInfo$DRSet))

for (i in uniqTransplNum) {
  # Find the indices for the primary transplants
  transplInd <- intersect(grep(paste('R', i, sep=''), sampleInfo$DRSet), 
    grep('Primary', sampleInfo$Site))
  fileName <- paste('Figures/dotplot_DonorTranspl_Series', i, '.png', sep='')
  png(file=fileName, width=1200, height=400*(length(transplInd) + 1))
    par(mfrow=c(length(transplInd) + 1, 1))
    donorName <-  paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[
      grep(paste('D', i, sep=''), sampleInfo$DRSet)]
    plotRawCghDotPlot(KCdataSet=allKC[,c('chrom', 'maploc', donorName)], 
      mirrorLocs=altMirrorLocs, doFilter=T, samples=1, plotTitle=paste('Donor tumor D', i))
    
    for (j in transplInd) {
      transplName <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[j]
      plotRawCghDotPlot(KCdataSet=allKC[,c('chrom', 'maploc', transplName)], 
        mirrorLocs=altMirrorLocs, doFilter=T, samples=1, 
        plotTitle=paste('Primary Transplant' ,sampleInfo$DRSet[j]))
    }
  dev.off()
}

# Primary transplants and their metastases

uniqPrimTranspl <- unique(sampleInfo$DRSet[grep('R', sampleInfo$DRSet)])

for (i in uniqPrimTranspl) {
  # Find the indices for the primary transplants and the metastases
  primInd <- intersect(grep(i, sampleInfo$DRSet), grep('Primary', sampleInfo$Site))
  metInd <- intersect(grep(i, sampleInfo$DRSet), grep('metastasis', sampleInfo$Type))
  fileName <- paste('Figures/dotplot_TransplMet_Series', i, '.png', sep='')
  png(file=fileName, width=1200, height=400*(length(metInd)+1))
    par(mfrow=c(length(metInd) + 1, 1))
    primName <-  paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[primInd]
    plotRawCghDotPlot(KCdataSet=allKC[,c('chrom', 'maploc', primName)], 
      mirrorLocs=altMirrorLocs, doFilter=T, samples=1, plotTitle=paste('Primary transplant', i))
    
    for (j in metInd) {
      metName <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[j]
      plotRawCghDotPlot(KCdataSet=allKC[,c('chrom', 'maploc', metName)], 
        mirrorLocs=altMirrorLocs, doFilter=T, samples=1, 
        plotTitle=paste('Metastasis' ,sampleInfo$DRSet[j], sampleInfo$Site[j]))
    }
  dev.off()
}

