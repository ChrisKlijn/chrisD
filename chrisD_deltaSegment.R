# chrisD_deltaSegment.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description: Segmentation of Delta profiles code.
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

# ----- FUNCTIONS ------

diffSeg <- function(KC) {

  CNA.KC <- CNA(as.matrix(KC[,3:ncol(KC)]), KC$chrom, KC$maploc, data.type=c("logratio"), 
    sampleid=colnames(KC)[3:ncol(KC)])
  CNA.KC.smoothed <- smooth.CNA(CNA.KC)
  KCseg <- segment(CNA.KC.smoothed, verbose=1, undo.splits='sdundo')
  
  return(KCseg)
  
}

diffPlot <- function(diffKC, sampleNo=1, squash=F, segKC=NULL, posThres=0, negThres=0, squashVal=2) {
  
  # Plots a collection of chromosomes with relative segments of gain/loss 
  
  # To Do:
  # - Check if chromosomes are numerical
  
  # Plotting
  sampleNo = sampleNo + 2
  
  if (is.null(segKC)) {
    posInd <- diffKC[,sampleNo] > posThres
    negInd <- diffKC[,sampleNo] < negThres
    zeroInd <- diffKC[,sampleNo] <= posThres & diffKC[,sampleNo] >= negThres
  }
  else {
    posInd <- segKC[,sampleNo] > posThres
    negInd <- segKC[,sampleNo] < negThres
    zeroInd <- segKC[,sampleNo] == 0
  }
  
  # Squashing 
  # This means scaling the data between -.5 and .5, by dividing max
  
  if (squash == T) {
    cat('Squashing...\n')
    diffKC[,sampleNo] <- diffKC[,sampleNo]/(squashVal*max(diffKC[,sampleNo]))
  }
  
  maxChrom <- max(diffKC$chrom)
  maxMaploc <- max(diffKC$maploc)
  
  # Set up the plot
  plot(type='n', x=0, y=0, xlim=c(0, maxMaploc + .1*maxMaploc), 
    ylim=c(0, maxChrom), 
    main=names(diffKC)[sampleNo], ylab=NA, xlab = 'Genomic Coordinate', 
    xaxt = 'n', yaxt='n', frame.plot=FALSE)
  abline(v=0, col='black')
  abline(h=seq(0.5, 20.5, by=1), col='lightgray')
  abline(h=seq(1, 20, by=1), col='lightgray', lty='dotted')
  text(x = maxMaploc + .1*maxMaploc, 
    y = seq(.70, maxChrom - .30, by=1),
    labels=paste('chr', seq(1,maxChrom, by=1)), pos=2, col='gray', cex=.5, offset=0)
    

  points(diffKC$maploc[zeroInd], 
    diffKC[zeroInd,sampleNo]+diffKC$chrom[zeroInd], 
    pch='.')
  points(diffKC$maploc[posInd], 
    diffKC[posInd,sampleNo]+diffKC$chrom[posInd], 
    cex=2, col='red', pch='.')
  points(diffKC$maploc[negInd], 
    diffKC[negInd,sampleNo]+diffKC$chrom[negInd], 
    cex=2, col='green', pch='.')


}



 # ----- END FUNCTIONS -----

# Quantile normalization
dataMatrix <- as.matrix(allKC[,3:ncol(allKC)])
dataMatrix <- normalize.quantiles(dataMatrix)
KCnorm <- allKC
KCnorm[,3:ncol(KCnorm)] <- dataMatrix

# Instantiate KC frame for the delta profiles
diffKCdonRec <- KCnorm[,c(1,2)]
diffKCrecMet <- KCnorm[,c(1,2)]

# Check if sample info is in the same order as the data
all.equal(paste(sampleInfo$Slide, sampleInfo$Spot, sep=''), colnames(KCnorm[,3:ncol(KCnorm)]))

# Make delta profiles for each of the pairs of donors with their primary tumors

# First define the donor hybs and name them
donorIndex <- grep('D[1-3]', sampleInfo$DRSet)
donorHybs <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[donorIndex]
names(donorHybs) <- sampleInfo$DRSet[donorIndex]

for (d in 1:length(donorHybs)) {
  
  # Match the donor 
  recepientIndex <- intersect(grep(gsub('D', 'R', names(donorHybs))[d], sampleInfo$DRSet),
    which(sampleInfo$Site == 'Primary'))
  recepientHybs <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[recepientIndex]
  names(recepientHybs) <- sampleInfo$DRSet[recepientIndex]
  
  for (r in 1:length(recepientHybs)) {
    # Calculate the probewise difference between the donor and the recepient
    diffKCdonRec$temp <- KCnorm[, recepientHybs[r]] - KCnorm[, donorHybs[d]]
    # Name the column containing the difference
    comparisonTitle <- paste(names(recepientHybs)[r], '-', names(donorHybs)[d], sep='')
    colnames(diffKCdonRec) <- gsub('temp', comparisonTitle, colnames(diffKCdonRec))
  }
  
}  

# Make delta profiles for each of the primary tumors and their metastases

# First define unique 
primaryIndex <- intersect(grep('R', sampleInfo$DRSet), which(sampleInfo$Site == 'Primary'))
primaryHybs <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[primaryIndex]
names(primaryHybs) <- sampleInfo$DRSet[primaryIndex]

for (d in 1:length(primaryHybs)) {
  
  # Find the metas
  hybName <- names(primaryHybs[d])
  metIndex <- which(sampleInfo$DRSet == hybName & sampleInfo$Type == 'metastasis')
  metHybs <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')[metIndex]
  names(metHybs) <- paste(sampleInfo$Site, sampleInfo$DRSet, sep=' - ')[metIndex]
  
  for (m in 1:length(metHybs)) {
    # Calculate the probewise difference between the primary and the metastasis
    diffKCrecMet$temp <- KCnorm[, metHybs[m]] - KCnorm[, primaryHybs[d]]
    # Name the column containing the difference
    colnames(diffKCrecMet) <- gsub('temp', names(metHybs)[m], colnames(diffKCrecMet))
  }
}  

# Segmentation of delta profiles

# input: KCdiff data frame

segRecMet <- diffSeg(diffKCrecMet)

save(file='chrisD_segmentedDelta.Rda', list=c('segRecMet'))

# Make a diffSeg object with 0 -> no sig seg, 1 -> sig seg
# Do it simple for now -- slow for loop
# IRanges might have the solution

segData <- segRecMet$output

# Only count segments with more than 5 probes and set everything above or below .1 to 0
segData$seg.mean[segData$num.mark < 5] <- 0
segData$seg.mean[abs(segData$seg.mean) < .1] <- 0

# Fix a weird renaming thing in DNAcopy
segData$ID <- gsub('[.]{3}', ' - ', segData$ID)

diffKCrecMet.Seg <- diffKCrecMet

for (s in 1:nrow(segData)) {
  colNum <- which(colnames(diffKCrecMet.Seg) == segData$ID[s])
  probeInd <- which(diffKCrecMet.Seg$chrom == segData$chrom[s] &
    diffKCrecMet.Seg$maploc >= segData$loc.start[s] &
    diffKCrecMet.Seg$maploc <= segData$loc.end[s])
  diffKCrecMet.Seg[probeInd, colNum] <- segData$seg.mean[s]
}

for (c in 3:ncol(diffKCrecMet.Seg)) {
  
  diffKCrecMet.Seg[diffKCrecMet.Seg[,c] > 0,c] <- .4
  diffKCrecMet.Seg[diffKCrecMet.Seg[,c] < 0,c] <- -.4

}

# Sum over probes to find multiple changed probes

a <- rowSums(diffKCrecMet.Seg[,3:ncol(diffKCrecMet.Seg)])

# Plot overview for all primary tumors

uniqTransplant <- 
  sampleInfo$DRSet[intersect(grep('R', sampleInfo$DRSet), which(sampleInfo$Site == 'Primary'))]

# Data Deltaplot

for (tr in uniqTransplant) {
  sampleNumbers <- grep(tr, colnames(diffKCrecMet)) - 2
  fileNamePNG <- paste('Figures/deltaPlot_', tr, '.png', sep='')
  png(file=fileNamePNG, width=400*length(sampleNumbers), height=600)
  par(mfrow=c(1,length(sampleNumbers)))
  for (s in sampleNumbers) {
    diffPlot(diffKCrecMet, sampleNo=s, segKC=diffKCrecMet.Seg, squash=T)
  }
  dev.off()
}

# Segment Deltaplot

for (tr in uniqTransplant) {
  sampleNumbers <- grep(tr, colnames(diffKCrecMet)) - 2
  fileNamePNG <- paste('Figures/deltaPlot_seg_', tr, '.png', sep='')
  png(file=fileNamePNG, width=400*length(sampleNumbers), height=600)
  par(mfrow=c(1,length(sampleNumbers)))
  for (s in sampleNumbers) {
    diffPlot(diffKCrecMet.Seg, sampleNo=s)
  }
  dev.off()
}

diffPlot(diffKCrecMet.Seg, sampleNo=2)
diffPlot(diffKCrecMet, sampleNo=2, segKC=diffKCrecMet.Seg, squash=T, squashVal=4)

# Scatterplots


altMirrorLocs <- mmMirrorLocs[-21]
attributes(altMirrorLocs) <- attributes(mmMirrorLocs)

tempPlotFrame <- cbind(allKC[,c(donorHybs[1], recepientHybs[1])], 
  KCnorm[,c(donorHybs[1], recepientHybs[1])])
tempPlotFrame.filter <- apply(tempPlotFrame, 2, function(x) {filter(x, rep(1, 10)/10)})

plotFrame <- cbind(.5 * (tempPlotFrame[,1] + tempPlotFrame[,2]),
  tempPlotFrame[,1] - tempPlotFrame[,2], 
  .5 * (tempPlotFrame[,3] + tempPlotFrame[,4]),
  tempPlotFrame[,3] - tempPlotFrame[,4])

plotFrame.filter <- cbind(.5 * (tempPlotFrame.filter[,1] + tempPlotFrame.filter[,2]),
  tempPlotFrame.filter[,1] - tempPlotFrame.filter[,2], 
  .5 * (tempPlotFrame.filter[,3] + tempPlotFrame.filter[,4]),
  tempPlotFrame.filter[,3] - tempPlotFrame.filter[,4])

layout(matrix(c(1,1,2,3,4,5,6,6), 4, 2, byrow = TRUE))
plotRawCghDotPlot(KCdataSet=diffKCdonRec, mirrorLocs=altMirrorLocs, doFilter=T, 
      samples=1)
plot(plotFrame[,1], plotFrame[,2], pch='.', main='non-qnorm')
chr8probes <- allKC$chrom == 8
points(plotFrame[chr8probes,1], plotFrame[chr8probes, 2], col=rainbow(n=20)[8], 
  cex=2, pch='.')
abline(h=0)
plot(plotFrame[,3], plotFrame[,4], pch='.', main='qnorm')
chr8probes <- allKC$chrom == 8
points(plotFrame[chr8probes,3], plotFrame[chr8probes, 4], col=rainbow(n=20)[8], 
  cex=2, pch='.')
abline(h=0)
plot(plotFrame.filter[,1], plotFrame.filter[,2], pch='.', main='non-qnorm filter')
chr8probes <- allKC$chrom == 8
points(plotFrame.filter[chr8probes,1], plotFrame.filter[chr8probes, 2], col=rainbow(n=20)[8], 
  cex=2, pch='.')
abline(h=0)
plot(plotFrame.filter[,3], plotFrame.filter[,4], pch='.', main='qnorm filter')
chr8probes <- allKC$chrom == 8
points(plotFrame.filter[chr8probes,3], plotFrame.filter[chr8probes, 4], col=rainbow(n=20)[8], 
  cex=2, pch='.')
abline(h=0)
plot(plotFrame.filter[,3], plotFrame.filter[,4], pch=19, cex=1, main='qnorm filter',
  col=rainbow(n=20, alpha=.1)[allKC$chrom]) 
abline(h=0)



plotCols <- col2rgb(rainbow(n = 20), alpha=T)
plotCols[4,] <- .1

plot(plotFrame.filter[,3], plotFrame.filter[,4], pch=21, cex=.4, main='qnorm filter', 
  bg=plotCols[,allKC$chrom], lwd=0)

plot(plotFrame.filter[,3], plotFrame.filter[,4], pch='.', cex=1, main='qnorm filter',col=rgb(0,0,0,.1)) 
plot(plotFrame.filter[,3], plotFrame.filter[,4], pch=19, cex=1, main='qnorm filter',
  col=rainbow(n=20, alpha=.1)[allKC$chrom]) 
abline(h=0)

smoothScatter(plotFrame.filter[,3], plotFrame.filter[,4], col=rainbow(n=20)[allKC$chrom])
