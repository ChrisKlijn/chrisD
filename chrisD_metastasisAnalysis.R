# chrisD_metastasisAnalysis.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description: Analysis of aCGH data to find metastasis differences
# -------------------------------------------------------------------

# Run on medoid

library(GenomicRanges)
library(gplots)
data(mmMirrorLocs)
setwd("/home/klijn/data/smallproj/chrisD/")
load('rawData_chrisD.Rda')
load('chrisD_segmentedDelta.Rda')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/smallProjects/chrisD/chrisD_functions.R')

# Only count segments with more than 5 probes and set everything above or below .1 to 0
segData <- segRecMet$output
segData$seg.mean[segData$num.mark < 5] <- 0
segData$seg.mean[abs(segData$seg.mean) < .1] <- 0
segData$seg.mean[segData$seg.mean < -.1] <- -1
segData$seg.mean[segData$seg.mean > .1] <- 1

# Set probevalues to segment values

diffKCrecMet.Seg <- diffKCrecMet

for (s in 1:nrow(segData)) {
  colNum <- which(colnames(diffKCrecMet.Seg) == segData$ID[s])
  probeInd <- which(diffKCrecMet.Seg$chrom == segData$chrom[s] &
    diffKCrecMet.Seg$maploc >= segData$loc.start[s] &
    diffKCrecMet.Seg$maploc <= segData$loc.end[s])
  diffKCrecMet.Seg[probeInd, colNum] <- segData$seg.mean[s]
}

recep <- c('R1', 'R2', 'R3')

for (r in recep) {
  png(file=paste('Figures/', r, 'priMetaDiff.png', sep=''), 
    width=1000, height=800)
  metMat <- as.matrix(diffKCrecMet.Seg[,
    grepl(r, colnames(diffKCrecMet.Seg))])
  plotMetMat(metMat, diffKCrecMet.Seg, 
    plotTitle=paste(r, 'primary-metastasis CNA diff'))
  dev.off()
}

png(file="Figures/allLungMetsDiff.png", width=1000, height=800)
metMat <- as.matrix(diffKCrecMet.Seg[,
   grepl('Lung', colnames(diffKCrecMet.Seg))])
plotMetMat(metMat, diffKCrecMet.Seg, 
    plotTitle='All lung primary-metastasis CNA diff')
dev.off()


chromEnds <- cumsum(unlist(lapply(mmMirrorLocs, max)))
probeChromEnds <- cumsum(tapply(diffKCrecMet.Seg$maploc, 
  diffKCrecMet.Seg$chrom, length))

# Order the matrix on tumor ID, and far or close metas
tumorSiteInfo <- unlist(strsplit(colnames(metMat), split = ' - '))
orderFrame <- data.frame(
  site=tumorSiteInfo[seq(1,length(tumorSiteInfo),2)],
  tumorID=tumorSiteInfo[seq(2,length(tumorSiteInfo),2)],
  stringsAsFactors=F)
orderFrame$close <- !grepl('Axillary', orderFrame$site)
orderVect <- order(orderFrame$tumorID, orderFrame$close, orderFrame$site)
orderFrame <- orderFrame[orderVect,]
metMat <- metMat[,orderVect]

# Set colors
gainLossCols <- colorpanel(3, low='green', mid='black', high='orange')
closeCols <- colorpanel(2, low='blue', high='red')
tumorCols <- colorpanel(length(unique(orderFrame$tumorID)), 
  low=colors()[1], high=colors()[24])
names(tumorCols) <- unique(orderFrame$tumorID)

# intialize plot
x11()
plot(0,0, xlim=c(0, nrow(metMat)+20000), ylim=c(0, ncol(metMat)+2),
  type='n',axes=F, xlab=NA, ylab=NA)

# Use run length encoding to draw probe-index based rectangle
# use additional rectangles for tumor site and tumor ID
for (j in 1:ncol(metMat)) {
  tumVect <- rle(metMat[,j])
  start <- 0
  for (i in 1:length(tumVect$length)) {
    rect(start, j-1, tumVect$length[i]+start, j, 
      col=cols[tumVect$value[i]+2], border=NA)
    start <- tumVect$length[i] + start + 1
  }
  rect(start+1000, j-1, start+4000, j, 
    col=closeCols[orderFrame$close[j]+1])
  rect(start+5000, j-1, start+8000, j, 
    col=tumorCols[orderFrame$tumorID[j]])
}
segments(x0=probeChromEnds, y0=0, y1=ncol(metMat), col='gray')
text(x=probeChromEnds-1000, y=ncol(metMat), labels=names(probeChromEnds),offset=0.1, pos=3, cex=.5)
text(x=max(probeChromEnds)+8000, y=seq(.5, ncol(metMat)-.5, 1),
  labels=colnames(metMat), pos=4, cex=.6)


# -----------


segDataFilter <- segData[segData$num.mark > 5 & 
  abs(segData$seg.mean) > .1,]

# Get a ranges object

filterSegs <- with(segDataFilter, GRanges(
  seqnames = Rle(paste('chr', chrom)),
  ranges=IRanges(start=loc.start,end=loc.end), 
  seg.mean=seg.mean, 
  num.mark=num.mark,
  ID=ID))

tumSegs <- filterSegs[grepl('R1', elementMetadata(filterSegs)$ID),]