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
source('~/gitCodeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/gitCodeChris/chrisD/chrisD_functions.R')

# Only count segments with more than 5 probes and set everything between the thresholds to 0

threshold <- .2

segData <- segRecMet$output
segData$seg.mean[segData$num.mark < 5] <- 0
segData$seg.mean[abs(segData$seg.mean) < threshold] <- 0
segData$seg.mean[segData$seg.mean < -threshold] <- -1
segData$seg.mean[segData$seg.mean > threshold] <- 1

# Fix a weird renaming thing in DNAcopy
# replace the three consequtive dots with a dash
segData$ID <- gsub('[.]{3}', ' - ', segData$ID)

# Set probevalues to segment values

diffKCrecMet <- diffKCrecMet[
  order(diffKCrecMet$chrom, diffKCrecMet$maploc),]
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
  png(file=paste('Figures/', r, 'priMetaDiff_02_black.png', sep=''), 
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

png(file="Figures/allLNMetsDiff.png", width=1000, height=800)
metMat <- as.matrix(diffKCrecMet.Seg[,
   grepl('LN', colnames(diffKCrecMet.Seg))])
plotMetMat(metMat, diffKCrecMet.Seg, 
    plotTitle='All lymph node primary-metastasis CNA diff')
dev.off()

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