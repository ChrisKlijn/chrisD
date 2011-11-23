# chrisD_metastasisAnalysis.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description: Analysis of aCGH data to find metastasis differences
# -------------------------------------------------------------------

# Run on medoid

setwd("/home/klijn/data/smallproj/chrisD/")
# Data
load('rawData_chrisD.Rda')
load('chrisD_segmentedDelta.Rda')
# Libraries and functions
library(multicore) # for parrallel preprocessing, works only on linux
source('~/gitCodeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/gitCodeChris/chrisD/chrisD_functions.R')

# -------------------------------------------
# Local functions

preprocessForDiffFigures <- function(threshold, segData, 
  diffKCrecMet, minMark=5) {
    
  # This function does all the plotting given a threshold on the delta
  # segments and a minimum number of markers

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

  returnList <- list(
    diffKCSeg=diffKCrecMet.Seg,
    threshold=threshold,
    minMark=minMark
    )

  return(returnList)

}

plotAllDiffFigures <- function(diffKCSegList) {
  
  # Make listitems in the the diffKCSegList local variables
  attach(diffKCSegList)
  
  recep <- c('R1', 'R2', 'R3')

  # --------
  # Recepients and their metas

  for (r in recep) {
    metMat <- as.matrix(diffKCSeg[,
      grepl(r, colnames(diffKCSeg))])
    png(file=
      paste('Figures/', r, 'priMetaDiff_', threshold, '.png', sep=''), 
      width=1000, height=800)
    plotMetMat(metMat, diffKCSeg, 
      plotTitle=paste(r, 'primary-metastasis CNA diff'))
    dev.off()
    postscript(file=
      paste('Figures/', r, 'priMetaDiff_', threshold, '.eps', sep=''), 
      width=10, height=8, paper='special', horizontal=F)
    plotMetMat(metMat, diffKCSeg, 
      plotTitle=paste(r, 'primary-metastasis CNA diff'))
    dev.off()
  }

  # ---------
  # All lung metas

  metMat <- as.matrix(diffKCSeg[,
     grepl('Lung', colnames(diffKCSeg))])
  
  png(file=paste('Figures/allLungMetsDiff', threshold, '.png', sep=''),
    width=1000, height=800)  
  plotMetMat(metMat, diffKCSeg, 
    plotTitle=
      paste('All lung primary-metastasis CNA diff -', threshold))
  dev.off()
  
  postscript(file=
    paste('Figures/allLungMetsDiff', threshold, '.eps', sep=''),
    width=10, height=8, paper='special', horizontal=F)
  plotMetMat(metMat, diffKCSeg, 
    plotTitle=
      paste('All lung primary-metastasis CNA diff -', threshold))
  dev.off()

  # ---------
  # All LN metas

  metMat <- as.matrix(diffKCSeg[,
     grepl('LN', colnames(diffKCSeg))])

  png(file=
    paste('Figures/allLNMetsDiff', threshold, '.png', sep=''), 
    width=1000, height=800)
  plotMetMat(metMat, diffKCSeg, 
      plotTitle=
        paste('All lymphnode primary-metastasis CNA diff -', threshold))
  dev.off()

  postscript(file=
    paste('Figures/allLNMetsDiff', threshold, '.eps', sep=''), 
    width=10, height=8, paper='special', horizontal=F)
  plotMetMat(metMat, diffKCSeg, 
      plotTitle=
        paste('All lymphnode primary-metastasis CNA diff -', threshold))
  dev.off()

}

# -------------------------------------------
# Actual code

# Get just the segments infomation
segData <- segRecMet$output

# Run preprocessing for two thresholds
thresholdList <- list(thres01=.1, thres02=.2)
KCSegList <- mclapply(thresholdList, preprocessForDiffFigures, segData,
  diffKCrecMet)

# Plot the delta figures using the preprocessed data
mclapply(KCSegList, plotAllDiffFigures)
