# chrisD_segmentKC.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description: segment the KC data
# -------------------------------------------------------------------

# Run on medoid

# Working dir
setwd("/home/klijn/data/smallproj/chrisD/")

# Libraries
library(DNAcopy)

# Data
load('rawData_chrisD.Rda')

# Check if sampleinfo and the data are ordered the same
all.equal(colnames(allKC[,3:ncol(allKC)]), paste(sampleInfo$Slide,
  sampleInfo$Spot, sep=''))

CNA.allKC <- CNA(as.matrix(allKC[,3:ncol(allKC)]), allKC$chrom, allKC$maploc, data.type=c("logratio"), 
  sampleid=colnames(allKC)[3:ncol(allKC)])
CNA.allKC.smoothed <- smooth.CNA(CNA.allKC)
allKCseg <- segment(CNA.allKC.smoothed, verbose=1, undo.splits='sdundo')

# Fix the segmented data, remove the appended X
colnames(allKCseg$data) <- gsub('X', '', colnames(allKCseg$data))
allKCseg$output$ID <- gsub('X','',allKCseg$output$ID)

save(file='chrisD_segmentedKC.Rda', list=c('allKCseg'))

