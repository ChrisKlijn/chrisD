# chrisD_DonRecDeltaLinear.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Code for the difference between the donors and 
#               recepients.
#               Use linear regression normalization to make the arrays
#               Comparable, then index the differences between the
#               paired samples.
# -------------------------------------------------------------------

# Run on medoid

# Working dir
setwd("/home/klijn/data/smallproj/chrisD/")

# Code
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_delta_functions.R')
library(KCsmart)
library(DNAcopy)
library(robustbase)
library(preprocessCore)

# Data
load('rawData_chrisD.Rda')
load('chrisD_segmentedKC.Rda')
data(mmMirrorLocs)
altMirrorLocs <- mmMirrorLocs[-21]
attributes(altMirrorLocs) <- attributes(mmMirrorLocs)

# Fix the segmented data, remove the appended X
colnames(allKCseg$data) <- gsub('X', '', colnames(allKCseg$data))
# Sort allKC on chromosome and maploc
allKC <- allKC[order(allKC$chrom, allKC$maploc),]

# Check if sampleinfo and the data are ordered the same
all.equal(colnames(allKC[,3:ncol(allKC)]), paste(sampleInfo$Slide,
  sampleInfo$Spot, sep=''))

# Remove the negative control sample

negativeControl <- grep('NegativeControl', sampleInfo$SampleID)
if (length(negativeControl) > 0) {
  allKC <- allKC[,-(negativeControl+2)]
  sampleInfo <- sampleInfo[-negativeControl,] 
}

# Assign tumor numbers to sample, NA warning for the negative control
sampleInfo$tumNum <- as.numeric(gsub('[A-Z|a-z]','', sampleInfo$DRSet))
# Assign CGHID
sampleInfo$CGHID <- paste(sampleInfo$Slide, sampleInfo$Spot, sep='')

# Deltas between donor tumors and recepient tumors
# Aggregate per tumor

# First define the donor hybs and name them

tumNums <- unique(sampleInfo$tumNum)
diffList <- vector(mode='list', length=length(tumNums))
names(diffList) <- paste('T', tumNums, sep='')

for (t in tumNums) {
  tempSampInfo <- subset(sampleInfo, tumNum == t)
  donorSample <- tempSampInfo$CGHID[grep('D', tempSampInfo$DRSet)]
  recepientSamples <- tempSampInfo$CGHID[grepl('R', tempSampInfo$DRSet) &
   tempSampInfo$Site == 'Primary']
  names(recepientSamples) <- tempSampInfo$DRSet[grepl('R', 
    tempSampInfo$DRSet) & tempSampInfo$Site == 'Primary']
  
  resultList <- vector(mode='list', length=length(recepientSamples))
  names(resultList) <- names(recepientSamples)

  for (r in 1:length(recepientSamples)) {
    tempKC <- allKC[,c('chrom', 'maploc', recepientSamples[r],
      donorSample)]
    tempSeg <- subset(allKCseg, 
      samplelist=c(recepientSamples[r], donorSample))
    resultList[[names(recepientSamples)[r]]] <-
      deltaLinear(comb= c(recepientSamples[r], donorSample), 
      tempKC, tempSeg, thres=.2)
  }

  diffList[[paste('T', t, sep='')]] <- 
    cbind(allKC[, c('chrom', 'maploc')], resultList)

}

par(mfrow=c(3,1))
plotRawCghDotPlot(KCdataSet=allKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=24, chromosomes=6, setcex=2)
plotRawCghDotPlot(KCdataSet=allKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=25, chromosomes=6, setcex=2)
plotRawCghDotPlot(KCdataSet=diffList[['T1']], mirrorLocs=altMirrorLocs, doFilter=T, samples=1, chromosomes=6, setcex=2)

par(mfrow=c(3,1))
plotRawCghDotPlot(KCdataSet=allKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=24, setcex=2)
plotRawCghDotPlot(KCdataSet=allKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=25, setcex=2)
plotRawCghDotPlot(KCdataSet=diffList[['T1']], mirrorLocs=altMirrorLocs, doFilter=T, samples=1, setcex=2)

par(mfrow=c(3,1))
for (i in 1:(length(diffList[['T3']])-2)) {
  plotRawCghDotPlot(KCdataSet=diffList[['T3']], mirrorLocs=altMirrorLocs, doFilter=T, samples=i, setcex=2)

}

par(mfrow=c(3,1))
for (i in 1:(length(diffList[['T3']])-2)) {
  plotRawCghDotPlot(KCdataSet=diffList[['T3']], mirrorLocs=altMirrorLocs, doFilter=T, samples=i, setcex=2, chromosomes=13)
}


# Tests

# Linear models

comb <- c('437371A02', '437371A01')
smallKC <- allKC[,c('chrom', 'maploc', comb)]
smallSeg <- subset(allKCseg, samplelist=comb)

smallFreq <- glFrequency(xout=smallSeg, threshold=1)
ind <- smallFreq$gain == 1 | smallFreq$loss == -1
fitrob <- lmrob(smallKC[ind,4] ~ smallKC[ind,3])
fitlm <- lm(smallKC[ind,4] ~ smallKC[ind,3])

plot(smallKC[ind,3], smallKC[ind, 4], pch='.', cex=2, 
  col=smallKC$chrom[ind], main='fitted lms on selected probes')
abline(a=coef(fitrob)[1], b=coef(fitrob)[2], col='red')
abline(a=coef(fitlm)[1], b=coef(fitlm)[2], col='blue')
abline(a=0, b=1, col='black', lty='dotted')
legend('topleft', legend=c('lm', 'robust lm', 'y=x'), 
  col=c('red', 'blue', 'black'), lty=c('solid', 'solid', 'dotted'), 
  horiz=T)


# Quantile Normalization

dataMatrix <- as.matrix(smallKC[,3:ncol(smallKC)])
dataMatrix <- normalize.quantiles(dataMatrix)
qnormKC <- smallKC
qnormKC[,3:ncol(qnormKC)] <- dataMatrix

par(mfrow=c(1,2))
plot(smallKC[,3],smallKC[,4], pch='.', cex=2, col=smallKC$chrom,
  main='Pre-qnorm')
abline(a=0, b=1, col='black', lty='dotted')
plot(qnormKC[,3],qnormKC[,4], pch='.', cex=2, col=qnormKC$chrom,
  main='qnorm')
abline(a=0, b=1, col='black', lty='dotted')

# set probes to segment values

segKC <- smallKC

for (i in 1:nrow(smallSeg$output)) {
  probesInSeg <- with(smallSeg$output, segKC$chrom == chrom[i] & 
    segKC$maploc >= loc.start[i] &
    segKC$maploc < loc.end[i])

  segKC[probesInSeg, smallSeg$output$ID[i]] <- 
    smallSeg$output$seg.mean[i]

}

# Fit on hard cutoff instead of the MAD based cutoff from
# glFrequency

ind2 <- abs(segKC[,3]) > .2 & abs(segKC[,4]) > .2
fitrob2 <- lmrob(smallKC[ind2,4] ~ smallKC[ind2,3])
fitlm2 <- lm(smallKC[ind2,4] ~ smallKC[ind2,3])
plot(smallKC[ind2,3], smallKC[ind2, 4], pch='.', cex=2, 
  col=smallKC$chrom[ind2], main='Selection on seg.mean, not MAD')
abline(a=coef(fitrob2)[1], b=coef(fitrob2)[2], col='red')
abline(a=coef(fitlm2)[1], b=coef(fitlm2)[2], col='blue')
abline(a=0, b=1, col='black', lty='dotted')
legend('topleft', legend=c('lm', 'robust lm', 'y=x'), 
  col=c('red', 'blue', 'black'), lty=c('solid', 'solid', 'dotted'), 
  horiz=T)

# Fit on probes set to their segmean. (so, many probes have equal values).
# This is quite probably not a good choice 

fitrob3 <- lmrob(segKC[ind2,4] ~ segKC[ind2,3])
fitlm3 <- lm(segKC[ind2,4] ~ segKC[ind2,3])
plot(segKC[ind2,3], segKC[ind2, 4], pch='.', cex=2, 
  col=smallKC$chrom[ind2], main='Fitted on probes set to seg.mean')
abline(a=coef(fitrob3)[1], b=coef(fitrob3)[2], col='red')
abline(a=coef(fitlm3)[1], b=coef(fitlm3)[2], col='blue')
abline(a=0, b=1, col='black', lty='dotted')
legend('topleft', legend=c('lm', 'robust lm', 'y=x'), 
  col=c('red', 'blue', 'black'), lty=c('solid', 'solid', 'dotted'), 
  horiz=T)

segKC$diffSegNorm <- (segKC[,3] + coef(fitrob3)[[1]]) * coef(fitrob3)[[2]] - segKC[,4]

# Visualize the segmean correlation:

par(mfrow=c(3,1))
plotRawCghDotPlot(KCdataSet=segKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=1, setcex=2)
plotRawCghDotPlot(KCdataSet=segKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=2, setcex=2)
plotRawCghDotPlot(KCdataSet=segKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=3, setcex=2, setylim=c(-1,1))

par(mfrow=c(3,1))
plotRawCghDotPlot(KCdataSet=segKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=1, setcex=6, chromosomes=6)
plotRawCghDotPlot(KCdataSet=segKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=2, setcex=6, chromosomes=6)
plotRawCghDotPlot(KCdataSet=segKC, mirrorLocs=altMirrorLocs, doFilter=T, samples=3, setcex=6, setylim=c(-1,1), chromosomes=6)

allKCsegProbe <- setProbeToSeg(allKC, allKCseg)

for (t in tumNums) {
  
  tempSampInfo <- subset(sampleInfo, tumNum == t)
  donorSample <- tempSampInfo$CGHID[grep('D', tempSampInfo$DRSet)]
  recepientSamples <- tempSampInfo$CGHID[grepl('R', tempSampInfo$DRSet) &
   tempSampInfo$Site == 'Primary']
  names(recepientSamples) <- tempSampInfo$DRSet[grepl('R', 
    tempSampInfo$DRSet) & tempSampInfo$Site == 'Primary']
  
  resultList <- vector(mode='list', length=3)
  names(resultList) <- c('donor', 'recepients', 'delta')

  resultList$donor <- allKCsegProbe[,c('chrom', 'maploc', donorSample)]
  resultList$recepients <- 
    allKCsegProbe[,c('chrom', 'maploc', recepientSamples)]
  resultKC <- allKCsegProbe[, c('chrom', 'maploc')]
  
  for (r in 1:length(recepientSamples)) {
    tempKC <- allKC[,c('chrom', 'maploc', recepientSamples[r],
      donorSample)]
    tempSeg <- subset(allKCseg, 
      samplelist=c(recepientSamples[r], donorSample))
    
    resultKC <- cbind(resultKC, 
        deltaLinearSeg(comb=c(recepientSamples[r], donorSample), 
        tempKC, tempSeg, thres=.2))
    colnames(resultKC)[ncol(resultKC)] <- names(recepientSamples)[r]

  }

  resultList$delta <- resultKC

  diffList[[paste('T', t, sep='')]] <- resultList

}

a <- diffList[['T2']]
par(mfrow=c(3,1))
plotRawCghDotPlot(KCdataSet=a$donor, mirrorLocs=altMirrorLocs, doFilter=T, samples=1, setcex=10)
plotRawCghDotPlot(KCdataSet=a$recepients, mirrorLocs=altMirrorLocs, doFilter=T, samples=2, setcex=10)
plotRawCghDotPlot(KCdataSet=a$delta, mirrorLocs=altMirrorLocs, doFilter=T, samples=2, setcex=10)


source("http://www.bioconductor.org/biocLite.R")
biocLite("CGHnormaliter", lib='~/lib/R')