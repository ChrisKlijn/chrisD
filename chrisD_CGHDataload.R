# chrisD_CGHDataload.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: aCGH data Chris D - matched primary and metastasis
# Description:  Code to load the aCGH data from the Nimblegen
#               source files
# -------------------------------------------------------------------

# Libraries, data and sources

# Run on medoid

setwd("/home/klijn/data/smallproj/chrisD/rawdata")

# Load sampleinformation

sampleInfo <- read.delim('../sampleinfo.csv', stringsAsFactors=F)
# Every row is duplicated in essence, so remove the odd rows
sampleInfo <- sampleInfo[seq(1, nrow(sampleInfo), by=2) ,]
# Order on hybridization number
sampleInfo <- sampleInfo[order(paste(sampleInfo$Slide, sampleInfo$Spot, sep='')),]


# Loading data into R

fileNamesData <- dir(pattern='reduce.txt') 

for (i in 1:length(fileNamesData)) {
	cat(fileNamesData[i], '\n')
	tempFrame <- read.delim(fileNamesData[i], stringsAsFactors=F)
	if (i == 1) {
		resultFrame <- tempFrame
		colnames(resultFrame)[4] <- fileNamesData[i]
	}
	resultFrame[,i+3] <- tempFrame[,4]
	colnames(resultFrame)[i+3] <- fileNamesData[i]
}

colnames(resultFrame) <- gsub('_CMFS01_segMNT_reduce.txt', '', colnames(resultFrame))

# Make a KC frame for plotting of the samples
allKC <- resultFrame[,2:ncol(resultFrame)]
rownames(allKC) <- resultFrame[,1]
colnames(allKC) <- c('chrom', 'maploc', colnames(allKC)[3:ncol(allKC)])
allKC <- allKC[-which(allKC$chrom == 'chrM' | allKC$chrom == 'chrY'),]
allKC$chrom <- gsub('chr', '', allKC$chrom)
allKC$chrom <- gsub('X', '20', allKC$chrom)
allKC$chrom <- gsub('Y', '21', allKC$chrom)
allKC$chrom <- as.numeric(allKC$chrom)

save(file='../rawData_chrisD.Rda', list=c('allKC', 'sampleInfo'))
