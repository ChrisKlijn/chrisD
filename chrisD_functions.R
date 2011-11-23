# chrisD_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: chrisD analyses
# Description: Functions used in the analysis of chrisD data
# -------------------------------------------------------------------

plotMetMat <- function (metMat, KC.Seg, plotTitle=NA) {

  # Function to plot a matrix of probe values with 0, 1 and -1 for
  # significant delta value yes/no
  # Also provide KC.Seg to calculate the probe-based chromosome ends

  probeChromEnds <- cumsum(tapply(KC.Seg$maploc, 
    KC.Seg$chrom, length))

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

  # Set colors, mid is the background color, low is the color of the
  # negative differences and high is the color of the positive 
  # differences
  gainLossCols <- colorpanel(3, low='green', mid='black', high='orange')
  closeCols <- colorpanel(2, low='blue', high='red')
  tumorCols <- colorpanel(length(unique(orderFrame$tumorID)), 
    low=colors()[1], high=colors()[24])
  names(tumorCols) <- unique(orderFrame$tumorID)

  # intialize plot
  plot(0,0, xlim=c(0, nrow(metMat)+30000), ylim=c(0, ncol(metMat)+1),
    type='n',axes=F, xlab=NA, ylab=NA, main=plotTitle)

  # Use run length encoding to draw probe-index based rectangle
  # use additional rectangles for tumor site and tumor ID
  for (j in 1:ncol(metMat)) {
    tumVect <- rle(metMat[,j])
    start <- 0
    for (i in 1:length(tumVect$length)) {
      rect(start, j-1, tumVect$length[i]+start, j, 
        col=gainLossCols[tumVect$value[i]+2], border=NA)
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
  legend('topleft', fill=c(gainLossCols, closeCols), cex=.7,
    legend=c('loss','none','gain','localMeta', 'distantMeta'), horiz=T)
}

plotCorMap <- function (dataMat, classLabels, labs, plotTitle='Temp',
  plotCol=NULL) {
  
  require(gplots) # For colorpanel and heatmap.2
     
  corMat <- cor(dataMat, use='na.or.complete')

  if (is.null(plotCol)) {
    plotCol <- colorpanel(n=length(unique(classLabels)),
      low=colors()[24],
      high=colors()[349])
  }
  
  colIndex <- as.numeric(as.factor(classLabels))

  heatmap.2(corMat, scale='none', trace='none',
    labRow=labs, labCol=labs, 
    col=colorpanel(265, low='blue', high='yellow'), margins=c(10,10),
    density.info='none',
    ColSideColors=plotCol[colIndex],
    RowSideColors=plotCol[colIndex],
    main=plotTitle)


}
