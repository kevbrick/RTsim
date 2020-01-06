### Camerini Group functions
### KB June 29 2018
#############################

scriptFolder    <- paste0(Sys.getenv('SHARE'),'/code/R/')
imgOutputFolder <- paste0(Sys.getenv('RFIGS'),'/')

## Fix double slashes in foldernames
scriptFolder    <- gsub('//','/',scriptFolder)
imgOutputFolder <- gsub('//','/',imgOutputFolder)

## Load bioconductor (## NOT USED CURRENTLY)
#source('http://bioconductor.org/biocLite.R')

## Load libraries
library(plyr)
library(dplyr)
library(data.table)
library(extrafont)
library(factoextra)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(grid)
library(gridExtra)
library(numform) ## Format numbers nicely
library(preprocessCore)
library(reshape2)
library(tictoc)

############################# 
### standardizeMNSD
### KB June 29 2018
### Normlize the rows of a data frame by mean & sd
### ARGS: 
# z     numeric data frame to normalize
### OUTPUTS:
# Normalized data frame
standardizeMNSD <- function(z) {
  if (is.null(dim(z))){
    rv <- (z-mean(z))/sd(z)
  }else{
    rowmean <- apply(z, 1, mean)
    rowsd <- apply(z, 1, sd)  
    
    rv <- sweep(z, 1, rowmean,"-")  #subtracting mean expression
    rv <- sweep(rv, 1, rowsd, "/")  # dividing by standard deviation
  }
  return(rv)
}

############################# 
### standardizeMedMad
### KB June 29 2018
### Normlize the rows of a data frame by median & mean absolute deviation
### ARGS: 
# z     numeric data frame to normalize
### OUTPUTS:
# Normalized data frame
standardizeMedMad <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

############################# 
### standardize0to1
### KB June 29 2018
### Normlize the rows of a data frame from 0 to 1
### ARGS: 
# z     numeric data frame to normalize
### OUTPUTS:
# Normalized data frame
standardize0to1 <- function(z) {
  ## For vector
  if (is.null(dim(z))){
    rv <- (z-min(z))/(max(z)-min(z))
  }else{ ## For DF / matrix
    rowmin <- apply(z, 1, min)
    rowmax <- apply(z, 1, max)  
    
    rv <- sweep(z, 1, rowmin,"-")  
    rv <- sweep(rv, 1, rowmax-rowmin, "/")  
  }
  return(rv)
}

############################# 
### theme7point
### KB June 29 2018
### Set the default theme to 7-point font
### suitable for publication
### ARGS: none
### OUTPUTS: Sets the theme
theme7point <- function() {
  theme_base_size <- 7
  theme_set(theme_bw(base_size = theme_base_size) %+replace% 
              theme(axis.text = element_text(size=theme_base_size,
                                             color='black'),
                    panel.grid = element_blank(),
                    panel.border=element_blank(),
                    axis.line=element_line(size=.2),
                    axis.ticks = element_line(size=.2)))
}

############################# 
### theme7point
### KB June 29 2018
### Set the default theme to 7-point font
### suitable for publication
### ARGS: none
### OUTPUTS: Sets the theme
theme28point <- function() {
  theme_base_size <- 28
  theme_set(theme_bw(base_size = theme_base_size) %+replace% 
              theme(axis.text = element_text(size=theme_base_size,
                                             color='black'),
                    panel.grid = element_blank(),
                    panel.border=element_blank(),
                    axis.line=element_line(size=.8),
                    axis.ticks = element_line(size=.8)))
}
############################# 
### getImgName
### KB June 29 2018
### Create a systematic image name with date and user
### ARGS: 
# fname         Output file name stem
# type          PNG/PDF/GIF/TIF (default: PDF)
# saveLocation  Folder (default uses system RFIGS path)
### OUTPUTS: File name
getIMGname <- function(fname=NULL,
                       type="PDF", ## 
                       saveLocation=NULL) {
  
  systemUserName <- as.character(Sys.info()[7])
  
  if (is.null(saveLocation)){
    saveLocation <- imgOutputFolder
  }
  
  if (type == 'PDF'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".pdf"))}
  if (type == 'PNG'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".png"))}
  if (type == 'GIF'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".gif"))}
  if (type == 'TIF'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".tiff"))}
}

############################# 
### ggCorMat
### KB July 03 2018
### Generate a ggplot correlation matrix
### ARGS: 
# mCC           correlation matrix
# newOrd1       Order of field 1
# newOrd2       Order of field 2
# flipIt        Flip matrix
# numOFF        Do NOT print CCs
# tileFontScale Scale tile text font (defauilt=1)
# varColor      Use different colors for CCs > 0.5 
# tileFontSize Like it says ... default = 7
### RETURNS: ggplot Geom
ggCorMat <- function(mCC,
                     newOrd1=NULL,
                     newOrd2=NULL,
                     flipIt=FALSE,
                     numOFF=FALSE,
                     scalePC=1,
                     tileFontScale=1,
                     varColor=FALSE,
                     tileFontSize=7,
                     keepLeadingZeros=FALSE,
                     decimalPlaces=1,
                     asPercentage=FALSE){
  library(ggplot2)
  library(reshape2)
  
  ## Set default order
  if (is.null(newOrd1)){
    newOrd1 <- names(mCC)
  }

  if (is.null(newOrd2)){
    newOrd2 <- names(mCC)
  }
  
  mCC <- round(mCC,2)
  melted_cormat <- melt(mCC)
  
  # Get lower triangle of the correlation matrix
  get_upper_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  # Get upper triangle of the correlation matrix
  get_lower_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  if (flipIt){
    lower_tri     <- get_lower_tri(mCC)
    melted_cormat <- melt(lower_tri)
  }else{
    lower_tri     <- get_upper_tri(mCC)
    melted_cormat <- melt(lower_tri)
  }
  
  
  if (length(newOrd1)>0){
    melted_cormat$Var1 <- factor(melted_cormat$Var1,newOrd1)
    if (length(newOrd2)>0){
      melted_cormat$Var2 <- factor(melted_cormat$Var2,newOrd2)
    }else{
      melted_cormat$Var2 <- factor(melted_cormat$Var2,newOrd1)
    }
  }
  
  ## Pad with required number of zeros
  if (asPercentage){
    melted_cormat$numz <- f_num(melted_cormat$value*100, digits=decimalPlaces)      
  }else{
    melted_cormat$numz <- f_num(melted_cormat$value ,digits=decimalPlaces)  
  }

    ## Remove leading zeros
  if (keepLeadingZeros){
    melted_cormat$numz <- f_pad_zero(melted_cormat$numz,
                                     width=(decimalPlaces+2))      
  }
  
  colScale <- colorRampPalette(c("green", "blue", "black", "#FF8811", "red"))
  
  if (varColor){
    melted_cormat$txtColor <- 0
    melted_cormat$txtColor[melted_cormat$value > 0.5] <- 1
  }else{
    melted_cormat$txtColor <- 0
  }
  
  if (asPercentage){
    myLbl <- 'Overlap (%)'
  }else{
    myLbl <- substitute(paste('Spearman ', R^"2"))
  }
  
  #melted_cormat <- mc
  
  ggheatmap <-  ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white",size=.3) +
    scale_fill_gradient2(low = "grey90", mid = 'orange', high = "orangered3", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name=myLbl,na.value='white') +
    theme(axis.text.x = element_text(vjust = 1, 
                                     hjust = 1,
                                     angle = 45,
                                     size = tileFontSize),
          axis.text.y = element_text(size = tileFontSize)) 
  coord_fixed()
  #    theme(axis.text = element_text(size = tileFontSize)) +
  
  gCInit <- ggheatmap + 
    theme(
      axis.title.x         = element_blank(),
      axis.title.y         = element_blank(),
      panel.grid.major     = element_blank(),
      panel.border         = element_blank(),
      panel.background     = element_blank(),      
      axis.ticks           = element_blank(),
      legend.position      = 'none',
      legend.justification = c(1,1),
      legend.text          = element_text(size = tileFontSize),
      legend.title         = element_text(size = tileFontSize),
      legend.key.size      = unit(.5,'cm'),
      legend.direction     = "vertical")
  #legend.key.size = unit(1*scalePC,'cm'),
  #legend.key.width = unit(0.6*scalePC,'cm'))
  
  if (numOFF){
    gCC <- gCInit
  }else{
    ## NOTE: Font sizes in geom_text are 14/5 X bigger than in theme 
    ## NO idea why !!!
    
    gCC <- gCInit + 
      geom_text(aes(Var2, Var1, label = numz, color = txtColor), 
                size = 5/14*tileFontSize,
                show_guide = FALSE) + 
      scale_color_gradient2(low = "black", mid = 'black', high = "white", 
                            midpoint = 0.1, limit = c(0,1))
  }

  return(gCC)
}

############################# 
### toTPM
### KB July 03 2018
### Convert a vector of strengths into a Tags (or Fragments) per million value
### ARGS: 
# x           vector of values
# noNeg       remove negative values by adding the minimum FPM to all 
## RETURNS: Normalized vector in T(F)PM
toTPM <- function (x,noNeg=FALSE){
  #x <- x+abs(min(x))+1;
  v <- (x/sum(x)*1000000)
  if (noNeg & min(v) < 0){
    v <- v - min(v) + 1
  }
  return(v)
}

############################# 
### convertToQuantiles
### KB July 03 2018
### Convert a vector of values into evenly size quantiles
### ARGS: 
# pData       vector of values
# nQ          number of quantiles (default = 5)
# plotMe      return a diagnostic plot
# revLbl      reverse order (default == smallest = 1)
# labelMaxMin change label of max and min bins to "max" and "min"
# maxLabel    max Label
# minLabel    min Label
## RETURNS: a factorized vector of quantiles
 
convertToQuantiles <- function (pData,nQ=5,plotMe=FALSE,revLbl=FALSE,
                                labelMaxMin=FALSE,
                                maxLabel=NULL,
                                minLabel=NULL){
  library(lsr)
  library(ggplot2)
  
  qvec <- quantile(jitter(pData), 
                   probs=seq(0,1,length=(nQ+1)), 
                   type=1,
                   na.rm=T)
  
  # Old alternate (KB)
  #asdc <- function(q){max(1,which(qvec<=q))}
  #qvals = apply(matrix(pData,nrow=1),2,asdc)
  
  if (revLbl){
    qvals <- cut(pData, breaks = qvec,include.lowest = TRUE, labels = nQ:1)
    qvals <- quantileCut(jitter(pData,0.001),nQ,labels=nQ:1)
  }else{
    qvals <- cut(pData, breaks = qvec,include.lowest = TRUE, labels = 1:nQ)
    qvals <- quantileCut(jitter(pData,0.001),nQ,labels=1:nQ)
  }
  
  
  
  if (plotMe){
    df <- data.frame("pD"=pData,"Q"=qvals)
    g <- ggplot(df,aes(x=Q,y=..count..,color=Q)) + geom_bar() 
    + theme(text=element_text(size=22)) 
    + xlab("Quantile") 
    + ylab("Count")
    print(g)
  }
  
  if (labelMaxMin){
    if (is.null(maxLabel)){maxLabel = 'Max'}
    if (is.null(minLabel)){minLabel = 'Min'}
    
    qV            <- qvals
    qV            <- factor(qV,levels=c(minLabel,2:(nQ-1),maxLabel))
    
    qV[qvals==nQ] <- maxLabel
    qV[qvals==1]  <- minLabel
    qvals <- qV
  }
  
  return(qvals)
}

############################# 
### getLegend
### KB July 10 2018
### get the legend of a ggplot
### ARGS: 
# g   grob
## RETURNS: a legend as a grob
getGGLegend<-function(g){
  tmp <- ggplot_gtable(ggplot_build(g))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

############################# 
### saveImage
### KB July 10 2018
### Save a plot as a PNG and PDF
### ARGS: 
# imgName      name of image
# img2Output   plot
# nH           height
# nW           width
## RETURNS: Nothing
saveImage <- function(imgName, img2Output, nH, nW){
  outPNG <- getIMGname(imgName,'PNG')
  outPDF <- getIMGname(imgName,'PDF')
  
  cat ('************* ERROR ***************\n')
  cat ("saveImage function NOT WORKING !!!!\n")
  cat ("***********************************\n")
  #### Currently NOT working
  return()
  
  ## Closes active graphics devices
  graphics.off()
  
  png(filename = outPNG, width = nH, height = nW, units='in', res=300)
  img2Output
  dev.off()
  
  pdf(file = outPDF, width = nH, height = nW)
  img2Output
  dev.off()
  
  print(paste0('Saved as ',outPNG))
  print(paste0('Saved as ',outPDF))
}

############################# 
### fancy_scientific
### KB July 10 2018
### Change number to scientific notation format:
### Most useful for scale_x[y]_log10(labels=fancy_scientific)
### ARGS: 
# l   number
## RETURNS: Formatted number as expression
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # remove +s
  l <- gsub("\\+", "", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

############################# 
### addThemeToStatCompare
### KB July 10 2018
### Fix font size and line width in stat_compare
### ARGS: 
# f   ggplot figure
## RETURNS: ggplot figure
addThemeToStatCompare <- function(f,fLineSize=0.2,fFontSize=7) {
  f$layers[[which_layers(f, "GeomSignif")]]$aes_params$size <- fLineSize
  f$layers[[which_layers(f, "GeomSignif")]]$aes_params$textsize <- fFontSize*5/14
  return(f)
}

############################# 
### darkenColor
### KB July 17 2018
### Darken a color
### ARGS: 
# color   color
# X       how much darker? (default=1.4x)
## RETURNS: color
darkenColor <- function(color, X=1.4){
  col <- col2rgb(color)
  col <- col/X
  col <- rgb(t(col), maxColorValue=255)
  col
}

############################# 
### lightenColor
### KB July 17 2018
### Lighten a color
### ARGS: 
# color   color
# X       how much lighter? (default=1.4x)
## RETURNS: color
lightenColor <- function(color, X=1.4){
  col <- col2rgb(color)
  col <- col*X
  col <- rgb(t(col), maxColorValue=255)
  col
}

############################# 
### plotHSstrength
### KB July 23 2018
### Make a density scatterplot of hotspot strength in two samples
### ARGS: 
# sA      vector of values for x-axis
# sB      vector of values for y-axis
# nameA   Name for x-axis
# nameB   Name for y-axis
# 
## RETURNS: geom
plotHSstrength <- function (sA,sB,nameA,nameB,
                            noGradient = FALSE,
                            sOK        = NULL,
                            lblAsIs    = FALSE,
                            topLeft    = FALSE,
                            xLim       = 0,
                            yLim       = 0,
                            stdStyle   = TRUE,
                            pointSz    = 2,
                            txtSz      = 7,
                            noLeg      = FALSE){
  
  sDataInit <- data.frame('tpmA'=toTPM(sA),'tpmB'=toTPM(sB))  
  
  if (length(sOK) > 0){
    sData <- sDataInit[sOK>0,]
  }else{
    sData <- sDataInit
  }
  
  ## Get max / min coords
  xMax <- max(sData$tpmA)
  xMin <- min(sData$tpmA)
  
  yMax <- max(sData$tpmB)
  yMin <- min(sData$tpmB)
  
  ## Set limits if we're using them
  if (xLim){
    xMax <- xLim
  }
  
  if (yLim){
    yMax <- yLim
  }
  
  sCC <- round(cor(sData[,1:2],method='spearman')[1,2]^2,2)
  
  if (lblAsIs){
    xLblName <- nameA
    yLblName <- nameB
  }else{
    xLblName <- paste0(nameA,' ','(FPM)')
    yLblName <- paste0(nameB,' ','(FPM)')
  }
  
  if (stdStyle){
    gInit <- ggplot(data=sData,aes(x=tpmA,y=tpmB)) + 
      geom_point(color='grey50',size=pointSz) +
      scale_x_log10(labels=fancy_scientific) + 
      scale_y_log10(labels=fancy_scientific) + 
      geom_smooth(method=lm,linetype=2,colour="NA",se=F) + 
      guides(alpha="none") + 
      annotation_logticks(sides='lb',
                          size=.2,
                          short=unit(0.050,'cm'),
                          mid=unit(0.075,'cm'),
                          long=unit(0.100,'cm')) +
      coord_cartesian(xlim=c(1,max(xMax,yMax)),
                      ylim=c(1,max(xMax,yMax))) +
      xlab(xLblName) + 
      ylab(yLblName) + 
      theme(legend.position=c(1,0),
            legend.justification=c(1,0)) 
  }else{
    gInit <- ggplot(data=sData,aes(x=tpmA,y=tpmB)) + 
      geom_point(color='grey50',size=pointSz) +
      scale_x_log10(labels=fancy_scientific) + 
      scale_y_log10(labels=fancy_scientific) + 
      geom_smooth(method=lm,linetype=2,colour="NA",se=F) + 
      guides(alpha="none") + 
      annotation_logticks(sides='lb') +
      coord_cartesian(xlim=c(1,max(xMax,yMax)),ylim=c(1,max(xMax,yMax))) +
      theme_MF() + 
      xlab(xLblName) + 
      ylab(yLblName) + 
      theme(plot.title=element_text(size=20),
            legend.position=c(1,0),
            legend.justification=c(1,0),
            legend.title=element_text(size=15)) 
  }
  #ggtitle(substitute(paste('Spearman ', R^"2"," = ", cc ,sep=''),list(cc = sCC)))
  
  if (noGradient){
    g <- gInit
  }else{
    g <- gInit + stat_density2d(aes(fill=..level..,
                                    alpha=..level..),
                                geom='polygon',
                                colour='NA')  
    #scale_fill_continuous(low="grey30",high="red") 
  }  
  
  g$labels$fill <- "Hotspot density"
  
  myLbl <- substitute(paste('Spearman ', R^"2"," = ", cc ,sep=''),list(cc = sCC));
  myLbl <- paste("R^2 == ", sCC)
  labelSize <- 10
  if (topLeft){
    if (stdStyle){
      
      thisTheme <- theme_get()
      txtSz <- thisTheme$axis.text$size
      
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = -Inf, y = Inf, 
                         colour = "black", hjust=0,
                         parse=TRUE,
                         size = 5/14*txtSz) +
        theme(axis.line = element_line(color='black'),
              legend.background=element_blank())
    }else{
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = -Inf, y = Inf, 
                         size = labelSize, 
                         colour = "black", hjust=0,
                         parse=TRUE) +
        theme(axis.line = element_line(color='black'),
              legend.background=element_blank())
    }
  }else{
    if (stdStyle){
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = 1, y = yMax-(yMax*.1),  
                         colour = "black", hjust=0,
                         parse=TRUE) +
        theme(axis.line = element_line(size=.5),
              legend.background=element_blank())
    }else{
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = 1, y = yMax-(yMax*.1),  
                         size = labelSize, 
                         colour = "black", hjust=0,
                         parse=TRUE) +
        theme(axis.line = element_line(size=.5),
              legend.background=element_blank())
    }
  }
  
  if (noLeg){
    gg <- gg + theme(legend.position='none')
  }
  
  return(gg)
}
