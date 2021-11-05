## simulateReplication
library(gganimate)
library(ggplot2)
library(zoo)
library(animation)
library(gganimate)
library(Metrics)
library(plyr)

######################## DEFINE FUNCTIONS ############################
######################################################################
### Generate a randomString
genRandomString <- function(n = 1, len = 10){
  a <- do.call(paste0, replicate(len, sample(LETTERS, n, TRUE), FALSE))
  rString <- paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
  substr(rString,1,len)
}
######################################################################
### Root mean square Error Function
getRMSE = function(m, o){
  NAok <- !is.na(m*o)
  NAok <- NAok & (m*o > 0)
  #retRMSE <- sqrt(mean((m[NAok] - o[NAok])^2))
  retRMSE <- Metrics:::rmse(m[NAok],o[NAok])
  return(retRMSE)
}

######################################################################
### Normalize -1 to 1
normalize1to1 <- function(x){((x-quantile(x,.01,na.rm=TRUE))/quantile((x-quantile(x,0.01,na.rm=TRUE)),0.99,na.rm=TRUE)*2)-1}

######################################################################
### Sum square Error Function
getSSE = function(m, o){
  NAok   <- !is.na(m*o)
  NAok <- NAok & (m*o > 0)
  #retSSE <- sum((m[NAok] - o[NAok])^2)
  retSSE <- Metrics:::sse(m[NAok],o[NAok])
  return(retSSE)
}

######################################################################
### Mean Absolute Square Error Function
getMASE = function(m, o){
  NAok <- !is.na(m*o)
  NAok <- NAok & (m*o > 0)
  retMASE <- Metrics:::mase(m[NAok],o[NAok],1)
  return(retMASE)
}

######################################################################
### Function to randomly subset origins
getRandomOrigins <- function (inDF,numForks, useStrength = TRUE, seed = NULL){
  
  if (!useStrength){
    inDF$strength <- 1
  }
  
  if (dim(inDF)[1] == numForks){
    return(inDF$rFrom)
  }
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  return(sample(inDF$rFrom,
                size=numForks, 
                replace=FALSE, 
                prob=inDF$strength))
}

######################################################################
### Function to calculate read enrichment and smooth with sliding window
smoothAndLog2Coverage <- function(dModel,winSz=25,stepSz=1,applySmoothing=FALSE,alreadyLog2=FALSE){
  
  if (alreadyLog2){
    l2 <- dModel
  }else{
    l2 <- log2(dModel/mean(dModel,na.rm=TRUE))
  }
  
  l2[is.infinite(l2)] <- 0
  l2[is.na(l2)] <- 0
  
  ## Smoothing
  if (!applySmoothing){
    return(l2)
  }else{
    smoothedData <- rollapply(data=l2, 
                              width = winSz, 
                              by = stepSz, 
                              FUN = mean, 
                              align = "center")
    
    return(c(rep(0,((winSz-1)/2)),smoothedData,rep(0,((winSz-1)/2))))
  }
  
  
}

'%!in%' <- function(x,y)!('%in%'(x,y))

######################################################################
## Get chromosome sizes
getChromSize <- function(myGenome,
                         myCS,
                         noSexCS=FALSE,
                         noX=FALSE,
                         noY=FALSE,
                         noM=TRUE){
  
  allCS <- c(paste0('chr',c(1:50)),'chrX','chrY','chrM')
  
  if (noSexCS){allCS <- allCS[allCS %!in% c('chrX','chrY','chrZ')]}
  if (noM){allCS <- allCS[allCS %!in% 'chrM']}
  if (noX){allCS <- allCS[allCS %!in% 'chrX']}
  if (noY){allCS <- allCS[allCS %!in% 'chrY']}
  
  genomeSz         <- read.table(paste0(Sys.getenv('GENOMES'),'/',myGenome,'/genome.fa.fai'),header=FALSE)
  names(genomeSz)  <- c('cs','size','cum','xa','xb')
  
  genomeSz$size    <- as.numeric(genomeSz$size)
  
  genomeSz$cs      <- factor(genomeSz$cs,levels=allCS)
  
  genomeSz         <- genomeSz[!is.na(genomeSz$cs),]
  
  genomeSz         <- genomeSz[order(genomeSz$cs),]
  
  genomeSz$min     <- c(0,cumsum(genomeSz$size[1:(length(genomeSz$size)-1)]))+1
  genomeSz$max     <- cumsum(genomeSz$size)
  
  genomeSz <- genomeSz[,c('cs','size','min','max')]
  
  ## Remove unused chromosomes
  z <- plyr:::join(data.frame(cs=levels(genomeSz$cs)),genomeSz,by='cs')  
  genomeSz$cs <- factor(genomeSz$cs, levels=z$cs[!is.na(z$size)])
  
  if (myCS %in% c('All','all','ALL')){
    return(genomeSz)
  }else{
    return(genomeSz$size[genomeSz$cs == myCS])
  }
}

######################################################################
## Import origins and select by CS
getOrigins <- function(oriBG,myCS,stepSz){
  
  halfSz       <- stepSz/2;
  
  gOri         <- read.table(oriBG,header=FALSE)
  names(gOri)  <- c("cs","from","to","strength")
  
  ## Remove Oris with -ve strength
  gOri         <- gOri[gOri$strength>0,] 
  
  gOriCS       <- subset(gOri, gOri$cs %in% myCS)
  
  gOriCS$mid   <- round((gOriCS$to + gOriCS$from)/2)
  
  ##KB oriCS$rFrom <- ceiling((oriCS$mid-5000)/10000)*10000 - 4999
  ##KB oriCS$rTo   <- ceiling((oriCS$mid-5000)/10000)*10000 + 5000
  gOriCS$rFrom <- ceiling((gOriCS$mid-halfSz)/stepSz)*stepSz - (halfSz-1)
  gOriCS$rTo   <- ceiling((gOriCS$mid-halfSz)/stepSz)*stepSz + (halfSz)
  
  return(gOriCS)
}

######################################################################
## Convert real coordinates to "bin" coordinates
markBinsByCS <- function(modDF, csDets, stepSz){
  
  halfSz <- stepSz/2
  
  retDF <- data.frame(from = ceiling((csDets$min - halfSz)/stepSz)*stepSz - (halfSz - 1),
                      to   = ceiling((csDets$max - halfSz)/stepSz)*stepSz + halfSz)
  
  
  
  return(retDF)
}
######################################################################
## Import origins and select by CS
getALLOrigins <- function(oriBG,csDets,stepSz){
  
  halfSz        <- stepSz / 2;
  
  gOri          <- read.table(oriBG,header=FALSE)
  names(gOri)   <- c("cs","from","to","strength")
  
  ## Remove Oris with -ve strength
  gOri          <- gOri[gOri$strength>0,] 
  
  gOriCS        <- plyr::join(gOri, csDets, by = 'cs', type = 'inner')
  
  gOriCS$from   <- gOriCS$from + gOriCS$min - 1
  gOriCS$to     <- gOriCS$to   + gOriCS$min - 1
  
  gOriCS$mid    <- round((gOriCS$to + gOriCS$from)/2)
  
  ##KB oriCS$rFrom <- ceiling((oriCS$mid-5000)/10000)*10000 - 4999
  ##KB oriCS$rTo   <- ceiling((oriCS$mid-5000)/10000)*10000 + 5000
  gOriCS$rFrom  <- ceiling((gOriCS$mid-halfSz)/stepSz)*stepSz - (halfSz-1)
  gOriCS$rTo    <- ceiling((gOriCS$mid-halfSz)/stepSz)*stepSz + (halfSz)
  
  return(gOriCS)
}

######################################################################
## import timing data
getRT <- function(inRT,myCS,myModel,lSmooth,lLog2,inRTdata=NULL){
  
  if (is.null(inRTdata)){
    RTdata        <- read.table(inRT,header=FALSE)
    names(RTdata) <- c("cs","from","to","cover")
  }else{
    RTdata <- inRTdata
  }
  
  RTcs <- subset(RTdata,RTdata$cs==myCS)[2:4]
  
  ## Lift sparse replication timing bedgraph to full representation
  if (sum(is.na(match(RTcs$from,table = myModel$from-1))) == 0){
    print ('Option A')
    myModel$rawstrength[match(RTcs$from,table = myModel$from - 1)] <- RTcs$cover  
  }else{
    print ('Option B')
    myModel$rawstrength[match(RTcs$from,table = myModel$from)] <- RTcs$cover
  }
  
  rawRT <- myModel$rawstrength
  retRT <- smoothAndLog2Coverage(myModel$rawstrength,
                                 applySmoothing = lSmooth,
                                 alreadyLog2 = lLog2)
  
  return(list(raw=rawRT, smooth=retRT))
}

######################################################################
## import timing data
getALLRT <- function(inRT,csDets,myModel,lSmooth,lLog2,inRTdata=NULL,stepSz=10000){
  
  halfSz <- stepSz/2;
  
  if (is.null(inRTdata)){
    RTdata        <- read.table(inRT,header=FALSE)
    names(RTdata) <- c("cs","from","to","cover")
  }else{
    RTdata <- inRTdata
  }
  
  RTdata <- subset(RTdata,RTdata$cs %in% unique(csDets$cs))
  
  RTdata$cs <- factor(RTdata$cs, 
                      levels=levels(csDets$cs))
  
  RTdata    <- RTdata[order(RTdata$cs,RTdata$from),]
  
  RTdata      <- plyr::join(x = RTdata, y = csDets, 
                            by = 'cs', 
                            type = 'inner')
  
  RTdata$binMin                  <- ceiling((RTdata$min - halfSz)/stepSz)*stepSz + 1
  #RTdata$binMin                  <- ceiling((RTdata$min)/stepSz)*stepSz - (halfSz - 1)
  RTdata$binMin[RTdata$binMin<0] <- 1
  
  RTdata$initFrom <- RTdata$from
  RTdata$initTo   <- RTdata$to
  
  RTdata$from  <- RTdata$from + RTdata$binMin 
  RTdata$to    <- RTdata$to   + RTdata$binMin 
  
  RTdata$cs    <- factor(RTdata$cs,levels(csDets$cs))
  
  #RTdata       <- RTdata[!is.na(RTdata$cs),2:4]
  #return(RTdata)
  #print(head(RTdata))
  
  ## Lift sparse replication timing bedgraph to full representation
  if (sum(is.na(match(RTdata$from,table = myModel$from))) == 0){
    print ('Option A')
    myModel$rawstrength[match(RTdata$from,table = myModel$from )] <- RTdata$cover  
  }else{
    print ('Option B')
    myModel$rawstrength[match(RTdata$from,table = myModel$from)] <- RTdata$cover
  }
  
  rawRT <- myModel$rawstrength
  retRT <- smoothAndLog2Coverage(myModel$rawstrength,
                                 applySmoothing = lSmooth,
                                 alreadyLog2 = lLog2)
  
  return(list(raw=rawRT, smooth=retRT))
}

######################################################################
## import hotspot data
getDSBs <- function(myGenome,myCS,stepSz){
  
  halfSz <- stepSz/2
  
  if (myGenome == 'mm10'){
    gDSBs         <- read.table('dsbs/B6_hotspots.SSDS_V_AffySeqPk1.bedgraph',header=FALSE)
    names(gDSBs)  <- c("cs","from","to","strength","affySeq")
  }
  
  if (myGenome == 'hg38'){
    gDSBs         <- read.table('dsbs/AA1_all.HSstrength.bedgraph',header=FALSE)
    #names(dsb)  <- c("cs","from","to","strength")
    names(gDSBs) <- c("cs","from","to","strength")
    gDSBs$affySeq  <- 1
  }
  
  retDSB      <- subset(gDSBs, gDSBs$cs==myCS)
  
  retDSB$mid  <- round((retDSB$to + retDSB$from)/2)
  ##KB dsbCS$rFrom <- ceiling((dsbCS$mid-5000)/10000)*10000 - 4999
  ##KB dsbCS$rTo   <- ceiling((dsbCS$mid-5000)/10000)*10000 + 5000
  
  retDSB$rFrom <- ceiling((retDSB$mid - halfSz)/stepSz)*stepSz- (halfSz - 1)
  retDSB$rTo   <- ceiling((retDSB$mid - halfSz)/stepSz)*stepSz + halfSz
  
  return(retDSB)
}

######################################################################
## import hotspot data
getALLDSBs <- function(myGenome, csDets, stepSz, myDSBs){
  
  halfSz <- stepSz/2
  
  gDSBs <- read.table(myDSBs, header=FALSE)  
  if (dim(gDSBs)[2] == 5){
    names(gDSBs)  <- c("cs","from","to","strength","affySeq")
  }else{
    names(gDSBs)  <- c("cs","from","to","strength")
  }
  # if (myGenome == 'mm10'){
  #   gDSBs         <- read.table('B6_hotspots.SSDS_V_AffySeqPk1.bedgraph',header=FALSE)
  #   names(gDSBs)  <- c("cs","from","to","strength","affySeq")
  #   }
  # 
  # if (myGenome == 'hg38'){
  #   gDSBs         <- read.table('dsbs/AA1_all.HSstrength.bedgraph',header=FALSE)
  #   #names(dsb)  <- c("cs","from","to","strength")
  #   names(gDSBs) <- c("cs","from","to","strength")
  #   gDSBs$affySeq  <- 1
  # }
  
  gDSBs         <- plyr::join(gDSBs, csDets, by = 'cs', type = 'inner')
  
  gDSBs$from   <- gDSBs$from + gDSBs$min - 1
  gDSBs$to     <- gDSBs$to   + gDSBs$min - 1
  
  #retDSB      <- subset(gDSBs, gDSBs$cs==myCS)
  retDSB       <- gDSBs
  
  retDSB$mid  <- round((retDSB$to + retDSB$from)/2)
  ##KB dsbCS$rFrom <- ceiling((dsbCS$mid-5000)/10000)*10000 - 4999
  ##KB dsbCS$rTo   <- ceiling((dsbCS$mid-5000)/10000)*10000 + 5000
  
  retDSB$rFrom <- ceiling((retDSB$mid - halfSz)/stepSz)*stepSz- (halfSz - 1)
  retDSB$rTo   <- ceiling((retDSB$mid - halfSz)/stepSz)*stepSz + halfSz
  
  return(retDSB)
}

######################################################################
## Compare modelled RT to experimental RT
correlationByGroup<- function(dd,col1,col2,grp,cType='pearson'){
  
  #-------------------------
  doCC <- function(xx){
    return(data.frame(R2 = round(cor(xx$x, xx$y, 
                                     use='complete.obs',
                                     method='pearson')^2,3)))
  }
  #-------------------------
  genLbl <- function(x){
    xx <- as.character(x[1]); 
    yy <- round(as.numeric(x[2]),2); 
    return(bquote(.(xx)*": "*R^2*" = "*.(yy)*"  "))
  }
  
  #-------------------------
  
  dF4CC <- data.frame(grp=dd[[grp]],x=dd[[col1]],y=dd[[col2]])
  
  retDF <- ddply(dF4CC, .(grp), doCC)
  names(retDF) <- c(grp,'R2')
  
  retDF$lbl <- apply(retDF,1,genLbl)
  
  return(retDF)
}

######################################################################
## Compare modelled RT to experimental RT
flipHiCByCS<- function(dd){
  
  #-------------------------
  doCC <- function(xx){
    return(data.frame(R2 = round(cor(xx$hiC, xx$RT, 
                                     use='complete.obs',
                                     method='spearman'),3)))
  }
  
  #-------------------------
  
  retDF <- ddply(dd, .(cs), doCC)
  names(retDF) <- c('cs','R2')
  
  return(retDF)
}

######################################################################
## Compare features to RT ; violins
plotRTviolins <- function(t,sample){
  
  sName     <- sub('.mm10',replacement = '',sample)
  t$simU    <- t[[paste0('simRT_',sample)]]
  t$expU    <- t[[paste0('expRT_',sample)]]
  
  t$RTexp   <- t$expU
  t$RTsim   <- -1*t$simU
  t$RTsimL2   <- (log2(abs(t$simU) / mean(abs(t$simU))))*-1
  
  tHS <- t[t$SSDST1 >  0 ,]
  tHS$typ <- 'HS'
  
  rHS <- t[(t$SSDST1 + t$SSDST2) <= 0,]
  rHS$typ <- 'rand'
  
  cHS <- t[t$COYin > 0,]
  cHS$typ <- 'CO'
  
  oHS <- t[t$origins > 0,]
  oHS$typ <- 'ori'
  
  affHS <- t[t$affySeq > 0,]
  affHS$typ <- 'affy'
  
  prHS <- t[t$Prdm9ChIPSeq > 0,]
  prHS$typ <- 'prdm9'
  
  r2HS <- t[sample(1:length(t$cs),nrow(tHS)*50,replace=TRUE),]
  r2HS$typ <- 'randN'
  
  tt <- rbind(oHS,affHS,prHS,tHS,cHS,rHS,r2HS)
  
  tt$typ <- factor(tt$typ,levels=c('rand','randN','ori','affy','prdm9','HS','CO'))
  sCols  <- scale_color_manual(values=c('grey50','grey50','magenta','orange','firebrick','forestgreen','dodgerblue2'))
  sFill  <- scale_fill_manual(values=c('grey50','grey50','magenta','orange','firebrick','forestgreen','dodgerblue2'))
  noLeg  <- theme(legend.position='none')
  lPos   <- theme(legend.position=c(1,1),legend.title=element_blank(),legend.justification=c(1,1))
  
  gE <- ggplot(tt,aes(x=typ,y=RTexp,fill=typ)) + 
    geom_violin(color='grey50',alpha=.3,trim = FALSE,lwd=.3) +
    geom_boxplot(fill='grey95',width=.1,notch=TRUE,lwd=.3,outlier.size=.1)+ 
    lPos +
    theme(legend.position='none') + 
    xlab('') + 
    ylab(paste0('Experimental RT\n',sName)) + 
    geom_hline(yintercept=0,lwd=0.3,lty='dashed') + 
    sCols + sFill
  
  gS <- ggplot(tt,aes(x=typ,y=RTsimL2,fill=typ)) + 
    geom_violin(color='grey50',alpha=.3,trim = FALSE,lwd=.3) + 
    geom_boxplot(fill='grey95',width=.1,notch=TRUE,lwd=.3,outlier.size=.1)+ 
    lPos +
    theme(legend.position='none') + 
    xlab('') + 
    ylab(paste0('Simulated RT)\n',sName)) + 
    geom_hline(yintercept=0,lwd=0.3,lty='dashed') + 
    sCols + sFill
  
  ggsave(getIMGname(paste0('recombViolins_expRT',sample),'PNG'),plot = gE,width=6*.8,height=3*.8)
  ggsave(getIMGname(paste0('recombViolins_simRT',sample),'PNG'),plot = gS,width=6*.8,height=3*.8)
  
  return(lViolins=list(exp=gE,
                       sim=gS))
}

######################################################################
## calculate quantiles by group
quantilesByGroup<- function(dd,grp,col2use,numQuantiles=25){
  
  #-------------------------
  getQ <- function(xx,q = 25){
    return(as.data.frame(convertToQuantiles(xx$x, nQ = q, numericQuantiles=TRUE)))
  }
  #-------------------------
  
  dF4CC <- data.frame(grp=dd[[grp]],x=dd[[col2use]])
  #print(head(dF4CC))
  retDF <- ddply(dF4CC, .(grp), getQ)
  names(retDF) <- c(grp,'q')
  
  return(retDF)
}
######################################################################
## Compare modelled RT to experimental RT
testModelVData <- function(oMod,
                           vRT,
                           nPCrep,
                           nReplicationTime,
                           noPlot=FALSE,
                           normBoth=FALSE,
                           style='new',
                           smoothLevel=10,
                           titleType='External',
                           upColor='#7fbf7b',
                           downColor='#af8dc3'){
  
  modelDF             <- oMod$model
  
  #############################################
  ## Get simulated RT and add requisite noise
  numReplicatingCells <- oMod$params$numCells 
  numNoiseCells2Add   <- round(numReplicatingCells/nPCrep*100) - numReplicatingCells
  numTotalCells       <- numReplicatingCells + numNoiseCells2Add
  
  simNoise            <- rep(2*numNoiseCells2Add,
                             length(modelDF$simByTime[,1])) 
  
  sim2use              <- smoothAndLog2Coverage(modelDF$simByTime[,nReplicationTime] + simNoise,
                                                applySmoothing = FALSE)
  
  #############################################
  ## Add strength and sim fields to model object
  modelDF$strength <- vRT
  modelDF$sim2use  <- sim2use
  
  ## FOR BLACKOUT REGIONS
  modelDF$preNormStr <- modelDF$strength
  modelDF$preNormSim <- modelDF$sim2use
  
  #############################################
  ## Normalize to account for median shifting
  if (normBoth){
    modelDF$strength       <- standardizeMNSD(modelDF$strength)
    modelDF$sim2use         <- standardizeMNSD(modelDF$sim2use)
  }
  
  #############################################
  
  ## Assure we only consider chroms that ARE NOT excluded
  csChk      <- as.data.frame(aggregate(modelDF$ori,by=list(cs=modelDF$whatCS),FUN=sum))
  goodCSList <- as.character(csChk$cs[csChk$x != 0])
  csOK       <- modelDF$whatCS %in% goodCSList
  
  ## Get correlations & other stats
  #print(sum(vRT[csOK]))
  #print(sum(sim2use[csOK]))
  
  nR2   <- round(cor(vRT[csOK], sim2use[csOK], use='complete.obs')^2,3)
  dfR2  <- correlationByGroup(data.frame(rt=vRT,sim=sim2use,cs=modelDF$whatCS),"rt","sim","cs")
  names(dfR2)[1] <- "whatCS"
  
  ## for single CS figs
  if (nrow(dfR2) == 1){
    dfR2$lbl <- dfR2$whatCS 
  }
  
  nRMSE  <- getRMSE(vRT[csOK], sim2use[csOK])
  nSSE   <- getSSE(vRT[csOK] , sim2use[csOK])
  nMASE  <- nRMSE
  
  #############################################
  ## Build title string
  if (oMod$params$useOriStrength){
    useStrengthString <- 'use strength'
  }else{
    useStrengthString <- 'no strength'
  }
  
  if (oMod$params$recycle){
    useRecyclingString <- 'useRec'
  }else{
    useRecyclingString <- 'noRec'
  }
  
  # titleText <- bquote(R^2*" = "*.(round(nR2,2))*"; "*
  #                       .(round(oMod$params$oriPerMb*2,2)) ~ replisomes^Mb*"; "*
  #                       .(nReplicationTime)*" cycles; "*
  #                       .(useStrengthString))
  
  titleText <- bquote(R^2*" = "*.(round(nR2,2))*"; "*
                        .(round(oMod$params$oriPerMb*2,2)) ~ replisomes^Mb)
  #############################################
  ## Plot figure if required
  if (noPlot){
    g <- 'No Plot'
    rName <- 'No Plot'
  }else{
    if (style == 'old'){
      g <- ggplot(modelDF,aes(x=trueCoord/1000000,y=strength)) +
        geom_hline(yintercept=0,color='grey10',alpha=.3,lwd=.2) +
        geom_line(color='grey40',alpha=.2,lwd=.2) +
        geom_line(color='grey20',alpha=.6,lwd=.2) +
        geom_line(aes(x=trueCoord/1000000,y=sim2use),
                  'color'='darkorange1') +
        ylab('Timing') + xlab('Position (Mb)') +
        ggtitle(label = titleText) +
        theme(plot.title = element_text(size = 7),
              panel.border=element_blank(),
              axis.line = element_line(),
              strip.background=element_blank(),
              strip.text=element_blank(),
              panel.grid = element_line(linetype='dotted',size=.3,color=alpha('grey60',alpha = 0.5)),
              legend.position = c(0,1)) + 
        facet_wrap(~whatCS,ncol=2)
    }else{
      #############################################
      ## Get yMax for scaling plots
      rangeData <- c(modelDF$strength,modelDF$sim2use)
      rangeQs   <- quantile(rangeData,probs = c(0.01,0.99),na.rm=TRUE)
      yMax      <- max(abs(rangeQs))
      
      #############################################
      ## Build DF for plot (speeds it up ... )
      plotDF <- data.frame(pos=modelDF$trueCoord/1000000,
                           RT=modelDF$strength,
                           sim=modelDF$sim2use,
                           whatCS=modelDF$whatCS,
                           preNormStr=modelDF$preNormStr,
                           preNormSim=modelDF$preNormSim)
      
      #############################################
      ## Get chromosome start and end coordinates
      r           <- rle(x = as.vector(modelDF$whatCS))
      csEnds      <- cumsum(r$lengths)
      csStarts    <- c(1,csEnds[1:length(csEnds)-1])
      valsToZero  <- sort(c(csEnds,csStarts))
      
      #############################################
      ## Build RT up and RT down segments
      plotDF$RTup                  <- plotDF$RT
      plotDF$RTup[plotDF$RT < 0]   <- 0
      plotDF$RTup[valsToZero]      <- 0
      
      plotDF$RTdown                <- plotDF$RT
      plotDF$RTdown[plotDF$RT > 0] <- 0
      plotDF$RTdown[valsToZero]    <- 0
      
      ## Downsample to plot in a reasonable time !!!
      plotDF <- plotDF[seq(1,dim(plotDF)[1],smoothLevel),]
      
      ## KB Oct 16 2019
      minByCS <- aggregate(1:length(plotDF$whatCS),by=list(cs=plotDF$whatCS),FUN=min)
      maxByCS <- aggregate(1:length(plotDF$whatCS),by=list(cs=plotDF$whatCS),FUN=max)
      
      skipVals <- sort(c(minByCS$x,maxByCS$x))
      
      ## SET UP "BLACKOUT BOXES" : Regions that coincide with GAPS & are not used for assessing fit
      naVals <- which(is.na(plotDF$preNormStr) | 
                        is.na(plotDF$preNormSim) | 
                        abs(plotDF$preNormStr) < 0.001 | 
                        abs(plotDF$preNormSim) < 0.001)
      
      naPos  <- naVals[!(naVals %in% skipVals)]
      
      naX <- data.frame(pos=naPos,
                        xfrom  = plotDF$pos[naPos-1],
                        xto    = plotDF$pos[naPos+1],
                        whatCS = plotDF$whatCS[naPos])
      
      #print(naX)
      
      ## ASSURE PRETTY Y-TICK VALS
      yTickVals <- round(seq(-yMax,yMax,length.out = 5),1)
      yTickVals[1] <- yTickVals[2]*2
      yTickVals[5] <- yTickVals[3]*2
      
      yTickLbl  <- yTickVals
      yTickLbl[c(1,5)]  <- ''
      
      if (normBoth){
        yLabel <- 'Normalized RT'
      }else{
        yLabel <- 'RT'
      }
      
      #print(yMax)
      #print(yTickVals)
      #print(yTickLbl)
      if (length(unique(plotDF$whatCS)) == 1){
        xCS <- gsub(pattern = 'chr','chromosome ',plotDF$whatCS[1])
        
        xLabel <- paste0('Position on ',xCS,' (Mb)')
      }else{
        xLabel <- 'Position (Mb)'
      }
      
      g <- ggplot(plotDF) +
        geom_hline(yintercept=0,color='grey10',alpha=.3,lwd=.2) +
        geom_area(aes(x=pos,y=RTup),fill=upColor,alpha=.9) + 
        geom_area(aes(x=pos,y=RTdown),fill=downColor,alpha=.9) + 
        geom_line(aes(x=pos,y=RT),color='grey40',alpha=.9,lwd=.15) +
        geom_line(aes(x=pos,y=sim),'color'='black',lwd=0.3) +
        geom_rect(data=naX,ymin=-Inf,ymax=Inf,aes(xmin=xfrom,xmax=xto),fill=alpha('white',1)) +
        geom_hline(yintercept=0,color='grey10',alpha=.3,lwd=.2) +
        ylab(yLabel) + 
        xlab(xLabel) +
        theme(plot.title = element_text(size = 7),
              panel.border=element_blank(),
              axis.line = element_line(),
              strip.background=element_blank(),
              strip.text=element_blank(),
              panel.grid.major = element_line(linetype='dashed',size=.1,color='grey80'),
              panel.grid.minor = element_blank(),
              legend.position = c(0,1)) + 
        facet_wrap(~whatCS,ncol=2) +
        coord_cartesian(ylim=c(-yMax,yMax)) + 
        scale_y_continuous(breaks = yTickVals, labels = yTickLbl)
      
      if (titleType == 'External'){
        g <- g + ggtitle(label = titleText) 
      }else{
        g <- g + annotate(geom='text',label=c(titleText),x=-Inf,y=Inf,hjust=-0.05,vjust=1,size=(7*5/14),parse=TRUE)
      }
      
      if (length(unique(plotDF$whatCS)) > 1){
        g <- g + geom_text(data=dfR2,
                           parse=TRUE,
                           aes(x=Inf, y=Inf, label=lbl),
                           hjust=1.2,
                           vjust=1,
                           size=7*5/14,check_overlap = TRUE)
      }
      
    }  
    rName <- paste0('repOriModelling_',oMod$params$genome,'_v_RT_',
                    oMod$params$chrom,
                    '_sim_',numReplicatingCells,'_of_',
                    numTotalCells,'cells_',
                    oMod$params$oriPerMb,'_OriperMb_',
                    nReplicationTime,'min_',
                    oMod$params$repRate,'KbPerMin_',
                    useStrengthString,'_',useRecyclingString,
                    "_RMSE_",round(nRMSE,3))  
    rName <- rName[1]
  }
  
  
  retList <- list(time              = nReplicationTime,
                  numNoiseCells2Add = numNoiseCells2Add,
                  pcRep             = nPCrep,
                  R2                = nR2,
                  MASE              = nMASE,
                  RMSE              = nRMSE,
                  SSE               = nSSE,
                  fig               = g,
                  name              = rName)
}

#################### END FUNCTIONS ###################################
