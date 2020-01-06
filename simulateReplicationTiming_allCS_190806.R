source(paste0(Sys.getenv('RTSCRIPTS'),'/RDCOstandardfunctions.R'))
source(paste0(Sys.getenv('RTSCRIPTS'),'/repliSim_loadModules.R'))
#source('repliSim_loadModules.R')
#source('RDCOstandardfunctions.R')

#######################################################
#######################################################
#######################################################
#######################################################
simRT_buildModel <- function (repRate         = 1, ## Kb per minute
                              originsFile     = '',
                              dsbFile         = NULL,
                              timePassed      = 1000, ## minutes
                              oriPerMb        = 0.2,
                              numCells        = 400,
                              plotAnimations  = FALSE,
                              AniPauseTime    = 0.5,
                              genome          = 'mm10',
                              chrom           = 'ALL',
                              excludeChrom    = NULL,
                              useOriStrength  = TRUE,
                              recycle         = TRUE,
                              plotMe          = FALSE,
                              simDSBs         = FALSE,
                              useSmoothing    = FALSE,
                              IODtime         = 50,
                              runToCompletion = FALSE,
                              noNormalize     = FALSE,
                              debug           = FALSE,
                              repStep         = 10,## Replication Step Size (Kb)
                              simOriSeq       = FALSE,
                              noSexChrom      = TRUE,
                              verbose         = FALSE,
                              okOut           = NULL){

  modelParams <- as.list(environment(), all=TRUE)
  
  if (is.numeric(recycle)){
    if (recycle == 1){recycle = TRUE}else{recycle=FALSE}
  }
  
  if (is.numeric(useOriStrength)){
    if (useOriStrength == 1){useOriStrength = TRUE}else{useOriStrength=FALSE}
  }
  
  if (useOriStrength){
    usingStrengthString = 'UsingStrength'
  }else{
    usingStrengthString = 'UniformStrength'
  }
  
  if (recycle){
    recyclingString = 'withRecycling'
  }else{
    recyclingString = 'noRecycling'
  }
  
  if (is.null(okOut)){
    okazakiOutFile <- paste0('okazaki_',genRandomString(len=15),'.tab')
  }else{
    okazakiOutFile <- okOut;
  }
  
  ##KBIODtime <- IODtime/10
  IODtime       <- IODtime/repStep

  ## SET Sizes
  repStepBP     <- repStep*1000
  repStepBPhalf <- repStep*1000/2
  
  ## Get chromosome sizes
  allCSsizes     <- getChromSize(genome,'ALL',
                                 noSexCS = noSexChrom,
                                 noM=TRUE)
    
  chromSize      <- max(allCSsizes$max)
  
  ## Import origins (with strength)
  oriCS         <- getALLOrigins(originsFile,allCSsizes,repStepBP)
  
  ## To exclude chromosomes, we will simply remove all origins from those chromosomes
  if (!is.null(excludeChrom)){
    for (eC in strsplit(excludeChrom,split = "chr")[[1]]){
      if (eC != ""){
        exclCS <- paste0("chr",eC)
        oriCS <- oriCS[oriCS$cs != exclCS,]
        print(paste0('Excluding ',exclCS,' ...'))
      }
    }
  }
  
  modelParams$csLength <- max(oriCS$max)
  
  ## Import DSBs (with strength)
  if (!is.null(dsbFile)){
    dsbCS         <- getALLDSBs(genome,allCSsizes,repStepBP,myDSBs = dsbFile)
  }
  
  numOri            <- min(round(oriPerMb*chromSize/1000000),length(oriCS$cs))

  #KB XX repLength    <- (repRate * timePassed)/10
  repLength           <- (repRate * timePassed)/repStep
  
  print(paste0(repLength,' cycles'))
  
  totalForks          <- numOri*2;
  ## KB 181115: This is silly
  ## Add noise later, ONLY when 
  #percentReplicating  <- pcReplicating
  #nReplicatingCells   <- round(numCells*pcReplicating/100)
  #nOtherCells         <- numCells - nReplicatingCells
  nReplicatingCells <- numCells
  ## END
  
  ## Make cs Nkb windows
  zStart              <- seq((repStepBPhalf+1)          ,chromSize,repStepBP)
  zEnd                <- seq((repStepBP+repStepBPhalf+1),chromSize,repStepBP)
  
  modelDF               <- data.frame("from"        = zStart[1:length(zEnd)],
                                      "end"         = zEnd,
                                      "ori"         = 0,
                                      "strength"    = 0,
                                      "rawstrength" = NA,
                                      "repModel"    = 2,
                                      "sim"         = FALSE,
                                      "boundaries"  = FALSE,
                                      "type"        = "Replication")
  
  modelDF$type <- factor(modelDF$type,levels=c("Replication","Origin firing"))
  
  ## Mark each chromosome in model DF
  modelDF$whatCS    <- 'unknown'
  modelDF$trueCoord <- modelDF$from
  
  for (nCS in 1:length(allCSsizes$cs)){
    myCO                    <- modelDF$from >= allCSsizes$min[nCS] & 
                               modelDF$from < allCSsizes$max[nCS]
    modelDF$whatCS[myCO]    <- as.character(allCSsizes$cs[nCS])
    modelDF$trueCoord[myCO] <- modelDF$from[myCO] - allCSsizes$min[nCS] + 1
  }
  
  allChromosomes <- levels(allCSsizes$cs)
  modelDF$whatCS <- factor(modelDF$whatCS,levels=allChromosomes)
  
  ## Calculate origin strength per window
  # Cannot use match because of cases where there is >1 origin per interval
  for (i in 1:dim(oriCS)[1]){
    nPos              <- modelDF$from == oriCS$rFrom[i]
    modelDF$ori[nPos] <- modelDF$ori[nPos] + oriCS$strength[i] 
  }

  # 2. get random Origins
  ## KB 181115: Do this later ... not here
  #modelDF$repModel    <- 2*nOtherCells
  modelDF$repModel    <- 0
  #END
  
  modelDF$dsbModel    <- 0
  modelDF$iod         <- 0
  modelDF$dsb         <- 0
  modelDF$dsbStrength <- 0
  modelDF$affySeq     <- 0
  
  if (!is.null(dsbFile)){
    for (i in 1:dim(dsbCS)[1]){
      nPos                      <- modelDF$from == dsbCS$rFrom[i]
      modelDF$dsbStrength[nPos] <- modelDF$dsbStrength[nPos] + dsbCS$strength[i] 
      modelDF$affySeq[nPos]     <- modelDF$affySeq[nPos] + dsbCS$affySeq[i] 
      modelDF$dsb[nPos]         <- 1 
    }
  }
  
  modelDF$affySeq <- modelDF$affySeq/max(modelDF$affySeq)
  
  fibers <- data.frame("pos"=modelDF$from)
  
  print(nReplicatingCells)
  print(repLength)
  
  oriTime    <- vector(length=nReplicatingCells*repLength*2) * 0 
  oriCell    <- vector(length=nReplicatingCells*repLength*2) * 0 
  oriCount   <- vector(length=nReplicatingCells*repLength*2) * 0
  oriForkDT  <- vector(length=nReplicatingCells*repLength*2) * 0 
  oriRepDone <- vector(length=nReplicatingCells*repLength*2) * 0 
  nLoop <- 1
  
  mReplicationPerMin   <- matrix(nrow = length(modelDF$from), ncol=repLength, data = 0)
  mReplicationPerCell  <- matrix(nrow = length(modelDF$from), ncol=nReplicatingCells, data = 0)
  mReplicationTiming   <- matrix(nrow = length(modelDF$from), data = 0)
  maxRT                <- 0

  #  ## KB 181125: no need for nOther HERE
  # KB 181113: 
  #  simByTime <- array(data = 2*nOtherCells, 
  #                     dim = c(length(modelDF$from), repLength*1.4),
  #                     dimnames = list(modelDF$from,1:(repLength*1.4)))
  
  repLenForArray <- round(repLength*1.4)
  simByTime <- array(data = 0, 
                     dim = c(length(modelDF$from), repLenForArray),
                     dimnames = list(modelDF$from,1:(repLenForArray)))
  
  #END
  
  cellByTime <- simByTime
  ## END
  
  mActiveForks <- array(data = 0, 
                        dim = c(nReplicatingCells, round(repLength*1.4,0), length(allChromosomes)),
                        dimnames = list(1:nReplicatingCells, 1:round(repLength*1.4,0), allChromosomes))
  
  ## Set interchromosomal boundaries so forks cannot elongate from one cs to next
  modelDF$boundaries[modelDF$from %in% aggregate(modelDF$from,
                                                 by=list(modelDF$whatCS),
                                                 max)$x] <- TRUE
  
  forksByCell <- list()
  
  for (nReps in 1:nReplicatingCells){
    
    mLen <- length(modelDF$iod)
    
    mDataSim  <- rep(FALSE,mLen)
    mDataIOD  <- rep(0,mLen)
    mDataDSB  <- rep(0,mLen)
    #mDataDSB  <- modelDF$dsbStrength
    mDataCyc  <- rep(0,mLen)
    mDataFrom <- modelDF$from
    mDataTo   <- modelDF$end
    mDataAffy <- modelDF$affySeq
    mDataType <- rep('',mLen)
    
    nDSBs          <- 0
    
    dsbsForThisCS  <- round(250 * chromSize / 2700000000) 
    
    oriForThisCell <- oriCS
    
    randOri <- getRandomOrigins(inDF = oriCS,
                                numForks = numOri,
                                useStrength = useOriStrength)
    
    randOriginPositions <- sort(match(randOri,table = mDataFrom))
    oriForThisCell <- oriForThisCell[!(oriCS$rFrom %in% mDataFrom[randOriginPositions]),]
    
    ## OK, so now each fork is an independent entity 
    #  We work on forks from now on, not origins
    forksR   <- data.frame('ori'         = randOriginPositions, 
                           'pos'         = randOriginPositions, 
                           'dir'         = 1 , 
                           timeFired     = 0, 
                           distTravelled = 0,
                           'ok'          = TRUE)
    
    forksL   <- data.frame('ori'         = randOriginPositions, 
                           'pos'         = randOriginPositions, 
                           'dir'         = -1 , 
                           timeFired     = 0,
                           distTravelled = 0,
                           'ok'          = TRUE)
    
    allForks <- rbind(forksL,forksR)
    
    ## select replication distance from distribution
    reps <- round(rnorm(1,repLength,2))
    reps[reps<0] <- 0

    ## for simOriSeq
    randomTimeToAbortSim <- sample(1:reps,1)
    
    if (verbose){
      print(paste0("Cell ",nReps," of ",nReplicatingCells,
                  " ... [OKabort = ",randomTimeToAbortSim,
                  ' of ',reps,']'))
    }

    ## Generate coverage @ origin
    #modelDF$sim[allForks$pos]  <- TRUE
    #modelDF$iod[allForks$pos]  <- 1
    #modelDF$type[allForks$pos] <- 'Origin firing'
    
    mDataSim[allForks$pos]  <- TRUE
    mDataIOD[allForks$pos]  <- 1
    mDataType[allForks$pos] <- 'Origin firing'

    modelPlot               <- modelDF;
    mDataType               <- factor(mDataType,levels=c("Replication","Origin firing"))
    
    ## for each replicative cycle
    #  1. Replicate repStep Kb (10Kb) @ each fork
    #  2. Recycle "activators" - require 2 fork collapses to fire new origin
    for (x in 1:reps){
      ## Move forks by 10 Kb
      
      allForks$pos[allForks$ok]           <- allForks$pos[allForks$ok] + allForks$dir[allForks$ok]
      allForks$distTravelled[allForks$ok] <- allForks$distTravelled[allForks$ok] + 1
      
      ## Remove forks @ replicated DNA and at start / end of CS
      #allForks$ok <- !(modelDF$sim[allForks$pos] == 2) & allForks$ok
      #allForks$ok <- !(modelDF$sim[allForks$pos] == 2) & allForks$pos > 0 & allForks$pos < length(modelDF$type) & allForks$ok
      
      #terSites <- allForks$pos[allForks$ok & mDataSim[allForks$pos]]
      #print(paste0(length(terSites),' forks will stop this cycle [',x,']'))
                               
      allForks$ok <- !(mDataSim[allForks$pos])           & 
                     !(modelDF$boundaries[allForks$pos]) &
                       allForks$pos > 0                  & 
                       allForks$pos < length(mDataType)  & 
                       allForks$ok
      
      ## Replicate DNA
      mDataSim[allForks$pos[allForks$ok]]  <- TRUE
      
      if (x <= IODtime){
        if (x <= (IODtime/2)){
          mDataIOD[allForks$pos[allForks$ok]]  <- 2
        }else{
          mDataIOD[allForks$pos[allForks$ok]]  <- 3
        }
      }
      
      if (plotAnimations){
        #modelDF$type <- 'Replication'
        mDataType <- rep('Replication',mLen)
      }

      ## License DSB hotspots
      if ((dsbsForThisCS > nDSBs) & simDSBs){
        forksPos    <- allForks$pos[allForks$ok]
        #forksAtDSBs <- forksPos[mDataDSB[forksPos]>0]
        forksAtDSBs <- forksPos[modelDF$dsbStrength > 0]
        
        #pDSB        <- mDataDSB[forksAtDSBs]
        pDSB        <- modelDF$dsbStrength[forksAtDSBs]
        
        DSBs2Make   <- forksAtDSBs[pDSB >= sample(1:1000000,length(forksAtDSBs))/500000]
        nDSBsNow    <- length(DSBs2Make)
        
        mDataDSB[DSBs2Make]  <- 1

        nDSBs       <- sum(mDataDSB>0)
        print(sprintf('Cell %1.0f; cycle %1.0f ... %1.0f DSBs selected [%1.0f total - %1.0f limit]',nReps,x,nDSBsNow,nDSBs,dsbsForThisCS))
      }
      
      ## Recycle ori activators
      numCurrentForks <- sum(allForks$ok)
      
      numNewOris <- floor((totalForks-numCurrentForks) / 2)
      
      ## Create new forks if enough old ones have run to completion
      if (sum(allForks$pos > length(mDataSim)) > 0){
        print('ERROR !!!!!!')
        print(paste0('totForks = ',totalForks,' ; curForks = ',numCurrentForks,' ; naForks = ',sum(is.na(allForks$ok))))
        print(paste0('R = ',recycle,
                   ' ; nO = ',numNewOris,
                   ' ; x = ',x,
                   ' ; cell# = ',nReps,
                   ' ; badForks = ',sum(allForks$pos > length(mDataSim)),
                   ' ; maxFork = ',max(allForks$pos),
                   ' ; allForks = ',sum(mDataSim[allForks$pos]),
                   ' ; okForks = ',length(allForks$pos<length(modelDF$type))))
      }
      
      if (numNewOris > 0 & recycle){
        
        ## Ensure that machinery is not recycled to origins within replicated domains
        oriForThisCell <- oriForThisCell[!(oriForThisCell$rFrom %in% mDataFrom[mDataSim]),]

        if (numNewOris > length(oriForThisCell$cs)){
          numNewOris <- length(oriForThisCell$cs)
        }        
        
        if (numNewOris > 0){
          ## Pick new origins
          randOri <- getRandomOrigins(inDF = oriForThisCell,
                                      numForks = numNewOris,
                                      useStrength = useOriStrength)
          
          randOriNew <- sort(match(randOri,table = mDataFrom))
          
          ## Remove chosen origins from future consideration
          oriForThisCell <- oriForThisCell[!(oriForThisCell$rFrom %in% mDataFrom[randOriNew]),]
          
          mDataSim[randOriNew]   <- TRUE
          mDataIOD[randOriNew]   <- 1
          mDataType[randOriNew]  <- 'Origin firing'
          
          forksR   <- data.frame('ori'         = randOriNew, 
                                 'pos'         = randOriNew, 
                                 'dir'         = 1 , 
                                 timeFired     = x, 
                                 distTravelled = 0,
                                 'ok'          = TRUE)
          
          forksL   <- data.frame('ori'         = randOriNew, 
                                 'pos'         = randOriNew, 
                                 'dir'         = -1 , 
                                 timeFired     = x, 
                                 distTravelled = 0,
                                 'ok'          = TRUE)
          
          allForks <- rbind(allForks,forksL,forksR)
        }
      }
      
      #modelDF$cycle <- x
      mDataCycle <- rep(x,mLen)
      
      nrow <- length(mDataFrom)
      
      if (plotAnimations && nReps == nReplicatingCells){
        if (x == 1){
          modFrom  <- integer(nrow)
          modSim   <- integer(nrow)
          modType  <- character(nrow)
          modCycle <- integer(nrow)
        }
        
        nFrom <- (1+(nrow*(x-1)))
        nTo   <- nFrom + nrow - 1  
        modFrom[nFrom:nTo]  <- mDataFrom; 
        modSim[nFrom:nTo]   <- mDataSim*2; 
        modCycle[nFrom:nTo] <- mDataCycle; 
        modType[nFrom:nTo]  <- mDataType; 
      }
      
      if (x <= repLength){
        mReplicationPerMin[,x]  <- mReplicationPerMin[,x] + mDataSim*2
        #noTimingYet <- ((mReplicationPerCell[,nReps] == 0) & modelDF$sim)
        mReplicationPerCell[((mReplicationPerCell[,nReps] == 0) & mDataSim),nReps] <- x
        
        if (length(modelDF$whatCS[allForks$ok])>0){
          dfForks      <- data.frame(fork=modelDF$whatCS[allForks$ok],
                                     cycle=x,
                                     cell=nReps)
  
          forkCountsDF <- plyr:::join(data.frame(Var1=allChromosomes),
                                      as.data.frame(table(modelDF$whatCS[allForks$pos[allForks$ok]])),
                                      by='Var1',type='left')
          
          forkCountsDF$Freq[is.na(forkCountsDF$Freq)] <- 0
          #print(sum(allForks$ok))
          #print(dim(mActiveForks))
          #print(dim(forkCountsDF))
          #print(as.data.frame(table(modelDF$whatCS[allForks$ok])))
          #print(forkCountsDF$Freq)
          mActiveForks[nReps,x,] <- forkCountsDF$Freq
        }
      }
      
      nLoop              <- nLoop + 1
      oriCount[nLoop]    <- sum(allForks$ok)
      oriTime[nLoop]     <- x
      oriCell[nLoop]     <- nReps
      #oriForkDT[nLoop]  <- mean(allForks$distTravelled[allForks$ok])
      oriRepDone[nLoop]  <- sum(mDataSim) / length(mDataSim) * 100
      mReplicationTiming <- mReplicationTiming + mDataSim
      maxRT <- maxRT+1
      
      simByTime[,x]      <- simByTime[,x] + mDataSim*2 + 2
      cellByTime[,x]     <- cellByTime[,x] + 1
      
      if (simOriSeq & x == randomTimeToAbortSim){
        ## Get fork position + random distance within this interval (allows coarse sim)
        #forkPos <- allForks$pos[allForks$ok] * 1000 * repStep + sample(1000*repStep,sum(allForks$ok))
        forkPos <- modelDF$trueCoord[allForks$pos[allForks$ok]] + sample(1000*repStep,sum(allForks$ok))
        forkDir <- allForks$dir[allForks$ok]
        forkCS  <- modelDF$whatCS[allForks$pos[allForks$ok]]
        
        for (nOkazaki in 1:5){
          okazakiFrom              <- forkPos + (150*(nOkazaki-1)) 
          okazakiFrom[forkDir>0]   <- forkPos[forkDir>0]-(150*nOkazaki)
          okazakiTo                <- okazakiFrom+150
          okazakiStrand            <- forkPos
          okazakiStrand[forkDir<0] <- '+'
          okazakiStrand[forkDir>0] <- '-'
          
          if (nOkazaki==1){
            dOkazaki <- data.frame(chrom=forkCS,
                                   from=format(okazakiFrom, scientific = FALSE),
                                   to=format(okazakiTo, scientific = FALSE),
                                   loop=x,
                                   cell=nReps,
                                   strand=okazakiStrand)
          }else{
            dOkazaki <- rbind(dOkazaki,data.frame(chrom=forkCS,
                                                  from=format(okazakiFrom, scientific = FALSE),
                                                  to=format(okazakiTo, scientific = FALSE),
                                                  loop=x,
                                                  cell=nReps,
                                                  strand=okazakiStrand))
          }
        }

        write.table(x = dOkazaki, 
                    append = TRUE, 
                    file = okazakiOutFile,
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE
                    )
        
        if (verbose){
          print(paste0('Exiting at loop ',x))
        }
        break
      }
      
    }
    
    if (plotAnimations && nReps == nReplicatingCells){
      print(nReps)
      modelPlot     <- data.frame(from  = modFrom,
                                  sim   = modSim,
                                  type  = modType,
                                  cycle = modCycle)
    }
    
    modelPlot$type <- factor(modelPlot$type,levels = c("Replication","Origin firing"))
    
    fibers[,paste0('R',nReps)]   <- mDataIOD
    
    modelDF$sim      <- mDataSim
    modelDF$repModel <- modelDF$repModel + mDataSim*2 + 2
    modelDF$dsbModel <- modelDF$dsbModel + mDataDSB
    
    allForks$cell <- nReps
    
    if (nReps == 1){
      returnForks <- allForks
    }else{
      returnForks <- rbind(returnForks,allForks)
    }
    forksByCell[[paste0('cell',nReps)]] <- allForks
  }
  
  # ### Output origins & ters
  # for (n in names(forksByCell)){
  #   forksByCell[[n]]$oriCS   <- modelDF$whatCS[forksByCell[[n]]$ori]
  #   forksByCell[[n]]$oriFrom <- modelDF$trueCoord[forksByCell[[n]]$ori]-1
  #   forksByCell[[n]]$oriTo   <- modelDF$trueCoord[forksByCell[[n]]$ori]
  #   
  #   forksByCell[[n]]$terCS   <- modelDF$whatCS[forksByCell[[n]]$pos]
  #   forksByCell[[n]]$terFrom <- modelDF$trueCoord[forksByCell[[n]]$pos]-1
  #   forksByCell[[n]]$terTo   <- modelDF$trueCoord[forksByCell[[n]]$pos]
  #   
  #   if (n == 'cell1'){
  #     allForkDetails <- forksByCell[[n]]
  #   }else{  
  #     allForkDetails <- rbind(allForkDetails,forksByCell[[n]])
  #   }
  # }
  
  #Duplicates ENV & uses waaaaay too much memory ... sort later ... 
  #if (simOriSeq){
  #  system(paste0('sort -k1,1 -k2n,2n -k3n,3n -o ',okazakiOutFile,' ',okazakiOutFile))
  #}
  
  if (plotAnimations){
    modelPlot$title <- factor(paste0(modelPlot$cycle*repStep," min (or Kb)"),
                              levels=paste0(unique(sort(modelPlot$cycle*10))," min (or Kb)"))
  }
  
  oriData    <- data.frame('nOri'=oriCount,
                           'time'=oriTime,
                           'cell'=oriCell,
                           "repPC"=oriRepDone)
  
  modelDF$rw <- smoothAndLog2Coverage(modelDF$repModel,applySmoothing = useSmoothing)
  
  simByTime <- simByTime[,apply(simByTime,2,sum)!=0]
  simByTime <- simByTime[,!is.na(apply(simByTime,2,sum))]
  
  modelDF$simByTime  <- simByTime
  modelDF$cellByTime <- cellByTime

  ## END
  
  modelParams$numSimCells    <- nReps
  modelParams$repIterations  <- repLength
  modelParams$repTime        <- timePassed
  
  dfForksPerCell <- as.data.frame.table(mActiveForks, 
                                        responseName = "nForks") 
  
  names(dfForksPerCell) <- c('cell','time','cs','nForks')
  
  aggForksPerCell1 <- aggregate(dfForksPerCell$nForks,
                               list(dfForksPerCell$time,
                                    dfForksPerCell$cs),
                               function(x){mean(x[x>0])}) ## Function ignores 0 values

  names(aggForksPerCell1) <- c('time','cs','nForks')
  
  aggForksPerCell <- join(aggForksPerCell1,allCSsizes,by='cs',type='inner')
  
  aggForksPerCell <- aggForksPerCell[!(is.na(aggForksPerCell$nForks)),]
  
  gForksPerCell   <- ggplot(aggForksPerCell,
                            aes(x=as.numeric(time)*10/60,
                                y=nForks/size*1000000/2,
                                group=cs)) + 
    geom_point(alpha=.3,size=.2,color='grey50') + 
    geom_smooth(se = FALSE, span=.1, lwd=.4, color='orange') +
    facet_wrap(~cs,ncol=4) + 
    geom_hline(yintercept=oriPerMb,lwd=.2) + 
    theme(legend.position='none',strip.background=element_blank(),
          strip.text=element_blank()) + 
    geom_text(aes(label=paste0(" ",cs)),
              x=-Inf,y=Inf,
              hjust=0,vjust=1,
              check_overlap = TRUE,
              size=7*5/14) + 
    xlab('Replication time (hr)') + 
    ylab('Oris per Mb') 
  
  gForksPerCellnoFacet   <- ggplot(aggForksPerCell,
                            aes(x=as.numeric(time)*10/60,
                                y=nForks/size*1000000/2,
                                group=cs)) + 
    geom_point(alpha=.3,size=.2,color='grey50') + 
    geom_smooth(se = FALSE, span=.1, lwd=.4, color='orange') +
    geom_hline(yintercept=oriPerMb,lwd=.2) + 
    theme(legend.position='none',strip.background=element_blank(),
          strip.text=element_blank()) + 
    xlab('Replication time (hr)') + 
    ylab('Oris per Mb')
    
  return(list(model=modelDF,
              params=modelParams,
              fibers=fibers,
              forks=allForks,
              forkDets=forksByCell,
              repPerMin=mReplicationPerMin,
              repPerCell=mReplicationPerCell,
              repTiming=mReplicationTiming,
              activeForks=list(df=dfForksPerCell,
                               agg=aggForksPerCell,
                               fig=gForksPerCell,
                               figNoFacet=gForksPerCellnoFacet),
              aniPlot=modelPlot))
}

#######################################################
#######################################################
#######################################################
#######################################################
simRT_addBGdata <- function(rtModel      = NULL,
                            RTbg         = '170223_0532_5_TGACCA_GC_021716_WGS_EarlyS_RT.oriSeq.mm10.log2.w500ks50k.WG.forModelling.bedgraph',
                            RTname       = NULL,
                            genome       = NULL,
                            isLog2       = TRUE,
                            noSexChrom   = TRUE,
                            useSmoothing = TRUE){

  csData <- getChromSize(myGenome = genome,
                         myCS     = 'ALL',
                         noM      = TRUE,
                         noSexCS  = noSexChrom)
  
  ## Get RT data
  fromRT <- getALLRT(inRT     = RTbg,
                     csDets   = csData,
                     myModel  = rtModel$model,
                     lSmooth  = useSmoothing,
                     lLog2    = isLog2,
                     inRTdata = NULL)
  ### FIX THIS !!!! StepSz MUST be passed !!
  
  rtModel[[RTname]]      <- fromRT
  
  ##NB ... Each time you add RT data, it becomes the default comparator for later !!!
  rtModel$experimentalRT <- fromRT
  
  return(rtModel)
}

#######################################################
#######################################################
#######################################################
#######################################################
simRT_compare2Data <- function(myModel      = NULL,
                               genome       = '',
                               chrom2use    = NULL,
                               RTdata       = NULL,
                               useSmoothing = TRUE,
                               isLog2       = TRUE,
                               noNormalize  = TRUE,
                               timepoint    = NULL,
                               pcRep        = NULL,
                               noSexChrom   = TRUE,
                               nMaxTime     = NULL,
                               normByCS     = FALSE,
                               plotMe       = TRUE){

  if (is.null(chrom2use)){
    modelDF              <- myModel$model
    if (!is.null(RTdata)){
      modelDF$strength     <- RTdata
    }
    
  }else{
    nCS                  <- myModel$model$whatCS == chrom2use
    myModel$params$chrom <- chrom2use
    myModel$model        <- myModel$model[nCS,]
    modelDF              <- myModel$model
    if (!is.null(RTdata)){
      modelDF$strength     <- RTdata[nCS]
    }
  }

  modelParams <- myModel$params
  
  ## Get time 2 finish replication
  nTime2Completion <- min(which(apply(myModel$model$simByTime < (myModel$params$numCells * 0.9 * 4),2,sum) / dim(myModel$model$simByTime)[1] <= 0.1))
  
  if (is.infinite(nTime2Completion)){
    nTime2Completion <- myModel$params$repTime*1.5
    print(paste0('## WARNING ## Model did not run to completion ... '))
  }
  
  print(paste0('Time to replicate (90% of genome in 90% of cells)= ',nTime2Completion,' cycles'))
  
  ## KB 181113: get stats for all timepoints
  if (is.null(nMaxTime)){
    maxTime  <- dim(modelDF$simByTime)[2]
  }else{
    maxTime  <- nMaxTime
  }
  
  noiseRange  <- seq(10,90,20)
  statsByTime <- array(data = 0, dim = c(maxTime,length(noiseRange),7)) ## KB 07/25/19: Add MASE
  allRW       <- array(data = 0, 
                       dim = c(length(modelDF$sim),
                               maxTime,
                               length(noiseRange)))
  
  if (normByCS){
    modelDF$strength <- log2(modelDF$strength/median(modelDF$strength))
  }  
  
  BGdata <- modelDF$strength
  
  ## KB : 07/25/19 : Get MAX/MIN RMSE & R2 values
  minRMSE <- 0
  minR2   <- 0  
  maxRMSE <- getRMSE(BGdata,rep(mean(BGdata),length(BGdata)))
  maxR2   <- 1
  
  ## Gen stats for ALL times and noise levels
  for (n in 1:maxTime){
    for (o in 1:length(noiseRange)){
      mStat <- testModelVData(myModel,BGdata,noiseRange[o],n,noPlot=TRUE,normBoth=TRUE)
      mStat$adjRMSE      <- round(1 - mStat$RMSE/maxRMSE,4)
      statsByTime[n,o,1] <- mStat$time
      statsByTime[n,o,2] <- mStat$pcRep
      statsByTime[n,o,3] <- mStat$R2
      statsByTime[n,o,4] <- round(mStat$RMSE,4)
      statsByTime[n,o,5] <- round(mStat$SSE,4)
      statsByTime[n,o,6] <- mStat$adjRMSE ## KB 07/25/19: Max adjusted RMSE
      statsByTime[n,o,7] <- mStat$R2 + 0.25*(mStat$adjRMSE)
      #statsByTime[n,o,7] <- mStat$R2 * (1+(round(1 - mStat$RMSE/maxRMSE,4)))
      
      print(paste0('N = ',mStat$time,
                   ' [',maxTime,'] ; ',
                   'o = ',mStat$pcRep,' ; ',
                   'add = ',mStat$numNoiseCells2Add, ' ; ',
                   'R2 = ',mStat$R2, ' ; ',
                   'score = ',statsByTime[n,o,7], ' ; ',
                   'RMSE = ',round(mStat$RMSE,4), ' ; ',
                   'adjRMSE = ',statsByTime[n,o,6]))
    }
  }
  
  ## Make DF of stats 
  dfStatsByTime <- as.data.frame(statsByTime)
  dfStatsByTime <- adply(statsByTime,
                         .margins = 1:2,
                         .fun = NULL)
  
  names(dfStatsByTime) <- c('V1','V2','time','pcRep','R2','RMSE','SSE','adjRMSE','score') ## KB 07/25/19: Add MASE
  
  dfStatsByTime <- dfStatsByTime[,3:9] ## KB 07/25/19: Add MASE
  
  #bestStat <- dfStatsByTime[min(which(dfStatsByTime$RMSE == min(dfStatsByTime$RMSE))),]
  #bestStat <- dfStatsByTime[min(which(dfStatsByTime$R2 == max(dfStatsByTime$R2))),]
  bestStat  <- dfStatsByTime[min(which(dfStatsByTime$score == max(dfStatsByTime$score,na.rm=TRUE))),]
  
  ## If none selected, use the BEST (for now)
  if (is.null(timepoint)){
    timepoint <- bestStat$time
  }
  
  if (is.null(pcRep)){
    pcRep     <- bestStat$pcRep
  }
  
  print(paste0('Using timepoint ',timepoint,' ...'))
  print(paste0('Using ',pcRep,'% replicating cells ...'))
  
  numNoiseCells2Add                  <- round(myModel$params$numCells/pcRep*100)
  
  simNoise                           <- rep(2*numNoiseCells2Add,length(modelDF$simByTime[,1])) 
  modelDF$rw2use                     <- smoothAndLog2Coverage(modelDF$simByTime[,timepoint] + simNoise,applySmoothing = TRUE)
  
  modelDF$strOK                      <- modelDF$strength
  modelDF$strength[modelDF$strOK==0] <- NA
  modelDF$rw2use[modelDF$strOK==0]   <- NA

  if (noNormalize){
    modelDF$strength                 <- as.numeric(modelDF$strength)
    modelDF$rw2use                   <- as.numeric(modelDF$rw2use)
  }else{
    modelDF$strength                 <- normalize1to1(as.numeric(modelDF$strength))
    modelDF$rw2use                   <- normalize1to1(as.numeric(modelDF$rw2use))
  }

  ok   <- abs(modelDF$strength) > 0.0001
  
  checkMe <- (modelDF$strength[ok] + modelDF$rw2use[ok])
  
  if (sum(!(is.nan(checkMe) | is.na(checkMe) | is.infinite(checkMe))) > 0){
    R2   <- round(cor(modelDF$strength[ok],modelDF$rw2use[ok],use='complete.obs'),2)
    RMSE <- round(getRMSE(modelDF$strength[ok],modelDF$rw2use[ok]),3)
    SSE  <- round(getSSE(modelDF$strength[ok],modelDF$rw2use[ok])/length(modelDF$rw2use[ok]),3)*100
  }else{
    R2   <- 0
    RMSE <- 1e99
    SSE  <- 1e99
  }

  modelDF$R2   <- 0
  modelDF$SSE  <- -999
  modelDF$RMSE <- -999
  
  for (cs in unique(modelDF$whatCS)){
    nP <- modelDF$whatCS == cs
    okCS <- ok & nP

    checkMe <- (modelDF$strength[okCS] + modelDF$rw2use[okCS])
      
    if (sum(!(is.nan(checkMe) | is.na(checkMe) | is.infinite(checkMe))) > 0){
      modelDF$R2[nP]   <- round(cor(modelDF$strength[okCS],modelDF$rw2use[okCS],use='complete.obs'),2)
      modelDF$RMSE[nP] <- round(getRMSE(modelDF$strength[okCS],modelDF$rw2use[okCS]),3)
      modelDF$SSE[nP]  <- round(getSSE(modelDF$strength[okCS],modelDF$rw2use[okCS])/length(modelDF$rw2use[okCS]),3)
    }else{
      modelDF$R2[nP]   <- 0
      modelDF$RMSE[nP] <- 1e99
      modelDF$SSE[nP]  <- 1e99
    }
  }
  
  theme7point()
  
  zStat <- reshape:::melt.data.frame(dfStatsByTime,
                                     id.vars = c('time','pcRep'))
  
  names(zStat) <- c('time','type','val')
  
  gStat <- ggplot(zStat,
                  aes(x=time*10,
                      y=val,
                      color=type)) + 
    geom_point(size=0.3) + 
    facet_grid(type~.,scales='free_y') + 
    theme(strip.background=element_blank(),
          strip.text=element_blank()) + 
    geom_text(x=-Inf,y=Inf,hjust=0,vjust=1,aes(label=paste0("  ",type)),
              check_overlap = TRUE,
              size=9*5/14) + 
    xlab('Time (min)') + 
    ylab('Value') + 
    theme(legend.position='none',
          panel.grid = element_line(size=.2,color='grey70'))
  
  gStatHM <- ggplot(dfStatsByTime,
                    aes(x=time,y=pcRep,fill=R2)) + 
    geom_tile(color='white') + 
    scale_fill_gradient2(low='black',
                         mid='red',
                         high='yellow',
                         midpoint=0.5)
  
  myModel$model <- modelDF
  
  sCS <- aggregate(modelDF[,c('R2','SSE','RMSE')],
                     by=list(modelDF$whatCS),
                     min)
  
  names(sCS)[1] <- 'cs'
  stats <- rbind(sCS,data.frame('cs'='ALL','R2'=R2,'SSE'=SSE,'RMSE'=RMSE))
  
  dfStatsByTime$oriPerMb    <- myModel$params$oriPerMb
  dfStatsByTime$useStrength <- myModel$params$useStrength
  dfStatsByTime$recycle     <- myModel$params$recycle
  
  myModel$params <- modelParams
  
  if (plotMe == TRUE){
    #mBest <- dfStatsByTime[dfStatsByTime$RMSE == min(dfStatsByTime$RMSE),][1,]
    mBest <- dfStatsByTime[dfStatsByTime$R2 == max(dfStatsByTime$R2),][1,]
    
    pc2use <- mBest$pcRep
    if (!is.null(pcRep)){
      #print (paste0('In 1 : ',pcRep))
      pc2use <- pcRep  
    }
    
    time2use <- mBest$time
    if (!is.null(timepoint)){
      #print (paste0('In 2 : ',timepoint))
      time2use <- timepoint 
    }
    
    print (paste0(pc2use,'% replicating ... for ',time2use,' cycles ... '))
    bestMod <- testModelVData(myModel,
                              BGdata,
                              pc2use,
                              time2use)
    
    g <- bestMod$fig
  
    myName <- bestMod$name  
    
    if (is.null(chrom2use)){
      nW <- 7
      nH <- 11
    }else{
      nW <- 6
      nH <- 3
    }
    
    ggsave(filename = getIMGname(myName,'PNG',saveLocation = ''),plot = g, width=nW, height=nH, units='in')
    ggsave(filename = getIMGname(myName,'PDF',saveLocation = ''),plot = g, width=nW, height=nH)
  }else{
    g <- ggplot()
  }
  
  dfStatsByTime$repTime <- nTime2Completion
  
  return(list("Score"=-1*RMSE,
              "Pred"=0,
              "pearson"=R2,
              "RMSE"=RMSE,
              "SSE"=SSE,
              'fig'=g,
              'stats'=stats,
              'bestStat'=bestStat,
              'model'=modelDF,
              'allRW'=allRW,
              "byTime" = list(raw=dfStatsByTime,
                              melt=zStat,
                              fig1=gStat,
                              fig2=gStatHM),
              'data'=myModel))
}

#######################################################
#######################################################
#######################################################
#######################################################
simRT_drawFit <- function(myModel      = NULL,
                          chrom2use    = NULL,
                          timepoint    = NULL,
                          pcRep        = NULL,
                          myName       = NULL,
                          noNormalize  = FALSE,
                          plotMe       = TRUE,
                          imgPath      = NULL){
  
  theme7point()
  
  if (is.null(myName)){
    myName <- paste0('RTmodel_',timepoint,'cycles_',pcRep,'pcRep')
  }
  
  if (is.null(chrom2use)){
    myModel$csPos                 <- 1:length(myModel$model$whatCS)
    modelDF                       <- myModel$model
    myName                        <- paste0(myName,'_ALLCS')
  }else{
    nCS                           <- myModel$model$whatCS == chrom2use
    myModel$csPos                 <- myModel$model$whatCS == chrom2use
    myModel$params$chrom          <- chrom2use
    myModel$model                 <- myModel$model[nCS,]
    modelDF                       <- myModel$model
    myName                        <- paste0(myName,'_',chrom2use)
    myModel$experimentalRT$smooth <- myModel$experimentalRT$smooth[nCS]
    myModel$experimentalRT$raw    <- myModel$experimentalRT$raw[nCS]
  }
  
  modelParams   <- myModel$params

  BGdata        <- myModel$experimentalRT$smooth
  
  mStat         <- testModelVData(myModel,
                                  BGdata,
                                  pcRep,
                                  timepoint,
                                  noPlot=FALSE,
                                  normBoth=!(noNormalize))
  
  #mStat$adjRMSE <- round(1 - mStat$RMSE/maxRMSE,4)
  #mStat$score   <- round(mStat$R2 + 0.25*(mStat$adjRMSE),4)

  mStat$fig
  
  print(paste0('N = ',mStat$time, ' ; ',
               'o = ',mStat$pcRep,' ; ',
               'add = ',mStat$numNoiseCells2Add, ' ; ',
               'R2 = ',mStat$R2))
  
  if (is.null(chrom2use)){
    nW <- 7; nH <- 6
  }else{
    nW <- 7; nH <- 2 
  }
  
  if(is.null(imgPath)){
    if (exists('imgOutputFolder')){
      imgPath <- imgOutputFolder
    }else{
      imgPath <- '.'
    }
  }
  
  if (!grepl('/$',imgPath)){
    imgPath <- paste0(imgPath,'/')
  }
  
  if (plotMe){
    print (paste0("Figs output to : ",getIMGname(myName,'PNG',saveLocation=imgPath)))
    
    ggsave(filename = getIMGname(myName,'PNG',saveLocation=imgPath),
           plot = mStat$fig,
           width=nW, height=nH, dpi=1000, units='in')
  
    ggsave(filename = getIMGname(myName,'PDF',saveLocation=imgPath),
           plot = mStat$fig,
           width=nW, height=nH)
  }
  
  return(list('fig'=mStat$fig,
              'stats'=mStat))
}

#######################################################
#######################################################
#######################################################
#######################################################
simRT_drawSingleCell <- function(myModel      = NULL,
                                 modelFile    = NULL,
                                 chroms2use   = 'chr11',
                                 myName       = NULL,
                                 outBG        = TRUE,
                                 ncells       = 100,
                                 imgPath      = NULL){
  
  if (!(is.null(modelFile))){
    load(verbose = TRUE,file = modelFile)
    myModel <- modelData$mod
  }
  
  theme7point()
  
  if (is.null(myName)){
    myName <- 'unnamed'
  }
  
  ## OUTPUT RTSim BEDGRAPH
  if (outBG){
    simRT_exportRTBG(myModel = modelData$mod, 
                     myName = myName,
                     type = 'sim')
    simRT_exportRTBG(myModel = modelData$mod, 
                     myName = myName,
                     type = 'exp')
  }
  
  retlist <- list()
  
  for (chrom2use in chroms2use){
    imgName <- paste0('RTmodel_SingleCellSimPlot_',chrom2use,'_',myName)
    
    nCS <- which(myModel$model$whatCS ==chrom2use)
    
    myModel$repPerCell[myModel$repPerCell==0] <- myModel$params$repIterations+5
    
    df <- as.data.frame(myModel$repPerCell[nCS,1:ncells])
    
    df$pos <- myModel$model$trueCoord[nCS]
  
    mDF <- reshape2:::melt.data.frame(df,id.vars=c('pos'))
    
    mDF$cell <- as.numeric(gsub('V','',mDF$variable))
    
    gSS <- ggplot(mDF,aes(x=pos/1000000,y=cell,fill=log(value))) +
      geom_tile() +
      coord_cartesian(expand = FALSE,
                      xlim=c(0,max(df$pos/1000000))) +
      scale_fill_gradient(high='dodgerblue4',low='yellow',
                          na.value = alpha('black',0.3)) + 
      xlab('') +
      theme(legend.position = 'none',
            panel.background = element_rect(fill='dodgerblue4'),
            axis.text.x=element_blank(),
            plot.margin = unit(c(-0.3,0,-0.3,0),'cm'))
    #xlab(paste0('Position on chromosome ',chrom2use,' (Mb)')) + 
  
    mnRep <- apply(abs(myModel$repPerCell[nCS,]),1,mean)
    sdRep <- apply(abs(myModel$repPerCell[nCS,]),1,sd)
    mnsd <- data.frame(pos=df$pos,mn=mnRep,mnMinus=mnRep-sdRep,mnPlus=mnRep+sdRep)
    
    mnsd$mnNormMNSD <- standardizeMNSD(mnsd$mn)
    mnsd$mnNorm01 <- standardize0to1(mnsd$mn)
    
    gMnSD <- ggplot(mnsd,aes(x=pos/1000000,y=mn)) + 
      geom_ribbon(aes(ymin=mnMinus,ymax=mnPlus),fill=alpha('grey40',.5)) + 
      geom_line(color='red',lwd=.3) + 
      coord_cartesian(expand = FALSE,xlim=c(0,max(df$pos/1000000)),
                      ylim=c(0,myModel$params$timePassed/10)) + 
      xlab(paste0('Position on chromosome ',chrom2use,' (Mb)')) + 
      ylab('RT (cycle)') +
      scale_y_continuous(breaks = seq(0,floor(myModel$params$repIterations/100)*100,
                                      length.out = 3)) +
      theme(legend.position = 'none',
            plot.margin = unit(c(0,0,0,0),'cm')) 
      
    
    gX2 <- ggarrange(gSS,gMnSD,align='v',heights=c(3,2),ncol=1,nrow=2)
      
    nW <- 3; nH <- 2 
  
    if(is.null(imgPath)){
      if (exists('imgOutputFolder')){
        imgPath <- imgOutputFolder
      }else{
        imgPath <- '.'
      }
    }
    
    print(getIMGname(imgName,'PNG',saveLocation=imgPath))
    
    ggsave(filename = getIMGname(imgName,'PNG',saveLocation=imgPath),
           plot = gX2,
           width=nW, height=nH, dpi=1000, units='in')
    
    ggsave(filename = getIMGname(imgName,'PDF',saveLocation=imgPath),
           plot = gX2,
           width=nW, height=nH)
    
    retlist[[chrom2use]]$figX2         <- gX2
    retlist[[chrom2use]]$figSingleCell <- gSS
    retlist[[chrom2use]]$figMNSD       <- gMnSD
    retlist[[chrom2use]]$dataMNSD      <- mnsd
    
  }
  
  return(retlist)
}

#######################################################
#######################################################
#######################################################
#######################################################
simRT_drawComparisonFigure <- function(mod       = NULL,
                                       modFile   = NULL,
                                       chroms2use = c('chr11','chr13'),
                                       myName     = NULL,
                                       ncells     = 100,
                                       imgPath    = NULL){
  
  if (!(is.null(modFile))){
    load(verbose = TRUE,file = modFile)
    myModel <- modelData$mod
    if (!exists('myModel$params$bestTime')){
      myModel$params$bestTime <- modelData$figData$bestStat$time
      myModel$params$bestPC   <- modelData$figData$bestStat$pcRep
    }
  }
  
  theme7point()
  
  if (is.null(myName)){
    myName <- 'unnamed'
  }
  
  retFigz <- list()
  
  for (chrom2use in chroms2use){
    
    z  <- simRT_drawFit       (myModel,
                               chrom2use,
                               myModel$params$bestTime, 
                               myModel$params$bestPC,
                               'KBTST1')
    
    z2 <- simRT_drawSingleCell(myModel,
                               chroms2use = chrom2use,
                               myName = 'KBTST2', 
                               ncells = ncells,
                               outBG = FALSE)
    
    noX <- theme(axis.text.x = element_blank())
    gAgg <- ggarrange(z$fig + coord_cartesian(expand = FALSE) + noX + xlab(''),
                      z2[[chrom2use]]$figSingleCell,
                      z2[[chrom2use]]$figMNSD + scale_y_reverse(),
              ncol=1,
              nrow=3,
              heights=c(2,3,2),
              align='v')
    
    retFigz[[chrom2use]]$figX3 <- gAgg 
    retFigz[[chrom2use]]$figA  <- z$fig + coord_cartesian(expand = FALSE) + noX + xlab('')
    retFigz[[chrom2use]]$figB  <- z2[[chrom2use]]$figSingleCell
    retFigz[[chrom2use]]$figC  <- z2[[chrom2use]]$figMNSD + scale_y_reverse()
  }

  return(retFigz)
}

#######################################################
#######################################################
#######################################################
#######################################################
simRT_exportRTBG    <- function(myModel      = NULL,
                                modelFile    = NULL,
                                myName       = NULL,
                                type         = 'sim',
                                imgPath      = NULL){
  
  if (!(is.null(modelFile))){
    load(verbose = TRUE,file = modelFile)
    myModel <- modelData$mod
  }
  
  theme7point()
  
  if (is.null(myName)){
    myName <- 'unnamed'
  }
  
  retlist <- list()
  
  myBGName <- paste0('./',myName,'.',type,'RT.WG.bedgraph')
  
  print (paste0("Output to ",myBGName," ..."))
  halfSz  <- round((myModel$model$trueCoord[2]-myModel$model$trueCoord[1])/2)
  
  ## Now exports either sim (default) or real RT
  if (type == 'sim'){
    myModel$repPerCell[myModel$repPerCell==0] <- max(myModel$repPerCell)+1
    
    mnRep <- apply(abs(myModel$repPerCell),1,mean)
    sdRep <- apply(abs(myModel$repPerCell),1,sd)
    
    dfForBG <- data.frame(cs=myModel$model$whatCS,
                          from=sprintf("%1i",myModel$model$trueCoord-halfSz),
                          to=sprintf("%1i",myModel$model$trueCoord+halfSz-1),
                          mn=(mnRep+1)*-1)
  }else{
    dfForBG <- data.frame(cs=myModel$model$whatCS,
                          from=sprintf("%1i",myModel$model$trueCoord-halfSz),
                          to=sprintf("%1i",myModel$model$trueCoord+halfSz-1),
                          mn=myModel$experimentalRT$smooth)
  }

  fwrite(dfForBG,
         file = myBGName,
         quote = FALSE, 
         sep = "\t",
         na = "0", 
         row.names = FALSE,
         col.names = FALSE
         )
}

#######################################################
# The next three funcitons are helper functions for ###
# simRT_drawGIF########################################
#######################################################
plotReplicationForAnimation <- function(dRep,rtLim,
                                        pScale=1,pSz=1,pLwd=1,
                                        pTxt=16){
  cs <- dRep$cs[1]
  
  dRep$RT[dRep$RT == 0] <- max(dRep$RT) + 1
  
  dRep$pos    <- dRep$pos/1000000
  dRep$RTok   <- (dRep$RT<rtLim)
  dRep$RTthis <- (dRep$RT==rtLim)
  csMax       <- max(dRep$pos+1)
  csMin       <- min(dRep$pos-1)
  
  dfOri   <- dRep[dRep$ori & dRep$RTok,]
  dfFire  <- dRep[dRep$ori & dRep$RTthis,]
  
  ## This is necessary to assure that dfRepli always has at least one entry.
  ## Otherwise dies once last replisome falls off ... 
  dfRepli           <- dRep[1,]
  dfRepli$pos[1]    <- csMin-5;
  dfRepli$RTthis[1] <- TRUE;
  dfRepli           <- rbind(dfRepli,dRep[dRep$RTthis,])
  
  pSz         <- pSz*pScale
  pLwd        <- pLwd*pScale
  textSize    <- pTxt*pScale
  
  g <- ggplot(dRep,aes(x=pos)) + 
    geom_line(aes(y=RTok*0.5),lwd=pLwd) + 
    geom_line(aes(y=RTok*-0.5),lwd=pLwd) +
    geom_point(data=dfRepli,aes(y=0),shape=21,fill='red',color='red',size=pSz*1.8) +
    geom_point(data=dfOri,
               y=-0.5,size=pSz*.4,
               color='green') + 
    geom_point(data=dfOri,
               y=0.5,size=pSz*0.4,
               color='green') + 
    geom_vline(data=dfFire,
               aes(xintercept=(pos)),
               color='forestgreen',lwd=0.6,alpha=.5) + 
    xlab(paste0('Position on ',cs,' (Mb)')) + 
    ylab('') + 
    scale_alpha_manual(values=c(0,1)) +
    theme(axis.text.x=element_text(size=textSize*.8),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y=element_blank(),
          legend.position='none',
          title = element_text(size=textSize),
          panel.grid = element_blank(),
          panel.background = element_rect(fill='grey97')) +
    coord_cartesian(ylim=c(-0.8,0.8),
                    xlim=c(csMin,csMax),
                    expand=FALSE)
  return(g)
}

# WHOLE GENOME ########################################
plotReplicationForAnimationWG <- function(dRep,rtLim,
                                          pScale=1,pSz=1,pLwd=1,
                                          pTxt=6){
  cs <- dRep$cs[1]
  
  dRep$csNum <- gsub("chr","",dRep$cs)
  
  dRep$RT[dRep$RT == 0] <- max(dRep$RT) + 1
  
  dRep$pos    <- dRep$pos/1000000
  dRep$RTok   <- (dRep$RT<rtLim)
  dRep$RTthis <- (dRep$RT==rtLim)
  csMax       <- max(dRep$pos+1)
  csMin       <- min(dRep$pos-1)
  
  dfOri   <- dRep[dRep$ori & dRep$RTok,]
  dfFire  <- dRep[dRep$ori & dRep$RTthis,]
  
  ## This is necessary to assure that dfRepli always has at least one entry.
  ## Otherwise dies once last replisome falls off ... 
  dfRepli           <- dRep[1,]
  dfRepli$pos[1]    <- csMin-5;
  dfRepli$RTthis[1] <- TRUE;
  dfRepli           <- rbind(dfRepli,dRep[dRep$RTthis,])
  #dfRepli <- dRep[dRep$RTthis,] 
  
  pSz         <- pSz*pScale
  pLwd        <- pLwd*pScale
  textSize    <- pTxt*pScale
  
  g <- ggplot(dRep,aes(x=pos)) + 
    geom_line(aes(y=RTok*0.5),lwd=pLwd) + 
    geom_line(aes(y=RTok*-0.5),lwd=pLwd) +
    geom_point(data=dfRepli,aes(y=0),shape=21,fill='red',color='red',size=pSz*1.8) +
    geom_point(data=dfOri,
               y=-0.5,size=pSz*.4,
               color='green') + 
    geom_point(data=dfOri,
               y=0.5,size=pSz*0.4,
               color='green') + 
    geom_vline(data=dfFire,
               aes(xintercept=(pos)),
               color='forestgreen',lwd=0.6,alpha=.5) + 
    xlab(paste0('Position (Mb)')) + 
    ylab('') + 
    scale_alpha_manual(values=c(0,1)) +
    theme(axis.text.x=element_text(size=textSize*.8),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y=element_blank(),
          legend.position='none',
          strip.placement = "outside",
          strip.text = element_text(size=textSize*0.8),
          title = element_text(size=textSize),
          panel.grid = element_blank(),
          panel.background = element_rect(fill='grey97')) +
    coord_cartesian(ylim=c(-0.8,0.8),
                    xlim=c(csMin,csMax),
                    expand=FALSE) +
    facet_wrap(~cs,ncol=2,
               strip.position = 'left')
  
  # geom_text(aes(label=csNum),
  #           x=-1,y=0,
  #           size=(textSize*0.8*5/14),
  #           check_overlap = TRUE,
  #           hjust=0.5,vjust=1) +
  
  return(g)
}

#######################################################
# The next three funcitons are helper functions for ###
# simRT_drawGIF########################################
#######################################################
plotRTForAnimation <- function(dMod,
                               rtLim,
                               cs=NULL,
                               pScale=1,
                               pSz=1,
                               pLwd=1,
                               pTxt=16){
  
  modRT <- data.frame(pos=dMod$mod$model$trueCoord/1000000,
                      cs=dMod$mod$model$whatCS,
                      RT=dMod$mod$repPerMin[,rtLim]/(dMod$mod$params$numCells*2)*100)
  
  modRT <- modRT[modRT$cs == cs,]
  
  ## Stupid end thing
  modRT <- modRT[modRT$pos < max(modRT$pos),]
  
  csMax       <- max(modRT$pos+1)
  csMin       <- min(modRT$pos-1)
  
  pSz         <- pSz*pScale
  pLwd        <- pLwd*pScale
  textSize    <- pTxt*pScale
  
  g <- ggplot(data=modRT) + 
    geom_line(aes(x=pos,y=RT),lwd=pLwd,color='darkred') + 
    xlab(paste0('Position (Mb)')) + ylab('') +
    scale_y_continuous(breaks=c(20,80),labels=c('20%','80%'),position='left')+
    theme(axis.text.x=element_text(size=textSize*.8),
          axis.text.y=element_text(size=textSize*.8),
          legend.position='none',
          title = element_text(size=textSize),
          panel.grid = element_blank(),
          panel.background = element_rect(fill='white'),
          strip.placement = "outside",
          strip.text = element_text(size=textSize*0.8)) +
    coord_cartesian(xlim=c(csMin,csMax),
                    ylim=c(0,105),
                    expand=FALSE)
  
  return(g)
}

# WHOLE GENOME ########################################
plotRTForAnimationWG <- function(dMod,
                                 rtLim,
                                 pScale=1,
                                 pSz=1,
                                 pLwd=1,
                                 pTxt=16){
    
    modRT <- data.frame(pos=dMod$mod$model$trueCoord/1000000,
                        cs=dMod$mod$model$whatCS,
                        RT=dMod$mod$repPerMin[,rtLim]/(dMod$mod$params$numCells*2)*100)
    
    csMax       <- max(modRT$pos+1)
    csMin       <- min(modRT$pos-1)
    
    pSz         <- pSz*pScale
    pLwd        <- pLwd*pScale
    textSize    <- pTxt*pScale
    
    g <- ggplot(data=modRT) + 
      geom_line(aes(x=pos,y=RT),lwd=pLwd,color='darkred') + 
      facet_wrap(~cs,ncol=2,strip.position = 'left') +
      xlab(paste0('Position (Mb)')) + 
      ylab('Rep (%)') + 
      scale_y_continuous(breaks=c(20,80),position='right')+
      theme(axis.text.x=element_text(size=textSize*.8),
            legend.position='none',
            title = element_text(size=textSize),
            panel.grid = element_blank(),
            panel.background = element_rect(fill='white'),
            strip.placement = "outside",
            strip.text = element_text(size=textSize*0.8)) +
      coord_cartesian(xlim=c(csMin,csMax),
                      ylim=c(0,105),
                      expand=FALSE)
    
    return(g)
  }

#######################################################
#######################################################
getModelDataForGIF <- function(chrom     = 'chr19',
                               from      = NULL,
                               to        = NULL,
                               cellNum   = 1,
                               modelPath = NULL,
                               model     = NULL){
  
  ## Get RTSim object
  if (is.null(model)){
    load('/data/RDCO/replicationPaper/pipe/accessoryFiles/bestRTsimModels/MeiS_ALL_S3_2to4C_STRA8_DMRT1n.Rdata')
  }else{
    modelData <- model
  }
  
  ## Make initial DF
  df <- data.frame(cs  = modelData$mod$model$whatCS,
                   pos = modelData$mod$model$trueCoord,
                   RT  = modelData$mod$repPerCell[,cellNum])
  
  
  ## Get region of interest
  if (is.null(from) | is.null(to)){
    dfSelect <- df[df$cs == chrom,]  
  }else{
    dfSelect <- df[df$pos > from & 
                   df$pos < to & 
                   df$cs == chrom, ]
    print(dim(dfSelect))
  }
  
  ## Melt DF
  dfMelted <- reshape2:::melt.data.frame(dfSelect,
                                         id.vars =c('cs','RT'),
                                         measure.vars = c('pos'))
  
  names(dfMelted) <- c('cs','RT','var','pos')
  
  ## Get origins from RT local minima
  dOri <- dfSelect[getLocalMinima(dfSelect$RT),]
  names(dOri) <- c('cs','pos','ori')
  
  ## Add oris to bigger DF
  dOri$ori                        <- TRUE
  dfFinal                         <- plyr::join(dfMelted,dOri,by=c('cs','pos'))
  dfFinal$ori[is.na(dfFinal$ori)] <- FALSE
  
  return(dfFinal) 
}



# WHOLE GENOME ########################################
getModelDataForGIFWG <- function(chrom     = 'chr19',
                                 from      = NULL,
                                 to        = NULL,
                                 cellNum   = 1,
                                 modelPath = '/data/RDCO/replicationPaper/pipe/accessoryFiles/bestRTsimModels/MeiS_ALL_S3_2to4C_STRA8_DMRT1n.Rdata',
                                 model     = NULL){
  
  ## Get RTSim object
  if (is.null(model)){
    load('/data/RDCO/replicationPaper/pipe/accessoryFiles/bestRTsimModels/MeiS_ALL_S3_2to4C_STRA8_DMRT1n.Rdata')
  }else{
    modelData <- model
  }
  
  ## Make initial DF
  df <- data.frame(cs  = modelData$mod$model$whatCS,
                   pos = modelData$mod$model$trueCoord,
                   RT  = modelData$mod$repPerCell[,cellNum])
  
  
  ## Get region of interest
  dfSelect <- df
  
  ## Melt DF
  dfMelted <- reshape2:::melt.data.frame(dfSelect,
                                         id.vars =c('cs','RT'),
                                         measure.vars = c('pos'))
  
  names(dfMelted) <- c('cs','RT','var','pos')
  
  ## Get origins from RT local minima
  dOri <- dfSelect[getLocalMinima(dfSelect$RT),]
  names(dOri) <- c('cs','pos','ori')
  
  ## Add oris to bigger DF
  dOri$ori                        <- TRUE
  dfFinal                         <- plyr::join(dfMelted,dOri,by=c('cs','pos'))
  dfFinal$ori[is.na(dfFinal$ori)] <- FALSE
  
  return(dfFinal) 
}

#######################################################
#######################################################
makeReplicationGIF <-function(dfWithOri,name='repMovie',cyc=550,nstep=5){
  
  plotReplicationForAnimation(dfWithOri,1,pScale=2,pSz=1,pLwd=.7,pTxt=12)
  
  ## chr19 will be 7 wide 
  ## calibrate others to mm10 chr19 (61.42119 Mb)... 
  relativeWidth <- (max(dfWithOri$pos)/1000000+1)/62.42119
  
  # Set animation options
  # Careful with interval and nmax combo !!
  ani.options(interval = 0.2,
              nmax=100,
              loop = 1,
              autobrowse=FALSE,
              ani.width=1200*relativeWidth,
              ani.height=130)
  
  saveGIF({
    rtEnd = cyc
    
    #Loop through time cycles
    for(i in seq(1,rtEnd,by=nstep)){
      
      p <- plotReplicationForAnimation(dfWithOri,i,pScale=2,pSz=1,pLwd=.7,pTxt=12)
      print(p)
      
    }#close the for loop
    
  }, convert = 'convert',
  img.name=paste0('allImg_',name),
  movie.name = paste0(name,'.gif')) 
  #close the animation builder
}

# WHOLE GENOME ########################################
makeReplicationGIFWG <-function(dfWithOri,name='repMovie',cyc=550,nstep=5){
  #plotReplication(dfWithOri,1,pScale=2,pSz=1,pLwd=.7,pTxt=12)
  plotReplicationForAnimationWG(dfWithOri,1,pScale=2,pSz=.6,pLwd=.4,pTxt=12)
  
  # Set animation options
  # Careful with interval and nmax combo !!
  ani.options(interval = 0.2,
              nmax=100,
              loop = 1,
              autobrowse=FALSE,
              ani.width=1700,
              ani.height=1000)
  
  saveGIF({
    rtEnd = cyc
    
    #Loop through time cycles
    for(i in seq(1,rtEnd,by=nstep)){
      
      #p <- plotReplication(dfWithOri,i,pScale=2,pSz=1,pLwd=.7,pTxt=12)
      p <- plotReplicationForAnimationWG(dfWithOri,i,pScale=2,pSz=.6,pLwd=.4,pTxt=12)
      
      print(p)
      
    }#close the for loop
    
  }, convert = 'convert',
  img.name=paste0('allImg_',name),
  movie.name = paste0(name,'.gif')) 
  #close the animation builder
}

# WHOLE GENOME ########################################
makeRTGIF <-function(dfModel,cs="",name='repMovie',cyc=550,nstep=5){
  
  plotRTForAnimation(dfModel,1,cs=cs,pScale=3,pSz=1,pLwd=.7,pTxt=20)
  
  ## chr19 will be 7 wide 
  ## calibrate others to mm10 chr19 (61.42119 Mb)... 
  relativeWidth <- (max(dfModel$mod$model$trueCoord/1000000)+1)/62.42119
  
  # Set animation options
  # Careful with interval and nmax combo !!
  ani.options(interval = 0.2,
              nmax=100,
              loop = 1,
              autobrowse=FALSE,
              ani.width=1200*relativeWidth,
              ani.height=300)
  
  saveGIF({
    rtEnd = cyc
    
    #Loop through time cycles
    for(i in seq(1,rtEnd,by=nstep)){
      
      #p <- plotReplication(dfWithOri,i,pScale=2,pSz=1,pLwd=.7,pTxt=12)
      #p <- plotRTForAnimation(dfModel,i,cs=cs,pScale=2,pSz=.6,pLwd=.4,pTxt=12)
      p <- plotRTForAnimation(dfModel,i,cs=cs,pScale=3,pSz=1,pLwd=.7,pTxt=20)
      
      print(p)
      
    }#close the for loop
    
  }, convert = 'convert',
  img.name=paste0('allRTModImg_',name),
  movie.name = paste0(name,'.gif')) 
  #close the animation builder
}

# WHOLE GENOME ########################################
makeRTGIFWG <-function(dfModel,name='repMovie',cyc=550,nstep=5){
  
  plotRTForAnimationWG(dfModel,1,pScale=2,pSz=.6,pLwd=.4,pTxt=12)
  
  # Set animation options
  # Careful with interval and nmax combo !!
  ani.options(interval = 0.2,
              nmax=100,
              loop = 1,
              autobrowse=FALSE,
              ani.width=1700,
              ani.height=1000)
  
  saveGIF({
    rtEnd = cyc
    
    #Loop through time cycles
    for(i in seq(1,rtEnd,by=nstep)){
      
      #p <- plotReplication(dfWithOri,i,pScale=2,pSz=1,pLwd=.7,pTxt=12)
      p <- plotRTForAnimationWG(dfModel,i,pScale=2,pSz=.6,pLwd=.4,pTxt=12)
      
      print(p)
      
    }#close the for loop
    
  }, convert = 'convert',
  img.name=paste0('allRTModImg_',name),
  movie.name = paste0(name,'.gif')) 
  #close the animation builder
}

#######################################################
#######################################################
simRT_drawGIF    <- function(myModel      = NULL,
                             myModelFile  = NULL,
                             chrom        = 'WG',
                             from         = NULL,
                             to           = NULL,
                             cellNum      = 1,
                             nCycles      = 550,
                             nStep        = 5,
                             gifName      = 'tmpGIFname'){
  
  if (is.null(myModel)){
    if (chrom == 'WG'){
      mData <- getModelDataForGIFWG(cellNum=cellNum,
                                    modelPath=myModelFile)
    }else{
      mData <- getModelDataForGIF(chrom=chrom, 
                                  from=from,
                                  to=to,
                                  cellNum=cellNum,
                                  modelPath=myModelFile)
      
      ## Fix stupid line @ end
      mData <- mData[mData$pos != max(mData$pos),]
    }
  }else{
    if (chrom == 'WG'){
      mData <- getModelDataForGIFWG(cellNum=cellNum,
                                    model=myModel)
    }else{
      mData <- getModelDataForGIF(chrom=chrom, 
                                  from=from,
                                  to=to,
                                  cellNum=cellNum,
                                  model=myModel)
      
      ## Fix stupid line @ end
      mData <- mData[mData$pos != max(mData$pos),]
    }
  }
  
  if (chrom == 'WG'){
    makeReplicationGIFWG(mData,gifName,nCycles,nStep)
  }else{
    makeReplicationGIF(mData,gifName,nCycles,nStep)
  }

}

#######################################################
#######################################################
simRT_drawAllGIFs    <- function(myModel      = NULL,
                                 myModelFile  = NULL,
                                 chrom        = 'WG',
                                 from         = NULL,
                                 to           = NULL,
                                 cellNum      = 1,
                                 nCycles      = NULL,
                                 nStep        = 5,
                                 gifName      = 'tmpGIFname'){
  
  if (is.null(myModel)){
    load(myModelFile)
  }else{
    modelData <- myModel
  }
    
  maxPerCell <- max(modelData$mod$repPerCell)
  maxPerMin  <- dim(modelData$mod$repPerMin)[2]
  maxCycles  <- min(maxPerCell,maxPerMin)
  
  if (is.null(nCycles) | nCycles > maxCycles){
    nCycles <- maxCycles
  }
  
  for (chrom2use in paste0("chr",c(19,13,11,3))){
    
    makeRTGIF(modelData,cyc = nCycles, nstep = nStep, 
              cs=chrom2use,name=paste0(gifName,"_",chrom2use,"_popRT"))
    
    for (cell in 1:20){
      simRT_drawGIF(myModel = modelData,
                    chrom   = chrom2use,
                    cellNum = cell,
                    nCycles = nCycles,
                    nStep   = nStep,
                    gifName = paste0(gifName,"_",chrom2use,"_cell-",cell))
    }
  }
  
  simRT_drawGIF(myModel = modelData,
                chrom   = 'WG',
                cellNum = cell,
                nCycles = nCycles,
                nStep   = nStep,
                gifName = paste0(gifName,"_WG_cell-",cell))
}


###########################################################
###########################################################
# plotAnimationAllCS <- function(myAniModel,
#                                chrom,
#                                simLen=20,
#                                varwidth=FALSE,
#                                name='animation'){
#   
#   theme7point() 
#   
#   modelPlot <- myAniModel
#   model20M  <- modelPlot[modelPlot$from > 25000000 & modelPlot$from < 45000000,]
#   modelZoom <- modelPlot[modelPlot$from > 36500000 & modelPlot$from < 38500000,]
#   
#   pCS   <- ggplot(modelPlot,
#                   aes(x=from/1000000,y=sim,
#                       size=type,
#                       shape=type,
#                       color=type,
#                       fill=type,
#                       frame=title)) + 
#     geom_hline(yintercept=2,
#                color='grey60',
#                alpha=.2,
#                lwd=1) +
#     geom_point() + 
#     coord_cartesian(ylim=c(1.85,2.025),expand = FALSE) + 
#     scale_size_manual(values=c(1.6,6)) + 
#     scale_shape_manual(values=c(21,21)) + 
#     scale_color_manual(values = c('black','red')) + 
#     scale_fill_manual(values = c('black','red')) + 
#     theme(axis.text.y = element_blank(),
#           axis.line.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           panel.border=element_blank(),
#           axis.line.x=element_line(color='black'),
#           legend.position='top',
#           legend.title=element_blank(),
#           legend.direction = 'horizontal',
#           plot.title = element_text(size=30,hjust=.5)) + 
#     ylab('') + xlab('Position (Mb)') + 
#     transition_time(cycle) + 
#     labs(title = 'Elapsed: {round(frame_time*10/60,0)} hr') 
#   
#   p20M   <- ggplot(model20M,
#                    aes(x=from/1000000,y=sim,
#                        size=type,
#                        shape=type,
#                        color=type,
#                        fill=type,
#                        frame=title)) + 
#     geom_hline(yintercept=2,
#                color='grey60',
#                alpha=.2,
#                lwd=1) +
#     geom_point() + 
#     coord_cartesian(ylim=c(1.85,2.025),expand = FALSE) + 
#     scale_size_manual(values=c(1.6,6)) + 
#     scale_shape_manual(values=c(21,21)) + 
#     scale_color_manual(values = c('black','red')) + 
#     scale_fill_manual(values = c('black','red')) + 
#     theme(axis.text.y = element_blank(),
#           axis.line.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           panel.border=element_blank(),
#           axis.line.x=element_line(color='black'),
#           legend.position='top',
#           legend.title=element_blank(),
#           legend.direction = 'horizontal',
#           plot.title = element_text(size=30,hjust=.5)) + 
#     ylab('') + xlab('Position (Mb)') + 
#     transition_time(cycle) + 
#     labs(title = 'Elapsed: {round(frame_time*10/60,0)} hr') 
#   
#   pZoom  <- ggplot(modelZoom,
#                    aes(x=from/1000000,y=sim,
#                        size=type,
#                        shape=type,
#                        color=type,
#                        fill=type,
#                        frame=title)) + 
#     geom_hline(yintercept=2,
#                color='grey60',
#                alpha=.2,
#                lwd=1) +
#     geom_point() + 
#     coord_cartesian(ylim=c(1.85,2.025),expand = FALSE) + 
#     scale_size_manual(values=c(1.6,6)) + 
#     scale_shape_manual(values=c(21,21)) + 
#     scale_color_manual(values = c('black','red')) + 
#     scale_fill_manual(values = c('black','red')) + 
#     theme(axis.text.y = element_blank(),
#           axis.line.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           panel.border=element_blank(),
#           axis.line.x=element_line(color='black'),
#           legend.position='top',
#           legend.title=element_blank(),
#           legend.direction = 'horizontal',
#           plot.title = element_text(size=30,hjust=.5)) + 
#     ylab('') + xlab('Position (Mb)') + 
#     transition_time(cycle) + 
#     labs(title = 'Elapsed: {round(frame_time*10/60,0)} hr') 
#   
#   if (varwidth){
#     myWidth <- 2400 * (max(modelPlot$from)/200000000) 
#   }else{
#     myWidth <- 1600
#   }
#   
#   xCS   <- animate(pCS,width=myWidth,height=240,duration = simLen, renderer = gifski_renderer(loop = F))
#   x20   <- animate(p20M,width=myWidth,height=240,duration = simLen, renderer = gifski_renderer(loop = F))
#   xZoom <- animate(pZoom,width=myWidth,height=240,duration = simLen, renderer = gifski_renderer(loop = F))
#   
#   anim_save(filename = getIMGname(paste0(name,'_',chrom,'_wholeCS'),type = 'GIF'),animation = xCS)
#   anim_save(filename = getIMGname(paste0(name,'_',chrom,'_20Mb'),type = 'GIF'),animation = x20)
#   anim_save(filename = getIMGname(paste0(name,'_',chrom,'_Zoom'),type = 'GIF'),animation = xZoom)
# }


