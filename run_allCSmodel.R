#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-b", "--RTbg"),         type="character", default=NULL,   help="RT BEDgraph file",        metavar="character"),
  make_option(c("-c", "--totCells"),     type="numeric",   default=500,    help="Total cells",             metavar="numeric"),
  make_option(c("-d", "--DSBs"),         type="character", default=NULL,   help="DSBs BEDgraph file",      metavar="character"),
  make_option(c("-e", "--optimalT"),     type="numeric",   default=NULL,   help="Optimal Model Time",      metavar="numeric"),
  make_option(c("-f", "--optimalPC"),    type="numeric",   default=NULL,   help="Optimal PC",              metavar="numeric"),
  make_option(c("-g", "--genome"),       type="character", default="mm10", help="Genome",                  metavar="character"),
  make_option(c("-i", "--id"),           type="character", default=3000,   help="Unique sampleID",         metavar="character"),
  make_option(c("-m", "--savemodel"),    type="logical",   default=FALSE,  help="save model",              metavar="logical"),
  make_option(c("-n", "--oriPerMb"),     type="numeric",   default=0.2,    help="Ori Per Mb",              metavar="numeric"),
  make_option(c("-o", "--oriBed"),       type="character", default=NULL,   help="Origins BED file",        metavar="character"),
	make_option(c("-q", "--scriptpath"),   type="character", default=getwd(),help="Path for support scripts",metavar="character"),
  make_option(c("-r", "--recycle"),      type="logical",   default=TRUE,   help="use recycling",           metavar="logical"),
  make_option(c("-s", "--strength"),     type="logical",   default=TRUE,   help="use strength",            metavar="logical"),
  make_option(c("-t", "--time"),         type="numeric",   default=3000,   help="Model run time (sec)",    metavar="numeric"),
  make_option(c("-u", "--outname"),      type="character", default=NULL,   help="Output file name prefix", metavar="character"),
  make_option(c("-v", "--excludeChrom"), type="character", default=NULL,   help="Exclude chrom",           metavar="character"), 
	make_option(c("-w", "--randseed"),     type="numeric",   default=NULL,   help="Random seed for picking origins", metavar="numeric"), 
  make_option(c("-x", "--tstChrom"),     type="character", default=NULL,   help="Test chrom",              metavar="character"),
  make_option(c("-y", "--imgOut"),       type="character", default=NULL,   help="Image output folder",     metavar="character"),
  make_option(c("-z", "--outRDs"),       type="logical",   default=FALSE,  help="Output RDS files",        metavar="logical")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Sys.setenv(RTSCRIPTS = getwd())

if (is.null(opt$imgOut)){
	imgOutputFolder <- paste0(getwd(),'/')
}else{
	imgOutputFolder <- paste0(opt$imgOut,'/')
}

source(paste0(Sys.getenv('RTSCRIPTS'),'/RDCOstandardfunctions.R'))
source(paste0(Sys.getenv('RTSCRIPTS'),'/repliSim_loadModules.R'))
source(paste0(Sys.getenv('RTSCRIPTS'),'/simulateReplicationTiming_allCS_190806.R'))

thisDate <- gsub(" ","_",gsub(":","",gsub("-",replacement = "",Sys.time())))

print(do.call(paste, list(names(opt), " = ", opt) ))

if (is.null(opt$excludeChrom)){
  chromStr <- 'allCS'
}else{
  chromStr <- paste0("no",opt$excludeChrom)
  
}

if (is.null(opt$outname)){
	myName <- paste0(paste(opt$id,
	                       thisDate,
	                       opt$time,
	                       opt$oriPerMb,
	                       opt$totCells,
	                       opt$genome,
	                       chromStr,
	                       "s",opt$strength,
	                       "r",opt$recycle,
	                       sep="_"),'.Rdata')

	myModel <- paste0(paste(opt$id,
		                thisDate,
		                opt$time,
		                opt$oriPerMb,
		                opt$totCells,
		                opt$genome,
		                chromStr,
		                "s",opt$strength,
		                "r",opt$recycle,
		                sep="_"),'.model.rds')

	myFigData <- paste0(paste(opt$id,
		                  thisDate,
		                  opt$time,
		                  opt$oriPerMb,
		                  opt$totCells,
		                  opt$genome,
		                  chromStr,
		                  "s",opt$strength,
		                  "r",opt$recycle,
		                  sep="_"),'.figdata.rds')

	myStatData <- paste0(paste(opt$id,
		                   thisDate,
		                   opt$time,
		                   opt$oriPerMb,
		                   opt$totCells,
		                   opt$genome,
		                   chromStr,
		                   "s",opt$strength,
		                   "r",opt$recycle,
		                   sep="_"),'.stats.txt')
}else{
	myName     <- paste0(opt$outname,'_',opt$genome,'.Rdata')
	myModel    <- paste0(opt$outname,'_',opt$genome,'.model.rds')
	myFigData  <- paste0(opt$outname,'_',opt$genome,'.figdata.rds')
	myStatData <- paste0(opt$outname,'_',opt$genome,'.stats.txt')
}

print(paste0("Stats file will output as: ",myStatData))

rtMod <- simRT_buildModel(repRate         = 1, ## Kb per minute
                          originsFile     = opt$oriBed,
                          timePassed      = opt$time,
                          oriPerMb        = opt$oriPerMb,
                          numCells        = opt$totCells,
                          genome          = opt$genome,
                          chrom           = 'ALL',
                          excludeChrom    = opt$excludeChrom,
                          useOriStrength  = opt$strength,
                          recycle         = opt$recycle,
                          dsbFile         = opt$DSBs,
                          plotAnimations  = FALSE,
                          useSmoothing    = TRUE,
                          IODtime         = 50,
                          noNormalize     = FALSE,
                          repStep         = 10,
                          simOriSeq       = FALSE,
                          randseed        = opt$randseed,
                          verbose         = TRUE)


rtMod <- simRT_addBGdata(rtModel = rtMod,
                         RTbg         = opt$RTbg,
                         genome       = opt$genome,
                         isLog2       = TRUE,
                         noSexChrom   = TRUE,
                         RTname       = 'rtInit',
                         useSmoothing = TRUE)

print(paste0("90% = ",quantile(abs(rtMod$rtInit$smooth),0.9)))

if (is.null(opt$optimalT) || is.null(opt$optimalPC)){
	rtComp <- simRT_compare2Data(myModel      = rtMod,
		                     genome       = opt$genome,
		                     RTdata       = rtMod$rtInit$smooth,
		                     useSmoothing = TRUE,
		                     isLog2       = TRUE,
		                     noNormalize  = TRUE,
		                     timepoint    = NULL,
		                     pcRep        = NULL,
		                     noSexChrom   = TRUE,
		                     chrom2use    = opt$tstChrom,
		                     plotMe       = FALSE)

	iWholeGenome <- simRT_drawFit(myModel     = rtMod,
		      chrom2use   = NULL,
		      timepoint   = rtComp$bestStat$time,
		      pcRep       = rtComp$bestStat$pcRep,
		      myName      = myName,
		      noNormalize = TRUE)

	iCS11 <- simRT_drawFit(myModel     = rtMod,
		      chrom2use   = 'chr11',
		      timepoint   = rtComp$bestStat$time,
		      pcRep       = rtComp$bestStat$pcRep,
		      myName      = myName,
		      noNormalize = TRUE)
}else{
	opt$optimalT <- round(opt$optimalT/10)
	print(paste0("Skipping parameter fits; Time = ",opt$optimalT,"; PC = ",opt$optimalPC))
	
	plotCS <- 'chr6'
	## Because sometimes, we might exclude the chromosome we want to plot
	if (!is.null(opt$excludeChrom)){
	  for (cs2Use in c(plotCS,levels(rtMod$whatCS))){
	    if (!grepl(cs2Use,opt$excludeChrom)){
	      plotCS <- cs2Use
	      break
	    }
	  }
	}
	
	print(paste0("Using ",plotCS," for plot ... [",opt$excludeChrom,"]"))
	
	iWholeGenome <- simRT_drawFit(myModel     = rtMod,
		      chrom2use   = NULL,
		      timepoint   = opt$optimalT,
		      pcRep       = opt$optimalPC,
		      myName      = myName,
		      noNormalize = FALSE,
		      plotMe      = FALSE,
		      imgPath     = '.')

	iCS11 <- simRT_drawFit(myModel     = rtMod,
		      chrom2use   = plotCS,
		      timepoint   = opt$optimalT,
		      pcRep       = opt$optimalPC,
		      myName      = myName,
		      noNormalize = FALSE,
		      plotMe      = FALSE,
		      imgPath     = '.')

	rtMod$params$bestTime <- opt$optimalT
	rtMod$params$bestPC   <- opt$optimalPC

	#outSimRT  <- simRT_exportSimRT(myModel = rtMod, myName = myName)
  #outRealRT <- simRT_exportRealRT(myModel = rtMod, myName = myName)
  
	outSimRT  <- simRT_exportRTBG(myModel = rtMod, myName = myName, type = 'sim')
	outRealRT <- simRT_exportRTBG(myModel = rtMod, myName = myName, type = 'exp')
  
	rtComp <- list(bestStat=data.frame(time=opt$optimalT,
				  	   pcRep=opt$optimalPC));

	ggsave(paste0(myName,'.WG.png')   ,iWholeGenome$fig,width=7,height=7)
	ggsave(paste0(myName,'.WG.pdf')   ,iWholeGenome$fig,width=7,height=7)
	ggsave(paste0(myName,'.',plotCS,'.png'),iCS11$fig  ,width=6,height=2)
	ggsave(paste0(myName,'.',plotCS,'.pdf'),iCS11$fig  ,width=6,height=2)
}

#print("90% = [3]")
#print(quantile(abs(rtMod$rtInit$smooth),0.9))

if (opt$savemodel){
	modelData <- list(mod=rtMod,
	                  figData=rtComp)
	save(modelData,file = myName)
}

rtComp$time$raw$bg        <- opt$RTbg
rtComp$time$raw$ori       <- opt$oriBed
rtComp$time$raw$useStr    <- opt$strength
rtComp$time$raw$recycling <- opt$recycle

rtComp$byTime$raw$bg        <- opt$RTbg
rtComp$byTime$raw$ori       <- opt$oriBed
rtComp$byTime$raw$useStr    <- opt$strength
rtComp$byTime$raw$recycling <- opt$recycle
rtComp$byTime$raw$nSimCells <- opt$totCells
rtComp$byTime$raw$genome    <- opt$genome
rtComp$byTime$raw$dateStamp <- thisDate

if (opt$outRDs){
  saveRDS(rtMod,  file = myModel)
  saveRDS(rtComp, file = myFigData)
}

print(paste0("Writing stats file: ",myStatData))

write.table(x         = rtComp$byTime$raw, file  = myStatData,
            row.names = FALSE,             quote = FALSE,
            sep       = "\t")
