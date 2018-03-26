#' Extract overlapping mutations from element annotation bed file
#'
#' @param elementBedfile The data frame for element annotated bed file 
#' @param reducedFunseqOutput reducedFunseq2 data frame
#' @param debugMode TRUE or FALSE
#'
#' @return reducedFunseq2 data frame
#'
#' @examples
#' #date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' selectedSampleId<-NA
#' #worDir<-getwd()
#' mutSig2CVthreshold<-0.1
#' rareMutationUpperLimit<-0.3
#' rareMutationLowerLimit<-0.1
#' rareMutationFreq<-0.02
#'
#' #runNetBox2(dataDir,cancerType,
#' #           mutationList,ampGeneList,delGeneList,epiSilencedList,
#' #           mutationFreq,ampGeneFreq,delGeneFreq,epiSilencedFreq,
#' #           pathwayCommonsDb,directed,
#' #           linkerPValThreshold,communityDetectionMethod,
#' #           keepIsolatedNodes,verbose=TRUE)
#'
#' @concept CNCDriver
#' @export
#' @importFrom data.table setkeyv
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table foverlaps



extractElementMutations<-function(elementBedfile,reducedFunseqOutput,debugMode=FALSE){

    testDT<-data.table(reducedFunseqOutput)
    setkeyv(testDT,c("chr","posStart","posEnd"))
    
    #######
    #
    # keep mutation on elementBedfile
    # 
    ######
    
    #elementBedfile<-read.table(elementBedfileName,sep="\t",header=FALSE,stringsAsFactors = FALSE)
    #colnames(elementBedfile)<-c("chr","posStart","posEnd","name")
    #colnames(dat)<-c("chr","posStart","posEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
    
    
    ##
    ## only take the first four columns in the bed file
    ##
    ## chr posStart posEnd       name                
    ## chr1    68891  69090      OR4F5 
    ## chr1   139310 139579 AL627309.1
    
    dat<-elementBedfile[,c(1:4)]
    
    #dat<-convertBED12format(dat,useCores)
    
    # add +1 to correct 0-index start position from bed format
    dat$posStart<-dat$posStart+1
    
    ####
    
    #feature<-strsplit(dat$elementID,"::")
    #featureDf<-do.call(rbind,feature)
    #dat$elementWholeName<-featureDf[,3]
    
    dacDT<-as.data.table(dat)
    
    setkeyv(dacDT,c("chr","posStart","posEnd"))
    
    #feature<-strsplit(dacDT$elementID,"::")
    #featureDf<-do.call(rbind,feature)
    #keepElementName<-featureDf[,2]
    
    mutKeepable<-foverlaps(testDT,dacDT,nomatch=0)
    #numOfWhitelistMutations<-nrow(mutKeepable)
    
    #cat(sprintf("there are %s mutations on the regions of element bed file\n",numOfWhitelistMutations))
    
    
    
      
      originalNames<-colnames(reducedFunseqOutput)
      
      tmpDF<-as.data.frame(mutKeepable)
      ncBlockDF<-tmpDF[,c(1:4)]
      ncBlockDF$elementName<-paste(ncBlockDF$chr,paste(ncBlockDF$posStart,ncBlockDF$posEnd,sep="-"),sep=":")
      
      chrDF<-tmpDF[,1]
      #sampleIDvector<-tmpDF[,7:8]
      secondDF<-tmpDF[,c(5:ncol(tmpDF))]
      
      subsetDF<-data.frame(chrDF,secondDF,stringsAsFactors = FALSE)
      colnames(subsetDF)<-originalNames
      subsetDF<-cbind(subsetDF,ncBlockDF[,c(4:5)])
      reducedFunseqOutput<-subsetDF
      
    
    
    remove(dacDT)
    remove(mutKeepable)
    remove(testDT)
    
    #dd<-reducedFunseqOutputNCDS
    #dd$index<-paste(dd$chr,paste(dd$posStart,dd$posEnd,sep="-"),sep=":")
    
    #####
    
    reducedFunseqOutput$index<-paste(reducedFunseqOutput$sampleID,reducedFunseqOutput$chr,paste(reducedFunseqOutput$posStart,reducedFunseqOutput$posEnd,sep="-"),sep=":")
    
    duplicatedLine<-duplicated(reducedFunseqOutput$index)
    
    reducedFunseqOutput<-reducedFunseqOutput[!duplicatedLine,]
    
    reducedFunseqOutput<-reducedFunseqOutput[,c(1:26)]
    
    numOfWhitelistMutations<-nrow(reducedFunseqOutput)
    cat(sprintf("there are %s mutations on the regions of element bed file\n",numOfWhitelistMutations))
    
    return(reducedFunseqOutput)
######

}

