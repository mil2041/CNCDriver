#' Make replication timing signal value for each regulatory element
#'
#' @param replicationTimingDF replicationTimingDF data frame
#' @param elementDF regulatory element data frame
#' @param useCores Default is one, number of cpu to use
#'
#' @return elementsRTdf data frame
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
#' @import IRanges
#' @import GenomicRanges
#' @importFrom parallel mclapply

makeElementReplicationTiming<-function(replicationTimingDF,elementDF,useCores=1){
  
  peaksGRanges=GRanges(seqnames=replicationTimingDF$chr,
                       ranges=IRanges(start=replicationTimingDF$start,end=replicationTimingDF$end),
                       name=replicationTimingDF$name,
                       signalValue=replicationTimingDF$signalValue)
  
  #####
  
  #elementsDF<-enhancerDF
  elementDF<-elementDF[,c(1:4)]
  colnames(elementDF)<-c("chr","posStart","posEnd","geneSymbol")
  
  elementsGRanges=GRanges(seqnames=elementDF$chr,
                            ranges=IRanges(start=elementDF$posStart,
                                           end=elementDF$posEnd),
                            geneSymbol=elementDF$geneSymbol )
  
  #######
  
  hits<-findOverlaps(elementsGRanges,peaksGRanges)
  
  #head(queryHits(hits))
  
  nameVector<-mcols(elementsGRanges)$geneSymbol[queryHits(hits)]
  signalVector<-mcols(peaksGRanges)$signalValue[subjectHits(hits)]
  tmpDF<-data.frame(nameVector,signalVector,stringsAsFactors = FALSE)
  
  nameVector<-unique(tmpDF$nameVector)
  
  replicationTimingSignal<-{}
  #for(name in nameVector){
  #for(i in 1:length(nameVector)){
  
  #useCores=6
  
  out<-mclapply(1:length(nameVector), function(x){
    
    
    #replicationTimingSignal[[name]]<-mean(tmpDF[tmpDF$nameVector %in% name,]$signalVector)
    cat(sprintf("process %s / %s \n",x,length(nameVector)))
    
    ## when a regulatory element region hit two bin of replitation timing in the genome, take mean value
    replicationTimingSignal[[nameVector[x]]]<-mean(tmpDF[tmpDF$nameVector %in% nameVector[x],]$signalVector)
    
  }, mc.cores=useCores)
  
  replicationTimingSignal<-unlist(out)
  names(replicationTimingSignal)<-nameVector
  
  elementsRTdf<-data.frame(names(replicationTimingSignal),replicationTimingSignal,stringsAsFactors = FALSE)
  colnames(elementsRTdf)<-c("elementName","signalValue")
  
  return(elementsRTdf)
  
}
