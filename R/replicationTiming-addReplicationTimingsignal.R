#' Add replication timing signal value
#'
#' @param reducedFunseqOutput reducedFunseqOutput data frame
#' @param replicationTimingDF replicationTimingDF data frame
#' @param useCores Default is one, number of cpu to use
#'
#' @return reducedFunseqOutput data frame
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

addReplicationTimingSignal<-function(reducedFunseqOutput,replicationTimingDF,useCores=1){
  
  tmpDF<-reducedFunseqOutput
  
  compositeScorePosGRanges<-GRanges(seqnames=tmpDF$chr,
                                    ranges=IRanges(start=(tmpDF$posStart),
                                                   end=(tmpDF$posEnd)))
                                    
                                    #geneSymbol=tmpDF$geneSymbol,
                                    #ref=tmpDF$ref,
                                    #alt=tmpDF$alt,
                                    #score=tmpDF$score,
                                    #occurence=tmpDF$occurence,
                                    #compositeScore=tmpDF$compositeScore)
  
  
  
  replicationTimingGRanges<-GRanges(seqnames=replicationTimingDF$chr,
                                    ranges=IRanges(start=(replicationTimingDF$start),
                                                   end=(replicationTimingDF$end)),
                                    name=replicationTimingDF$name,
                                    signalValue=replicationTimingDF$signalValue)
  
  
  hits<-findOverlaps(compositeScorePosGRanges,replicationTimingGRanges)
  
  matchDF<-data.frame(queryHits(hits),subjectHits(hits),stringsAsFactors = FALSE)
  colnames(matchDF)<-c("query","subject")
  
  rtSignal<-rep(NA,length(compositeScorePosGRanges))
  signalVector<-mcols(replicationTimingGRanges)$signalValue[matchDF$subject]
  tmpDF$signalValue<-rtSignal
  tmpDF$signalValue[matchDF$query]<-signalVector
  
  return(tmpDF)
  
}