#' mergeVariantPosition
#'
#' @param geneDF list of geneDF data frame
#' @param sizeFactor sizeFactor data frame
#' @param countsCutOff default is two
#'
#' @return collapsedDF data frame
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

mergeVariantPosition<-function(datDF,sizeFactor,countsCutOff=2){
   recurrenceVector<-table(datDF$posIndex)
   positionVector<-names(recurrenceVector)

   ##
   ## amplify count >=2
   ## reduce count < 2  
   ##
    
   mergedDF<-lapply(1:length(positionVector), function(x) {
     
     #cat(sprintf("%s",x))
     tumorName<-datDF[datDF$posIndex==positionVector[x],]$tumorType
     score<-datDF[datDF$posIndex==positionVector[x],]$score
     
     
     
     if(length(score)>=countsCutOff){
     
       sizeFactorNormalized<-sizeFactor/max(sizeFactor)
       #sizeFactorNormalized<-sizeFactor
     
       recurrenceNormalizedFactor<-1/sizeFactorNormalized[tumorName]
     
       rescaledScore<-recurrenceNormalizedFactor*score
     
     }else{
       rescaledScore<-score/sizeFactor[tumorName]
     }
     
     #triMutFreq<-datDF[datDF$posIndex==positionVector[x],]$triMutFreq
     #reScaleFactor<-length(score)
     #rescaleScore<-score*(1-triMutFreq)^reScaleFactor
     
     
     
     triMutIndex<-datDF[datDF$posIndex==positionVector[x],]$triMutIndex
     signalValue<-datDF[datDF$posIndex==positionVector[x],]$signalValue[1]
     
     tumorNameMerged<-paste(tumorName,triMutIndex,sep=":")
     
     
     compositeScore<-sum(score)   
     compositeScoreScaled<-sum(rescaledScore)  
     tumorCounts<-table(tumorNameMerged)
     tumorName<-names(tumorCounts)
     tumorCounts<-as.numeric(tumorCounts)
     
     strName<-sapply(1:length(tumorName), function(x){
        paste(tumorName[x],"[",tumorCounts[x],"]",sep="")      
     })
     
     strName<-paste(strName,collapse=",")
     
     data.frame(positionVector[x],recurrenceVector[x],compositeScore,compositeScoreScaled,strName,signalValue,stringsAsFactors = FALSE)
     
     
   })
   
   mergedDF<-do.call(rbind.data.frame,mergedDF)
   colnames(mergedDF)<-c("posIndex","occurence","compositeScore","compositeScoreScaled","categoryCounts","signalValue")
   
   return(mergedDF)
}