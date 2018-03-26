#' Calculate recurrence number
#'
#' @param datDF geneDFunique data frame
#'
#' @return results data frame
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
mergeVariantPositionSimple<-function(datDF){
  recurrenceVector<-table(datDF$posIndex)
  positionVector<-names(recurrenceVector)

  mergedDF<-lapply(1:length(positionVector), function(x) {

    #cat(sprintf("%s",x))
    tumorName<-datDF[datDF$posIndex==positionVector[x],]$tumorType
    score<-datDF[datDF$posIndex==positionVector[x],]$score
    #triMutIndex<-datDF[datDF$posIndex==positionVector[x],]$triMutIndex
    #signalValue<-datDF[datDF$posIndex==positionVector[x],]$signalValue[1]

    #tumorNameMerged<-paste(tumorName,triMutIndex,sep=":")


    compositeScore<-sum(score)
    tumorCounts<-table(tumorName)
    #tumorCounts<-table(tumorNameMerged)
    tumorName<-names(tumorCounts)
    tumorCounts<-as.numeric(tumorCounts)

    #strName<-sapply(1:length(tumorName), function(x){
    #  paste(tumorName[x],"[",tumorCounts[x],"]",sep="")
    #})

    #strName<-paste(strName,collapse=",")

    #data.frame(positionVector[x],recurrenceVector[x],compositeScore,strName,signalValue,stringsAsFactors = FALSE)
    data.frame(positionVector[x],recurrenceVector[x],compositeScore,stringsAsFactors = FALSE)

  })

  mergedDF<-do.call(rbind.data.frame,mergedDF)
  colnames(mergedDF)<-c("posIndex","occurence","compositeScore")

  return(mergedDF)
}
