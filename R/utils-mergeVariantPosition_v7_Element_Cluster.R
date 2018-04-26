#' mergeVariantPositionElementClusterV2
#'
#' @param datDF list of geneDF data frame
#' @param sizeFactor sizeFactor vector, a parameter only used in pancancer calculation. A sizeFactor for total mutations counts normalization among multiple cancer types
#' @param countsCutOff default is two, a parameter only used in pancancer calculation. Scale the mutation occurence value when the counts is equal or above countsCutOff
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

mergeVariantPositionElementClusterV2<-function(datDF,sizeFactor,countsCutOff=2){
  
   nonClusterKeyword<-"0"
   datDFbeyondCluster<-datDF[datDF$cluster %in% nonClusterKeyword,]
   
   mergedDFbeyondCluster<-{}
   
   if(FALSE){
   #if(nrow(datDFbeyondCluster)){
     
       recurrenceVectorBeyondCluster<-table(datDFbeyondCluster$posIndex)
       positionVectorBeyondCluster<-names(recurrenceVectorBeyondCluster)
       
       mergedDFbeyondCluster<-lapply(1:length(positionVectorBeyondCluster), function(x) {
         
         #cat(sprintf("%s",x))
         
         datDFtmp<-datDFbeyondCluster[datDFbeyondCluster$posIndex==positionVectorBeyondCluster[x],]
         
         tmp2DF<-datDFtmp
         tmp3DF<-tmp2DF[,c(1:16,18:19)]
         #tmp3DF<-cbind(tmp2DF[,c(1,3:12)],tmp2DF[,15:24])
         
         
         tumorName<-datDFtmp$tumorType
         score<-datDFtmp$score
         
         
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
         
         posIndex<-datDFtmp$posIndex
         pickedIdx<-ceiling(median(order(posIndex)))
         posIndex<-posIndex[pickedIdx]
         
         #posIndexWithClusterID<-paste(posIndex,nonClusterKeyword,sep=":")
         posIndexWithClusterID<-nonClusterKeyword
         
         # merged into one line
         tmp3DFmerged<-tmp3DF[tmp3DF$posIndex %in% posIndex,][1,]
         
         triMutIndex<-datDFtmp$triMutIndex
         signalValue<-datDFtmp$signalValue[1]
         
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
         
         #posIndex<-posIndexWithClusterID
         occurence<-recurrenceVectorBeyondCluster[x]
         categoryCounts<-strName
         data.frame(tmp3DFmerged,posIndexWithClusterID,occurence,compositeScore,compositeScoreScaled,categoryCounts,signalValue,stringsAsFactors = FALSE)
         
         
       })
       
       mergedDFbeyondCluster<-do.call(rbind.data.frame,mergedDFbeyondCluster)
       #colnames(mergedDFbeyondCluster)<-c("posIndex","occurence","compositeScore","compositeScoreScaled","categoryCounts","signalValue")
   }
          
   ######
   #datDFinCluster<-datDF[!(datDF$cluster %in% nonClusterKeyword),]
   
   #mergedDFinCluster<-{}
   
   if(FALSE){
   #if(nrow(datDFinCluster)>0){
     
       recurrenceVectorInCluster<-table(datDFinCluster$cluster)
       clusterVectorInCluster<-names(recurrenceVectorInCluster)
       
       mergedDFinCluster<-lapply(1:length(clusterVectorInCluster), function(x) {
         
         datDFtmp<-datDFinCluster[datDFinCluster$cluster==clusterVectorInCluster[x],]
         
         tmp2DF<-datDFtmp
         tmp3DF<-tmp2DF[,c(1:16,18:19)]
         #tmp3DF<-cbind(tmp2DF[,c(1,3:12)],tmp2DF[,15:24])
         
         tumorName<-datDFtmp$tumorType
         score<-datDFtmp$score
         
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
         
         posIndex<-datDFtmp[datDFtmp$cluster==clusterVectorInCluster[x],]$posIndex
         pickedIdx<-ceiling(median(order(posIndex)))
         posIndex<-posIndex[pickedIdx]
         
         # merged into one line
         tmp3DFmerged<-tmp3DF[tmp3DF$posIndex %in% posIndex,][1,]
         
         #posIndexWithClusterID<-paste(posIndex,clusterVectorInCluster[x],sep=":")
         posIndexWithClusterID<-clusterVectorInCluster[x]
         
         triMutIndex<-datDFtmp[datDFtmp$cluster==clusterVectorInCluster[x],]$triMutIndex
         signalValue<-datDFtmp[datDFtmp$cluster==clusterVectorInCluster[x],]$signalValue[1]
         
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
         
         #posIndex<-posIndexWithClusterID
         occurence<-recurrenceVectorInCluster[x]
         categoryCounts<-strName
         data.frame(tmp3DFmerged,posIndexWithClusterID,occurence,compositeScore,compositeScoreScaled,categoryCounts,signalValue,stringsAsFactors = FALSE)
         
         
       })
       
       mergedDFinCluster<-do.call(rbind.data.frame,mergedDFinCluster)
       #colnames(mergedDF)<-c("posIndex","occurence","compositeScore","compositeScoreScaled","categoryCounts","signalValue")
       
   }
         
   ######
   datDFinCluster<-datDF[!(datDF$cluster %in% nonClusterKeyword),]
   
   mergedDFinCluster<-{} 

   if(nrow(datDFinCluster)){
     
     recurrenceVectorInCluster<-table(datDFinCluster$posIndex)
     positionVectorInCluster<-names(recurrenceVectorInCluster)
     
     mergedDFinCluster<-lapply(1:length(positionVectorInCluster), function(x) {
       
       #cat(sprintf("%s",x))
       
       datDFtmp<-datDFinCluster[datDFinCluster$posIndex==positionVectorInCluster[x],]
       
       clusterID<-datDFtmp$cluster[1]
       tmp2DF<-datDFtmp
       tmp3DF<-tmp2DF[,c(1:16,18:19)]
       #tmp3DF<-cbind(tmp2DF[,c(1,3:12)],tmp2DF[,15:24])
       #tmp3DF$oldScore<-tmp2DF$oldScore
       
       
       tumorName<-datDFtmp$tumorType
       score<-datDFtmp$score
       
       
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
       
       posIndex<-datDFtmp$posIndex
       pickedIdx<-ceiling(median(order(posIndex)))
       posIndex<-posIndex[pickedIdx]
       
       #posIndexWithClusterID<-paste(posIndex,clusterID,sep=":")
       posIndexWithClusterID<-clusterID
       
       # merged into one line
       tmp3DFmerged<-tmp3DF[tmp3DF$posIndex %in% posIndex,][1,]
       
       
       triMutIndex<-datDFtmp$triMutIndex
       signalValue<-datDFtmp$signalValue[1]
       
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
       
       #posIndex<-posIndexWithClusterID
       occurence<-recurrenceVectorInCluster[x]
       categoryCounts<-strName
       data.frame(tmp3DFmerged,posIndexWithClusterID,occurence,compositeScore,compositeScoreScaled,categoryCounts,signalValue,stringsAsFactors = FALSE)
       
       
     })
     
     mergedDFinCluster<-do.call(rbind.data.frame,mergedDFinCluster)
     #colnames(mergedDFbeyondCluster)<-c("posIndex","occurence","compositeScore","compositeScoreScaled","categoryCounts","signalValue")
   }
   
   ######
   
   #mergedDF<-rbind(mergedDFbeyondCluster,mergedDFinCluster)
   mergedDF<-mergedDFinCluster
   
   
   ######
   
   return(mergedDF)
}