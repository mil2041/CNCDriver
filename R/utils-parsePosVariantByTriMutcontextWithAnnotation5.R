#' Parse Pos Variant By TriMutContext With Annotation
#'
#' @param geneDFunique geneDFunique data frame
#' @param mutationDistMatrix mutationDistMatrix data frame
#' @param useCore default is one
#'
#' @return variantTriMutCategoryParsed data frame
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
#' @importFrom stringr str_extract_all
#' @importFrom parallel mclapply
#' @importFrom plyr rbind.fill

parsePosVariantByTriMutContextWithAnnotation5<-function(geneDFunique,mutationDistMatrix,useCores=1){
  
  #stringVector<-a1$categoryCounts
  stringVector<-geneDFunique$categoryCounts
  
  categoryMatch<-gsub("[[0-9]+]","",stringVector)
  
  counts<-str_extract_all(stringVector,"([0-9]+)")
  countsRatio<-sapply(1:length(counts), function(x){paste(counts[[x]],collapse=":")})
  counts<-sapply(1:length(counts), function(x){sum(as.numeric(counts[[x]]))})
  
  tmpStr<-strsplit(stringVector,"\\,")
  
  tmpStr<-mclapply(1:length(tmpStr), function(x){
    tmp<-strsplit(tmpStr[[x]],":")
    tmp2<-sapply(1:length(tmp),function(y){tmp[[y]][1]})
    tumorName<-paste(unique(tmp2),collapse=",")
    #tumorName<-paste(tmp2,collapse=",")
    numOfTumorType<-length(unique(tmp2))
    tmp3<-str_extract_all(tmp,"[ACGT][ACGT]@[ACGT]+.[ACGT]+")
    categoryName<-paste(unique(tmp3),collapse=",")
    numOfCategory<-length(unique(tmp3))
    
    data.frame(tumorName,numOfTumorType,categoryName,numOfCategory,stringsAsFactors = FALSE)
  
  },mc.cores=useCores)
  
  tmpStr<-rbind.fill(tmpStr)
  triMutContextAnnotation<-data.frame(tmpStr,countsRatio,counts,categoryMatch,stringsAsFactors = FALSE)
  
  geneDFuniqueSimple<-data.frame(geneDFunique$compositeScore,geneDFunique$compositeScoreScaled,geneDFunique$posIndex,geneDFunique$signalValue,geneDFunique$geneSymbol,stringsAsFactors = FALSE)
  colnames(geneDFuniqueSimple)<-c("compositeScore","compositeScoreScaled","posIndex","signalValue","geneSymbol")
  
  str<-categoryMatch
  
  splitedDat<-mclapply(1:length(str), function(x){
  #splitedDat<-lapply(1:length(str), function(x){
    #cat(sprintf("iter %s\n",x))
    tmpStr<-unlist(strsplit(str[x],","))
    tmpStr2<-strsplit(tmpStr,":")
    tmpCounts<-strsplit(countsRatio[x],":")
    tmpDF3<-data.frame(do.call(rbind,tmpStr2),tmpCounts,stringsAsFactors = FALSE)
    colnames(tmpDF3)<-c("tumorType","categoryName","counts")
    return(tmpDF3)
  #})
  },mc.cores=useCores)  
      
  posCategoryFreq<-mclapply(1:length(splitedDat), function(y){
  #posCategoryFreq<-lapply(1:length(splitedDat), function(y){
  #  cat(sprintf("iter %s\n",y))
    categoryFreq<-sapply(1:nrow(splitedDat[[y]]), function(z){
          selectedCol<-which(colnames(mutationDistMatrix) %in% splitedDat[[y]][z,1])
          freq<-mutationDistMatrix[splitedDat[[y]][z,2],selectedCol]
    })
    
    splitedDat[[y]]$prob<-categoryFreq
    dat1<-splitedDat[[y]]
    weightedFreq<-sum(as.numeric(dat1$counts)*dat1$prob)/sum(as.numeric(dat1$counts))
    return(weightedFreq)
  #})  
  },mc.cores=useCores)
    
  posCategoryFreq<-unlist(posCategoryFreq)
  
  #result<-data.frame(triMutContextAnnotation,geneDFuniqueSimple,posCategoryFreq,stringsAsFactors = FALSE)
  result<-data.frame(triMutContextAnnotation,geneDFunique,posCategoryFreq,stringsAsFactors = FALSE)
  
  return(result)
}