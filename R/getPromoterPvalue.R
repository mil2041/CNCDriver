#' Run CNCDriver promoter p-value calculation
#'
#' @param inputDir Where is annotated funseq2 result
#' @param tumorType study name
#' @param mutationType User provided mutated gene list
#' @param reSampleIter User provided re-sapling iteration numbers
#' @param seedNum User provided random number seed, default is 42
#' @param debugMode TRUE or FALSE
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
#' @importFrom plyr rbind.fill
#' @importFrom stats p.adjust
#' @importFrom utils write.table
getPromoterPvalue<-function(inputDir,tumorType,mutationType,reSampleIter=1000,seedNum=42,debugMode=FALSE){


#library(data.table)
#library(plyr)
######

#inputDir<-"~/work/Ekta_lab/Priyanka_project/fsigsnv_results"
#tumorType<-"Prostate"

workDir<-file.path(inputDir,tumorType)

#mutationType<-"CDS"
#seedNum<-42
#reSampleNum<-1000
reSampleNum<-reSampleIter
set.seed(seedNum)

if( !file.exists(paste(workDir,sep="/")) ){
  dir.create(paste(workDir,sep=""),recursive=TRUE)
}

########

filePath<-file.path(workDir,"input")
fileName<-paste("reducedFunseqOutputNCDS_",tumorType,".Rd",sep="")
fileName<-file.path(filePath,fileName)
load(fileName)

#####

cat(sprintf("Processing %s - %s\n",tumorType,mutationType))
reducedFunseqOutputNCDS<-reducedFunseqOutputNCDS[grep("Promoter",reducedFunseqOutputNCDS$GENE),]

#####
parseFunseqGeneField<-function(aa,keyword){
  geneSymbol<-{}
  #a1<-strsplit(aa,",")
  a1<-strsplit(aa,"\\)")
  #a1<-unlist(a1)

  for(i in 1:length(a1)){
    #cat(sprintf("gene:%s\t",i))
    idx <- which(grepl(keyword, a1[[i]]))
    a1[[i]]<-gsub("\\,","",a1[[i]])
    t3 <- strsplit(a1[[i]][idx], "\\(")
    #cat(sprintf("%s\n",t3[1][1]))
    geneSymbol[i]<-t3[[1]][1]

  }
  return(geneSymbol)
}

#####

reducedFunseqOutputNCDS$GENEparsed<-parseFunseqGeneField(aa=reducedFunseqOutputNCDS$GENE,keyword="Promoter")

gg<-strsplit(as.character(reducedFunseqOutputNCDS$NCDS),":",fixed=TRUE)
cc<-data.frame(do.call("rbind",gg),stringsAsFactors=FALSE)
reducedFunseqOutputNCDS$score<-as.numeric(cc[,1])

####

mergeDF<-reducedFunseqOutputNCDS

#fileName<-paste(tumorType,"_",mutationType,"_merge_variant_details.txt",sep="")
#fileName<-file.path(workDir,fileName)
#write.table(mergeDF,fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

groupType<-c("allSamples")

groupDF<-{}
groupDF[[groupType[1]]]<-mergeDF

#tmpDF<-reducedFunseqOutputNCDS

for(i in 1:length(groupType)){

  tmpDF<-groupDF[[groupType[i]]]

  groupName<-groupType[i]
  cat(sprintf("Processing %s\n",groupName))

  #tmpDF<-reducedFunseqOutputCDS

  posIndex<-paste(tmpDF$chr,tmpDF$posStart,tmpDF$posEnd,sep="@")
  geneDF<-data.frame(tmpDF$sampleID,posIndex,tmpDF$GENEparsed,tmpDF$ref,tmpDF$alt,tmpDF$score,stringsAsFactors = FALSE)
  colnames(geneDF)<-c("sampleID","posIndex","geneSymbol","ref","alt","score")

  if( nrow(geneDF[is.na(geneDF$score),])>0 ){
    geneDF[is.na(geneDF$score),]$score<-0
  }

  geneNameVector<-unique(geneDF$geneSymbol)
  geneDF<-split(geneDF,geneDF$geneSymbol)

  geneDFpatient<-{}

  for(geneName in geneNameVector){
    #tmpDat<-geneDF[[geneName]][!(duplicated(geneDF[[geneName]]$posIndex)),]
    tmpDat<-geneDF[[geneName]]
    npat<-length(unique(tmpDat$sampleID))
    geneDFpatient[[geneName]]<-npat
  }

  geneDFunique<-{}
  compositeScoreVector<-{}
  uniqueVariantPos<-{}

  for(geneName in geneNameVector){
    geneDFunique[[geneName]]<-geneDF[[geneName]][!(duplicated(geneDF[[geneName]]$posIndex)),]

    recurrenceVector<-table(geneDF[[geneName]]$posIndex)
    geneDFunique[[geneName]]$occurence<-recurrenceVector[as.character(geneDFunique[[geneName]]$posIndex)]
    geneDFunique[[geneName]]$compositeScore<-geneDFunique[[geneName]]$score*geneDFunique[[geneName]]$occurence
    compositeScoreVector[[geneName]]<-sum(geneDFunique[[geneName]]$compositeScore)
    uniqueVariantPos[[geneName]]<-nrow(geneDFunique[[geneName]])

  }

  geneDF<-rbind.fill(geneDF)

  geneDFunique<-rbind.fill(geneDFunique)
  compositeScoreDF<-data.frame(uniqueVariantPos,compositeScoreVector,stringsAsFactors = FALSE)
  rownames(compositeScoreDF)<-names(compositeScoreVector)
  colnames(compositeScoreDF)<-c("uniqueVariantPos","compositeScore")

######

  ncdsMutationFreq<-sort(table(compositeScoreDF$uniqueVariantPos),decreasing=TRUE)
  ncdsMutationCheckList<-rownames(compositeScoreDF)

######
#set.seed(42)
#mutationType<-"CDS"
#reSampleNum<-1000000

compositeFunseqScore<-{}
compositeFunseqScoreResample<-{}
numOfAlterationPos<-{}
numOfAlteration<-{}
numOfPatient<-{}
numOfAboveCFscore<-{}
pValue<-{}

outputDf<-{}

reSampleDistributionSize<-sort(as.numeric(names(ncdsMutationFreq)),decreasing = TRUE)
compositeFunseqScoreResample<-matrix(rep(0,length(reSampleDistributionSize)*reSampleNum),nrow=length(reSampleDistributionSize),ncol=reSampleNum)

cat(sprintf("Generate reSampling distritubtion\n"))

for(i in 1:length(reSampleDistributionSize)){
  #i<-21
  cat(sprintf("%s/%s reSampling distribution\n",i,length(reSampleDistributionSize)))

  if(reSampleDistributionSize[i]==1 && nrow(geneDFunique) < reSampleNum){
    tmpVector<-c(geneDFunique$compositeScore,rep(0,(reSampleNum-nrow(geneDFunique))))
    compositeFunseqScoreResample[i,]<-tmpVector

  }else{
    for(j in 1:reSampleNum){
      #i<-1
      index<-sample(1:nrow(geneDFunique),size=reSampleDistributionSize[i],replace=FALSE)
      compositeFunseqScoreResample[i,j]<-sum(geneDFunique[index,]$compositeScore)
      #compositeFunseqScoreResample[1]

    }


  }


}

#####

#tumorType<-"PRAD"
#dim(compositeFunseqScoreResample)

outputDir<-file.path(workDir,"result",mutationType)

if( !file.exists(outputDir) ){
  dir.create(outputDir,recursive=TRUE)
}

fileName<-paste(tumorType,"_",mutationType,"_",groupName,"_compositeFunseqScoreResample_iter_",reSampleNum,".Rd",sep="")
fileName<-file.path(outputDir,fileName)
save(compositeFunseqScoreResample,file=fileName)

fileName<-paste(tumorType,"_",mutationType,"_",groupName,"_ncdsMutationFreqTable.txt",sep="")
fileName<-file.path(outputDir,fileName)
write.table(geneDFunique,file=fileName,sep="\t",quote=FALSE,row.names =FALSE,col.names = TRUE)

reSampleSize<-{}
numOfgeneCheck<-length(unique(geneDFunique$geneSymbol))

cat(sprintf("Start calculating p-value\n"))

for(k in 1:numOfgeneCheck){
#for(k in 1:10){

  cat(sprintf("%s/%s\t",k,numOfgeneCheck))
  cat(sprintf("type:%s\tgene:%s\n",mutationType,ncdsMutationCheckList[k]))

  numOfAlterationPos[k]<-compositeScoreDF[rownames(compositeScoreDF) %in% ncdsMutationCheckList[k],]$uniqueVariantPos
  numOfAlteration[k]<-sum(geneDFunique[geneDFunique$geneSymbol %in% ncdsMutationCheckList[k],]$occurence)
  numOfPatient[k]<-geneDFpatient[[ncdsMutationCheckList[k]]]
  compositeFunseqScore[k]<-compositeScoreDF[rownames(compositeScoreDF) %in% ncdsMutationCheckList[k],]$compositeScore
  numOfAboveCFscore[k]<-sum(compositeFunseqScoreResample[which(reSampleDistributionSize==numOfAlterationPos[k]),]>=compositeFunseqScore[k])

  if(numOfAlterationPos[k]==1 && nrow(geneDFunique) < reSampleNum){
     if(compositeFunseqScore[k]==0){
        pValue[k]<-1
        numOfAboveCFscore[k]<-nrow(geneDFunique)
     }else{
        pValue[k]<-(numOfAboveCFscore[k]+1)/(nrow(geneDFunique)+1)
     }
    reSampleSize[k]<-nrow(geneDFunique)
  }else{
    pValue[k]<-(numOfAboveCFscore[k]+1)/(reSampleNum+1)
    reSampleSize[k]<-reSampleNum
  }


}

cat(sprintf("Finish calculating p-value\n"))

#outputDf<-data.frame(ncdsMutationCheckList[1:10],numOfAlteration,compositeFunseqScore,numOfAboveCFscore,rep(reSampleNum,length(ncdsMutationCheckList[1:10])),pValue)
outputDf<-data.frame(ncdsMutationCheckList,numOfAlterationPos,numOfAlteration,numOfPatient,compositeFunseqScore,numOfAboveCFscore,reSampleSize,pValue)
outputDf<-outputDf[order(outputDf$pValue),]
outputDf$qValue<-p.adjust(outputDf$pValue,method = "BH")
colnames(outputDf)<-c("geneSymbol","numOfAlterationPos","numOfAlteration","numOfPatient","fsigScore","numOfAbovefsigScore","reSampleNum","pValue","qValue")

fileName<-paste(tumorType,"_outputDf_",mutationType,"_",groupName,"_",reSampleNum,".txt",sep="")
fileName<-file.path(outputDir,fileName)
write.table(outputDf,file=fileName,sep="\t",quote=FALSE,row.names =FALSE,col.names = TRUE)

#removeIdx<-outputDf$numOfAlteration<3
#bb<-outputDf[!removeIdx,]
#bb$qValue<-p.adjust(bb$pValue,method="BH")
#bb<-bb[order(bb$pValue),]

}

return(outputDf)

}
