#' calculate Pvalue With Mut CategoryDistribution Matched
#'
#' @param variantTriMutCategory variantTriMutCategory data frame
#' @param backgroundPosVariant backgroundPosVariant data frame
#' @param mutationDistMatrix mutationDistMatrix data frame
#' @param featureDF featureDF data frame
#' @param reSampleNum number of re-sampling iterations
#' @param replaceFlag default is FALSE, whether to replace or not in the sampling process
#' @param replicationTimingCutOff, Default is 0.2, a numeric value ranging from 0 to 1
#' @param debugFileName fileName to output debug information
#' @param debugMode default is FALSE, TRUE or FALSE
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
calculatePvalueWithMutCategoryDistributionMatched7<-function(variantTriMutCategory,
                                                             backgroundPosVariant,mutationDistMatrix,
                                                             featureDF,reSampleNum,replaceFlag=FALSE,
                                                             replicationTimingCutOff,
                                                             debugFileName,debugMode=FALSE){
  
  ###
  ## use compositeDriverScoreScaled
  ##
  
  gname<-variantTriMutCategory$geneSymbol[1]
  
  tumorTypeInElement<-unique(unlist(strsplit(variantTriMutCategory$tumorName,",")))
  
  upperSignalbound<-featureDF$signalUpper
  lowerSignalbound<-featureDF$signalLower
  
  for(iterNum in 0:10){
    
    # 7 comes from 70/10, 70 is max value of replication timing value, 10 as 10 portions
    relaxFactor<-7*iterNum
    
    upperSignalbound<-min(upperSignalbound+relaxFactor,100)
    lowerSignalbound<-max(lowerSignalbound-relaxFactor,0)
    
    a2<-backgroundPosVariant
    a3<-a2[a2$signalValue<=upperSignalbound & a2$signalValue>=lowerSignalbound,]
    backgroundPosVariantTmp<-a3
    
    drawSize<-nrow(variantTriMutCategory)
    
    #criteria1<-grepl(paste(unique(tumorTypeInElement),collapse="|"),backgroundPosVariantTmp$tumorName)
    
    ## where there is no avaialbe position 
    #backgroundPosVariantMatched<-backgroundPosVariantTmp[criteria1,]
    #cat(sprintf("backgroundPosVariantMatched: %s\n",nrow(backgroundPosVariantMatched)))
    #if(nrow(backgroundPosVariantMatched)>0){break}
    
    ## where there is no avaialbe position 
    backgroundPosVariantMatched<-backgroundPosVariantTmp
    cat(sprintf("backgroundPosVariantMatched: %s\n",nrow(backgroundPosVariantMatched)))
    if(nrow(backgroundPosVariantMatched)>0)
    {
      
      
      backgroundSize<-nrow(backgroundPosVariantMatched)
      
      if(backgroundSize>=drawSize){
        
        break
      }
      
    }
    
  }
  
  backgroundPosVariantMatched$index<-seq(1,nrow(backgroundPosVariantMatched),by=1)
  rownames(backgroundPosVariantMatched)<-backgroundPosVariantMatched$index
  
  #backgroundCandidate<-{}
  
  #backgroundPosVariant$probSelected<-makeTriMutProbAdjustment(backgroundPosVariant,mutationDistMatrix)
  
  #drawSize<-nrow(variantTriMutCategory)
  backgroundSize<-nrow(backgroundPosVariantMatched)
  
  combinationNumEstimate<-choose(backgroundSize,drawSize)
  
  
  ####
  #  start of random shuffle
  ####
  
  ## Legacy code
  if(FALSE){
    
        if(combinationNumEstimate < reSampleNum){
          
          # work on this for panCancer
          #nearestCategoryName<-getNearestCategory(mutCategoryName, mutFreqLocal,replicationTimingCutOff)
          #mutCategoryList[[i]]<-nearestCategoryName
          # do we have enough combinations after category matching expansion ?
          #expandedGroupCounts[i]<-sum(mutFreqLocal[mutFreqLocal$categoryName %in% mutCategoryList[[i]],]$counts)
          #expandedCombnNum[i]<-choose(expandedGroupCounts[[i]],mutCategoryCandidate$variantCounts)
          
          cat(sprintf("[ Warning ]: %s has variant need to extend matching space\n",gname),file=debugFileName,append=TRUE)
          
        }else{
          
          cat(sprintf("[ Good ]: element %s has combnNum: %s that is bigger than reSampleNum: %s\n",gname,combinationNumEstimate,reSampleNum),file=debugFileName,append=TRUE)
          
        }
  
  }  
  
  #####
  # geneDFuniqueMatched<-geneDFuniqueSelected[geneDFuniqueSelected$triMutIndex %in% mutCategorySelected,]
  #####
  
  numOfAlterationPos<-nrow(variantTriMutCategory)
  #compositeScore<-sum(variantTriMutCategory$compositeScore)  
  compositeScoreScaled<-sum(variantTriMutCategory$compositeScoreScaled)  
  #compositeScoreScaled<-sum(variantTriMutCategory$compositeScore)  
  
  if(numOfAlterationPos==1 && combinationNumEstimate < reSampleNum){
    
    ######
    ## We can calculate exact p-value without re-sampling
    ######
    
    if(compositeScoreScaled==0){
      
      #pValue<-1
      #numOfAboveCScore<-nrow(backgroundPosVariantMatched)
      compositeScoreResample<-backgroundPosVariantMatched$compositeScoreScaled
      numOfAboveCScore<-sum(compositeScoreResample>=compositeScoreScaled)
      
      pValue<-(numOfAboveCScore+1)/(nrow(backgroundPosVariantMatched)+1)
      
      
    }else{
      
      ######
      #
      # when numOfAlterationPos==1 && ( nrow(geneDFuniqueMatched) < reSampleNum ) && compositeScore !=0
      #
      ######
      
      #indexMatrix<-makeReSampleIndexMatrixParallel(nrow(geneDFuniqueSelected),numOfAlterationPos[k],reSampleNum,6,seedNum = 123)
      #compositeScoreResample<-{}
      #compositeScoreResample<-sapply(1:nrow(indexMatrix), function(x) sum(geneDFuniqueSelected[indexMatrix[x,],]$compositeScore))
      compositeScoreResample<-backgroundPosVariantMatched$compositeScoreScaled
      numOfAboveCScore<-sum(compositeScoreResample>=compositeScoreScaled)
      
      pValue<-(numOfAboveCScore+1)/(nrow(backgroundPosVariantMatched)+1)
      
    }
    
    reSampleSize<-nrow(backgroundPosVariantMatched)
    
  }else{
    
    ######
    ## We need to do re-sampling to calculate p-value
    ######
    
    time_start<-proc.time()
    #indexMatrix<-makeReSampleIndexMatrixParallel(nrow(geneDFuniqueSelected),numOfAlterationPos[k],reSampleNum,useNodes,seedNum = 123)
    #variantGroupSize<-nrow(geneDFuniqueMatched)
    #probSelected<-geneDFuniqueSelected$prob
    
    if(combinationNumEstimate < reSampleNum){
      
      drawIteration<-combinationNumEstimate
    
    }else{
      
      drawIteration<-reSampleNum
      
    }
    
    if(numOfAlterationPos>1){
      
      # each column is a re-sampleing index that repeats # of column times 
      #indexMatrix<-replicate(reSampleNum,sample(1:variantGroupSize,size=numOfAlterationPos,replace=FALSE),simplify = "matrix")
      
      indexMatrix<-replicate(drawIteration,sample(backgroundPosVariantMatched$compositeScoreScaled,size=numOfAlterationPos,replace=replaceFlag,prob=backgroundPosVariantMatched$probSelected),simplify = "matrix")
      compositeScoreResample<-colSums(indexMatrix)
      
      #reSampleResult<-makeReSampleIndex(variantTriMutCategory,backgroundPosVariantMatched,combinationNumEstimate,reSampleNum,debug)
      #indexMatrix<-reSampleResult$indexMatrix
      combinationNum<-drawIteration   
      
    }else{
      
      ######  
      # numOfAlterationPos is 1 and nrow(backgroundPosVariantMatched) >= reSampleNum
      ######
      
      #indexMatrix<-matrix(replicate(reSampleNum,sample(1:variantGroupSize,size=1,replace=FALSE),simplify = "matrix"),ncol=reSampleNum)
      ## next work from here
      #indexMatrix<-matrix(replicate(reSampleNum,sample(1:variantGroupSize,size=1,replace=FALSE,prob=probSelected),simplify = "matrix"),ncol=reSampleNum)
      
      indexMatrix<-replicate(drawIteration,sample(backgroundPosVariantMatched$compositeScoreScaled,size=numOfAlterationPos,replace=replaceFlag,prob=backgroundPosVariantMatched$probSelected),simplify = "matrix")
      compositeScoreResample<-indexMatrix
      
      combinationNum<-drawIteration 
    }
    
    #time_elapesed<-proc.time()-time_start
    #cat (sprintf ("Time for build reSampleing matrix: %.2f sec\n", time_elapesed[3]) )
    
    if(debugMode){
    
        if(combinationNum < reSampleNum){ 
          
          cat(sprintf("[ Warning ]: %s combnNum: %s is less than reSampleNum: %s\n",gname,combinationNum,reSampleNum),file=debugFileName,append=TRUE)
          #cat(sprintf("Warning: %s combnNum: %s is less than reSampleNum: %s\n",gname,combnNum,reSampleNum))
          
          # speacial awareness!!
          reSampleNum<-combinationNum
          
        }else{
          
          cat(sprintf(" [ Good ]: %s combnNum: %s is equal or bigger than reSampleNum: %s\n",gname,combinationNumEstimate,reSampleNum),file=debugFileName,append=TRUE)
          
        }
    
    }
    
    if(!debugMode){
      #compositeScoreResample<-{}
      
      # this step takes time
      #system.time(
      #compositeScoreResample<-sapply(1:ncol(indexMatrix), function(x) sum(backgroundPosVariantMatched[ backgroundPosVariantMatched$index %in% indexMatrix[,x],]$compositeScore))
      #)
      
      #compositeScoreResample<-reSampleResult$compositeScoreResample
      
      
      numOfAboveCScore<-sum(compositeScoreResample>=compositeScoreScaled)
      
      #hist(compositeScoreResample,breaks=100,col="grey",xlab="compositeScore",xlim=c(0,compositeScore+5))
      #abline(v=compositeScore,col="red")
      
      pValue<-(numOfAboveCScore+1)/(drawIteration+1)
      
    }else{
      
      
      
      numOfAboveCScore<-sum(compositeScoreResample>=compositeScoreScaled)
      pValue<-(numOfAboveCScore+1)/(drawIteration+1)
      
      #####
      ## In debugMode, skip p-value calculations
      #####
      
      #numOfAboveCScore<-0
      #pValue<-1
      
    }
    
    reSampleSize<-drawIteration
    
  }
  
  
  result<-list(pValue=pValue,numOfAboveCScore=numOfAboveCScore,reSampleSize=reSampleSize)
  return(result)
  
  
}