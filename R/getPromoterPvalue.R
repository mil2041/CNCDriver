#' Run CNCDriver promoter p-value calculation
#'
#' @param inputFileDir The path for prepared R object file [for example,reducedFunseqOutputNCDS_GBM.Rd]
#' @param outputFileDir The path for output files
#' @param promoterRegionBedFile The defintion of promoter regions in bed file format
#' @param elementKeyWord Default is "Promoter", the keyword of promoter annotation in FunSeq2
#' @param triNucleotideDistributionFile Cancer type specific mutation counts in 96 trinucleotide category  
#' @param filterOutBlacklistMutations TRUE or FALSE
#' @param mutationBlacklistFile The file for mutation blacklist regions
#' @param replicationTimingGenomeBinnedFile Replication timing value at 1MB bin
#' @param replicationTimingElementBinnedFile Replication timing value of each element within 1MB bin
#' @param tumorType Study name
#' @param mutationType User provided mutated gene list
#' @param cellType Cell type of the study, [ BJ, GM12878, HeLaS3, HepG2, IMR90, K562, MCF7, SKNSH or average ]
#' @param replicationTimingCutOff Default is 0.2, a numeric value ranging from 0 to 1
#' @param seedNum User provided random number seed, default is 42
#' @param reSampleIterations User provided re-sapling iteration numbers
#' @param reRunPvalueCutOff Default is 0.1, a numeric value ranging from 0 to 1
#' @param useCores Default is 1, number of cpu to use
#' @param taskNum Default is 0, 0 means to run all elements
#' @param unitSize  Number of elements to run per task
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
#' @import doRNG
#' @importFrom plyr rbind.fill
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom data.table data.table
#' 
getPromoterPvalue<-function(inputFileDir,outputFileDir,
                            promoterRegionBedFile,elementKeyWord="Promoter",
                            triNucleotideDistributionFile,
                            filterOutBlacklistMutations,mutationBlacklistFile,
                            replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFile,
                            tumorType,mutationType,cellType,replicationTimingCutOff,
                            seedNum=42,reSampleIterations,reRunPvalueCutOff=0.1,
                            useCores=1,taskNum=0,unitSize,debugMode=FALSE){


    #library(data.table)
    #library(plyr)
    ######
    
    #inputDir<-"~/work/Ekta_lab/Priyanka_project/fsigsnv_results"
    #tumorType<-"Prostate"
    
    
    replicationTimingCutOff<-as.numeric(replicationTimingCutOff)
    replicationTimingParm<-paste("RT_",replicationTimingCutOff,sep="")
    
    reSampleIterations<-as.numeric(reSampleIterations)
    reRunPvalueCutOff<-as.numeric(reRunPvalueCutOff)
    
    useCores<-as.numeric(useCores)
    taskNum<-as.numeric(taskNum)
    unitSize<-as.numeric(unitSize)
    
    workDir<-file.path(outputFileDir,"output_result_triMut_match",mutationType,tumorType,replicationTimingParm)
    
    #mutationType<-"CDS"
    #seedNum<-42
    #reSampleNum<-1000
    reSampleNum<-reSampleIterations
    set.seed(seedNum)
    
    if( !file.exists(paste(workDir,sep="/")) ){
      dir.create(paste(workDir,sep=""),recursive=TRUE)
    }
    
    cat(sprintf("Start processing %s - %s\n",tumorType,mutationType))
    
    
    ########
    # Load pre-processed Rd file
    ########
    filePath<-file.path(inputFileDir)
    fileName<-paste("reducedFunseqOutputNCDS_",tumorType,".Rd",sep="")
    fileName<-file.path(filePath,fileName)
    cat(sprintf("Load pre-processed Rd file\n"))
    load(fileName)
    
    #####
    ## load promoter bed file
    #####
    
    cat(sprintf("Load elementBedfile\n"))
    #filePath<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode"
    #fileName<-"gencode.v19.promoter.bed"
    fileName<-promoterRegionBedFile
    
    elementBedfileName<-file.path(promoterRegionBedFile)
    #reducedFunseqOutput<-reducedFunseqOutputNCDS
    
    
    elementBedfile<-read.table(elementBedfileName,sep="\t",header=FALSE,stringsAsFactors = FALSE)
    colnames(elementBedfile)<-c("chr","posStart","posEnd","name")
    #colnames(dat)<-c("chr","posStart","posEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
    
    
    ## extract mutations overlap with pre-defined element regions
    reducedFunseqOutputNCDS<-extractElementMutations(elementBedfile,reducedFunseqOutputNCDS)
    
    ######
    # load triNucleotideDistribution file
    ######
    cat(sprintf("Load triNucleotideDistribution file\n"))
    #filePath<-"~/work/Ekta_lab/compositeDriver_data/triNucleotidedistribution"
    #fileName<-"compositeDriver_mutationDistributionMatrix_all_kmer_3_counts_Mar_2018_v0.txt"
    #fileName<-file.path(filePath,fileName)
    
    fileName<-triNucleotideDistributionFile
    
    mutationDistMatrix<-read.table(fileName,sep="\t",header=TRUE,check.names = FALSE,stringsAsFactors = FALSE)
    
    categoryName<-rownames(mutationDistMatrix)
    categoryName<-gsub(" ","@",categoryName)
    rownames(mutationDistMatrix)<-categoryName
    
    ######
    # This section is for pan-cancer calculation
    # It is no required for single cancer type calculatin
    ######
    if(FALSE){
        filePath<-file.path("~/work/Ekta_lab/cncdriver_analysis_Mar_2018/sample_summary")
        fileName<-"Eric_sampleDist_Mar_2018.txt"
        fileName<-file.path(filePath,fileName)
        tumorDistTable<-read.table(fileName,header=TRUE,sep="\t",stringsAsFactors=FALSE)
        #tumorDistTable<-tumorDistTable[order(-tumorDistTable$totalMutCounts),]
        
        mutationBurden<-sort(colSums(mutationDistMatrix),decreasing = TRUE)
        mutationBurdenDF<-data.frame(names(mutationBurden),mutationBurden,stringsAsFactors = FALSE)
        colnames(mutationBurdenDF)<-c("histologyName","totalMutationCounts")
        
        tumorDistTable<-merge(tumorDistTable,mutationBurdenDF,by="histologyName")
        
        tumorDistTable$mutationBurdenAverage<-tumorDistTable$totalMutationCounts/tumorDistTable$histologySize
        
        tumorDistTable<-tumorDistTable[order(-tumorDistTable$mutationBurdenAverage),]
        
        
        ######
        
        #tumorDistTable<-tumorDistTable[tumorDistTable$histologyName %in% tumorType,]
        
        #sampleSizeFactor<-tumorDistTable$histologySize
        #names(sampleSizeFactor)<-tumorDistTable$histologyName
        
        #sampleSizeFactor<-tumorDistTable$totalMutationCounts
        #names(sampleSizeFactor)<-tumorDistTable$histologyName
        
        #sampleSizeFactor<-tumorDistTable$mutationBurdenAverage
        #names(sampleSizeFactor)<-tumorDistTable$histologyName
        
        #sampleSizeFactor<-rep(1,nrow(tumorDistTable))
        #names(sampleSizeFactor)<-tumorDistTable$histologyName
        
        sampleSizeFactor<-(log10(tumorDistTable$totalMutationCounts/min(tumorDistTable$totalMutationCounts)))+1
        names(sampleSizeFactor)<-tumorDistTable$histologyName
        
        #countsSizeFactor<-tumorDistTable$histologySize
        #names(countsSizeFactor)<-tumorDistTable$histologyName
        
        countsSizeFactor<-tumorDistTable$totalMutationCounts
        names(countsSizeFactor)<-tumorDistTable$histologyName
        
        #countsSizeFactor<-tumorDistTable$mutationBurdenAverage
        #names(countsSizeFactor)<-tumorDistTable$histologyName
    }
    ######
    
    
    
    ######
    
    # normalize tri-nucleotide context per cancer type
    if(TRUE){
      
      tmp_mm<-mutationDistMatrix
      
      #tmp_mm<-tmp_mm[,c(1:2)]
      tumorTypeSum<-colSums(tmp_mm)
      #overallSum<-sum(tumorTypeSum)
      
      # normalization on columns for each tumorType
      tmp_mm_t<-lapply(1:ncol(tmp_mm), function(x){
        tmp_mm[,x]/tumorTypeSum[colnames(tmp_mm)[x]]
      })
      
      tmp_mm_t<-do.call(rbind,tmp_mm_t)
      
      tmp_mm2<-t(tmp_mm_t)
      rownames(tmp_mm2)<-rownames(tmp_mm)
      colnames(tmp_mm2)<-colnames(tmp_mm)
      
      tmp_mm<-tmp_mm2
      
      tmp_mm_normalized<-tmp_mm
      
      mutationDistMatrix<-tmp_mm_normalized
      
    }
    
    cat(sprintf("Finish processing triNucleotideDistribution matrix\n"))
    
    ######
    
    parseFunseqGeneField<-function(field,keyword,useCores=1){
      geneSymbol<-{}
      #a1<-strsplit(aa,",")
      a1<-strsplit(field,"\\)")
      #a1<-unlist(a1)
      
      tmpStr<-mclapply(1:length(a1), function(x){
        #cat(sprintf("gene:%s\t",i))  
        idx <- which(grepl(keyword, a1[[x]]))
        a1[[x]]<-gsub("\\,","",a1[[x]])
        t3 <- strsplit(a1[[x]][idx], "\\(")
        #cat(sprintf("%s\n",t3[1][1]))
        tmpStr<-t3[[1]][1]
        
      },mc.cores=useCores)
      
      geneSymbol<-unlist(tmpStr)
      
      return(geneSymbol)  
    }
    
    cat(sprintf("Start parsing GENE field annotations\n"))
    reducedFunseqOutputNCDS$GENEparsed<-parseFunseqGeneField(field=reducedFunseqOutputNCDS$GENE,keyword=elementKeyWord,useCores=1)
    
    tmpString<-strsplit(as.character(reducedFunseqOutputNCDS$NCDS),":",fixed=TRUE)
    tmpStringFrame<-data.frame(do.call("rbind",tmpString),stringsAsFactors=FALSE)
    reducedFunseqOutputNCDS$score<-as.numeric(tmpStringFrame[,1])
    
    
    #######
    
    #######
    # filter out mutations on blacklist
    ######
    
    #filePath<-"~/work/Ekta_lab/compositeDriver_data/genome_blacklist"
    #fileName<-"wgEncodeDacMapabilityConsensusExcludable.bed"
    #fileName<-file.path(filePath,fileName)
    
    fileName<-mutationBlacklistFile
    
    testDT<-data.table(reducedFunseqOutputNCDS)
    setkeyv(testDT,c("chr","posStart","posEnd"))
    
    dacDT<-fread(fileName)
    colnames(dacDT)<-c("chr","posStart","posEnd","annotation","score","strand")
    setkeyv(dacDT,c("chr","posStart","posEnd"))
    
    mutExcludable<-foverlaps(testDT,dacDT,nomatch=0)
    numOfBlacklistMutations<-nrow(mutExcludable)
    
    cat(sprintf("there are %s mutations on blacklist\n",numOfBlacklistMutations))
    
    removeIdx<-mutExcludable[,c("chr","i.posStart","i.posEnd"),with=FALSE]
    
    
    if(filterOutBlacklistMutations){
      woBlacklistDT<-testDT[ !list(removeIdx),]
      reducedFunseqOutputNCDS<-setDF(woBlacklistDT)
      remove(woBlacklistDT)
      
    }
    
    remove(dacDT)
    remove(mutExcludable)
    remove(testDT)
    remove(removeIdx)
    
    
    #####
    # sanity check
    if(FALSE){
      
      reducedFunseqOutput<-reducedFunseqOutputNCDS  
      
      reducedFunseqOutput$index<-paste(reducedFunseqOutput$sampleID,reducedFunseqOutput$chr,paste(reducedFunseqOutput$posStart,reducedFunseqOutput$posEnd,sep="-"),sep=":")
      
      duplicatedLine<-duplicated(reducedFunseqOutput$index)
      
      educedFunseqOutput<-reducedFunseqOutput[!duplicatedLine,]
      
      reducedFunseqOutput<-reducedFunseqOutput[,c(1:28)]
      
      numOfWhitelistMutations<-nrow(reducedFunseqOutput)
      cat(sprintf("there are %s mutations on the regions of element bed file\n",numOfWhitelistMutations))
      
    }
    
    #####
    # Add triNucleotideDistribution
    #####
    
    cat(sprintf("Add triNucleotideDistribution information\n"))
    
    triNucleotideDat<-addTriNucleotideDistribution(reducedFunseqOutputNCDS)
    
    reducedFunseqOutputNCDS<-cbind(reducedFunseqOutputNCDS,triNucleotideDat)
    
    tmpDat<-reducedFunseqOutputNCDS
    
    ####
    # Parse motifBR, motifG, and NCENC field
    ####
    
    cat(sprintf("Parse MOTIFBR, MOTIFG, and NECNC annotations\n"))
    
    motifBreakFrame<-parseMOTIFBR(tmpDat,useCores)
    motifBreakFrameSimple<-motifBreakFrame[,c(2,7,9,8)]
    #motifBreakFrameSimple[!is.na(motifBreakFrameSimple$motifBreakMotifTFName),]
    
    motifGainFrame<-parseMOTIFG(tmpDat,useCores)
    motifGainFrameSimple<-motifGainFrame[,c(1,6,8,7)]
    #motifGainFrameSimple[!is.na(motifGainFrameSimple$motifGainMotifTFName),]
    
    motifAnnotation<-cbind(motifBreakFrameSimple,motifGainFrameSimple)
    
    #tmpDatSimple<-tmpDat[,c(1:6,9,12,19,20,22)]
    
    ncencFrame<-parseNCENC(tmpDat,useCores)
    outDf<-cbind(ncencFrame,motifAnnotation)
    
    reducedFunseqOutputNCDS<-cbind(tmpDat,outDf)
    
    ####
    
    ####
    ## associate replication timing data with each variant
    ## note: there is no chromosome Y or chromosome M data in the replication timing
    ####
    #cellType<-"BJ"
    
    cat(sprintf("Add replication timing signal for cellType %s\n",cellType))
    
    filePath<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed"
    fileName<-"UW_RepliSeq_wavelet_1Mb_windows.txt"
    fileName<-file.path(filePath,fileName)
    
    fileName<-replicationTimingGenomeBinnedFile
    
    replicationTimingDF<-read.table(fileName,header=TRUE,sep="\t",stringsAsFactors=FALSE)
    
    selectedColumn<-which(colnames(replicationTimingDF) %in% cellType)
    replicationTimingDF<-replicationTimingDF[,c(1:4,selectedColumn)]
    colnames(replicationTimingDF)<-c("chr","start","end","name","signalValue")
    replicationTimingDF$index<-paste(replicationTimingDF$chr,paste(replicationTimingDF$start,replicationTimingDF$end,sep="-"),sep=":")
    
    reducedFunseqOutputNCDS<-addReplicationTimingSignal(reducedFunseqOutputNCDS,replicationTimingDF)
    reducedFunseqOutputNCDS<-reducedFunseqOutputNCDS[!is.na(reducedFunseqOutputNCDS$signalValue),]
    
    #####
    # Set element-wise replication timing signal boundary
    #####
    cat(sprintf("Set element-wise replication timing signal boundary at %s similarity [range: 0-1]\n",replicationTimingCutOff))
    
    #filePath<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed"
    #fileName<-"UW_RepliSeq_wavelet_1Mb_windows_promoter.txt"
    
    fileName<-replicationTimingElementBinnedFile
    
    elementReplicationTimingDF<-fread(fileName,header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE)
    
    selectedColumn<-which(colnames(elementReplicationTimingDF) %in% cellType)
    elementReplicationTimingDF<-elementReplicationTimingDF[,c(1,selectedColumn)]
    colnames(elementReplicationTimingDF)<-c("elementName","signalValue")
    
    #replicationTimingCutOff<-0.3
    
    
    if(FALSE){
      elementMat<-elementReplicationTimingDF
      featureMat<-data.frame(elementMat[,2:ncol(elementMat)])
      colnames(featureMat)<-c("replicationTiming")
      rownames(featureMat)<-elementMat[,1]
      
      rangeMat<-getNearestElement(featureMat,replicationTimingCutOff,useCores=useCores)
      rangeMat$elementName<-rownames(rangeMat)
      
      elementReplicationTimingDF<-merge(elementReplicationTimingDF,rangeMat,by="elementName")
    }
    
    if(TRUE){
      #signalConst<-elementReplicationTimingDF$signalValue
      #names(signalConst)<-elementReplicationTimingDF$elementName
      
      maxValue<-max(elementReplicationTimingDF$signalValue)
      minValue<-min(elementReplicationTimingDF$signalValue)
      interval<-0.5*(maxValue-minValue)*replicationTimingCutOff
      
      #signalUpper<-rep(maxValue,nrow(elementReplicationTimingDF))
      #signalLower<-rep(minValue,nrow(elementReplicationTimingDF))
      
      signalUpper<-elementReplicationTimingDF$signalValue+interval
      signalLower<-elementReplicationTimingDF$signalValue-interval
      
      rangeValue<-data.frame(elementReplicationTimingDF,signalLower,signalUpper,stringsAsFactors = FALSE)
      colnames(rangeValue)<-c("elementName","signalValue","signalLower","signalUpper")
      
      idx<-which(rangeValue$signalUpper>=maxValue)
      deltaUpper<-rangeValue[idx,]$signalUpper-maxValue
      rangeValue[idx,]$signalLower<-rangeValue[idx,]$signalLower-deltaUpper
      rangeValue[idx,]$signalUpper<-rep(maxValue,length(idx))
      
      idx<-which(rangeValue$signalLower<=minValue)
      deltaLower<-minValue-rangeValue[idx,]$signalLower
      rangeValue[idx,]$signalUpper<-rangeValue[idx,]$signalUpper+deltaLower
      rangeValue[idx,]$signalLower<-rep(minValue,length(idx))
      
      rangeValue<-rangeValue[,-which(colnames(rangeValue) %in% "signalValue")]
    }
    
    elementReplicationTimingDF<-merge(elementReplicationTimingDF,rangeValue,by="elementName")
    
    # debug for plotting
    if(FALSE){
      
      cc<-elementReplicationTimingDF
      
      dd<-cc[order(cc$signalValue),]
      dd$index<-seq(1,nrow(cc),by=1)
      
      #name<-"TERT"
      #name<-"HOXC5"
      name<-"POTEB2"  
      #name<-"PCDH15"
      
      ss<-dd[dd$elementName %in% name,]
      #nameList<-dd[dd$signalValue<=ss$dmax & dd$signalValue>=ss$dmin,]$elementName
      nameList<-dd[dd$signalValue<=ss$signalUpper & dd$signalValue>=ss$signalLower,]$elementName
      
      kk<-dd[dd$elementName %in% nameList,]
      
      plot(dd$index,dd$signalValue,col="black")
      points(kk$index,kk$signalValue,col="red")
      points(ss$index,ss$signalValue,col="green",pch=19)
      
      plot(geneDFunique$signalValue,geneDFunique$occurence)
      
      kk<-geneDFunique[geneDFunique$signalValue<=ss$signalUpper & geneDFunique$signalValue>=ss$signalLower,]
      points(kk$signalValue,kk$occurence,col="red")
      #points(ss$signalValue,ss$signalValue,col="green",pch=19)
      
      plot(geneDFunique$signalValue,geneDFunique$compositeScore)
      points(kk$signalValue,kk$compositeScore,col="red")
      
    }
    
    cat(sprintf("Finish processing replication Timing file \n"))
    
    
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
      
      posIndex<-paste(tmpDF$chr,paste(tmpDF$posStart,tmpDF$posEnd,sep="-"),sep=":")
      geneDF<-data.frame(tmpDF$sampleID,posIndex,tmpDF$chr,tmpDF$posStart,
                         tmpDF$posEnd,tmpDF$GENEparsed,tmpDF$ref,tmpDF$alt,
                         tmpDF$TFPname,tmpDF$TFMname,tmpDF$TFMfullName,
                         tmpDF$motifBreakMotifTFName,tmpDF$motifGainMotifTFName,
                         tmpDF$score,tmpDF$RECUR,tmpDF$index,tmpDF$signalValue,
                         tmpDF$tumorType,stringsAsFactors = FALSE)
      #                     tmpDF$tumorType,stringsAsFactors = FALSE)
      colnames(geneDF)<-c("sampleID","posIndex","chr","posStart",
                          "posEnd","geneSymbol","ref","alt",
                          "TFPname","TFMname","TFMfullName",
                          "motifBreakMotifTFName","motifGainMotifTFName",
                          "score","recur","triMutIndex","signalValue","tumorType")
      
      # debug
      # checkScoreDF<-geneDF[is.na(geneDF$score),]
      # dim(checkScoreDF)
      
      if( nrow(geneDF[is.na(geneDF$score),])>0 ){
        geneDF[is.na(geneDF$score),]$score<-0
      }
      
      geneNameVector<-unique(geneDF$geneSymbol)
      geneDF<-split(geneDF,geneDF$geneSymbol)
      
      geneDFpatient<-{}
      
      geneDFpatient<-mclapply(1:length(geneNameVector), function(x){
        #tmpDat<-geneDF[[geneName]][!(duplicated(geneDF[[geneName]]$posIndex)),]
        geneName<-geneNameVector[x]
        tmpDF<-geneDF[[geneName]]
        npat<-length(unique(tmpDF$sampleID))
        geneDFpatient[[geneName]]<-npat
      },mc.cores=useCores)
      
      geneDFpatient<-unlist(geneDFpatient)
      names(geneDFpatient)<-geneNameVector
      
      compositeScoreVector<-{}
      compositeScoreScaledVector<-{}
      uniqueVariantPos<-{}
      
      ####
      #filePath<-file.path("~/work/Ekta_lab/cncdriver_analysis_Mar_2018/compositeDriver_input",mutationType)
      filePath<-file.path(inputFileDir,mutationType)
      fileName<-paste("geneDFunique_",tumorType,"_",mutationType,".Rd",sep="")
      fileName<-file.path(filePath,fileName)  
      
      if(file.exists(fileName)){
        
        cat(sprintf("Load %s\n",fileName))
        load(fileName)
        
      }else{  
        
        cat(sprintf("Processing geneDFunique object for calculation\n"))
        
        geneDFunique<-{}
        
        ## this number is 1 for single cancer type calculation
        ## A correction factor only used in pancancer calculation
        if(TRUE){
         sampleSizeFactor<-1
         names(sampleSizeFactor)<-tumorType
        }
        
        geneDFunique<-mclapply(1:length(geneNameVector), function(x){
          #geneDFunique<-mclapply(1:6, function(x){
          #for(x in 1:length(geneNameVector)){
          #cat(sprintf("%s\n",x))  
          geneName<-geneNameVector[x]
          geneDFunique[[geneName]]<-geneDF[[geneName]][!(duplicated(geneDF[[geneName]]$posIndex)),]
          geneDFunique[[geneName]]<-geneDFunique[[geneName]][order(geneDFunique[[geneName]]$posIndex),]
          recurrenceVector<-table(geneDF[[geneName]]$posIndex)
        
          collapsedDF<-mergeVariantPosition(geneDF[[geneName]],sampleSizeFactor,countsCutOff = 2)
          
          #geneDFunique[[geneName]]$occurence<-mergedDF[geneDFunique[[geneName]]$posIndex,]$occurence
          #geneDFunique[[geneName]]$compositeScore<-mergedDF[geneDFunique[[geneName]]$posIndex,]$compositeScore
          #geneDFunique[[geneName]]$tumorCounts<-mergedDF[geneDFunique[[geneName]]$posIndex,]$tumorCounts
          tmpDF<-geneDFunique[[geneName]]
          geneDFunique[[geneName]]<-cbind(tmpDF[,2:13],collapsedDF[,2:6])
          
          #}
          
        },mc.cores=useCores)
        
        names(geneDFunique)<-geneNameVector
        
        if( !file.exists(paste(filePath,sep="/")) ){
          dir.create(paste(filePath,sep=""),recursive=TRUE)
        }
        
        save(geneDFunique,file=fileName)
      }
      
      #####
      
      compositeScoreVector<-mclapply(1:length(geneNameVector), function(x){ 
        geneName<-geneNameVector[x]
        compositeScoreVector[[geneName]]<-sum(geneDFunique[[geneName]]$compositeScore)
      },mc.cores=useCores)
      
      compositeScoreVector<-unlist(compositeScoreVector)
      names(compositeScoreVector)<-geneNameVector
      
      #######
      
      compositeScoreScaledVector<-mclapply(1:length(geneNameVector), function(x){ 
        geneName<-geneNameVector[x]
        compositeScoreScaledVector[[geneName]]<-sum(geneDFunique[[geneName]]$compositeScoreScaled)
      },mc.cores=useCores)
      
      compositeScoreScaledVector<-unlist(compositeScoreScaledVector)
      names(compositeScoreScaledVector)<-geneNameVector
      
      #######
      
      uniqueVariantPos<-mclapply(1:length(geneNameVector), function(x){
        geneName<-geneNameVector[x]
        uniqueVariantPos[[geneName]]<-nrow(geneDFunique[[geneName]])
      },mc.cores=useCores)
      
      uniqueVariantPos<-unlist(uniqueVariantPos)
      names(uniqueVariantPos)<-geneNameVector
      
      #######
      
      numOfAlteration<-mclapply(1:length(geneNameVector), function(x){
        geneName<-geneNameVector[x]
        numOfAlteration<-sum(geneDFunique[[geneName]]$occurence)
      },mc.cores=useCores)
      
      numOfAlteration<-unlist(numOfAlteration)
      names(numOfAlteration)<-geneNameVector
      
      #######
      #geneDF<-rbind.fill(geneDF)
      geneDFunique<-rbind.fill(geneDFunique)
      
      #####
      
      #filePath<-file.path("~/work/Ekta_lab/cncdriver_analysis_Mar_2018/compositeDriver_input",mutationType)
      filePath<-file.path(inputFileDir,mutationType)
      fileName<-paste("variantTriMutCategoryParsed_",tumorType,"_",mutationType,".Rd",sep="")
      fileName<-file.path(filePath,fileName)  
      
      if(!file.exists(fileName)){
        
        cat(sprintf("Processing variantTriMutCategoryParsed object for calculation\n"))
        
        variantTriMutCategoryParsed<-parsePosVariantByTriMutContextWithAnnotation5(geneDFunique,mutationDistMatrix,useCores=useCores)
        
        if( !file.exists(paste(filePath,sep="/")) ){
          dir.create(paste(filePath,sep=""),recursive=TRUE)
        }
        
        save(variantTriMutCategoryParsed,file=fileName)
        
      }else{
        
        cat(sprintf("Load %s\n",fileName))
        load(fileName)
      }
      
      ####
      
      compositeScoreDF<-data.frame(uniqueVariantPos,numOfAlteration,geneDFpatient,compositeScoreVector,compositeScoreScaledVector,stringsAsFactors = FALSE)
      rownames(compositeScoreDF)<-names(compositeScoreVector)
      colnames(compositeScoreDF)<-c("uniqueVariantPos","numOfAlteration","numOfPatient","compositeScore","compositeScoreScaled")
      
      cat(sprintf("Finish data preparation\n"))
      
    
    ######
    #set.seed(42)
    #mutationType<-"CDS"
    #reSampleNum<-1000000
    
      compositeScore<-{}
      compositeScoreResample<-{}
      numOfAlterationPos<-{}
      numOfAlteration<-{}
      numOfPatient<-{}
      numOfAboveCScore<-{}
      reSampleSize<-{}
      pValue<-{}
      
      outputDf<-{}
    
    #######
    
      ncdsMutationCheckList<-unique(geneDFunique$geneSymbol)
      
      # remove some genes not have replication timing data
      ncdsMutationCheckList<-intersect(ncdsMutationCheckList,elementReplicationTimingDF$elementName)
      
      #ncdsMutationCheckList<-c("TERT","WDR74","PLEKHS1","LEPROTL1","TBC1D12","MYH14","TRPM4","CHODL","ZSCAN5B","ZNF595","XKR4","MMS19")
      #addList<-c("SDCCAG8","LRRC4C","ALG10","ALG10B","HIST1H2AJ","HDAC9","TTI2","EMC2","MGST3","PAN2","HIST1H2AM","ZNF143")
      
      #ncdsMutationCheckList<-c(ncdsMutationCheckList,addList)
      
      numOfgeneCheck<-length(ncdsMutationCheckList)
      
      cat(sprintf("%s elements have alterations\n",numOfgeneCheck))
      
      #numOfgeneCheck<-12
      
      #####
      # further split tasks into small units if large-scale parallelization calculations
      ####
      
      #unitSize<-12
      totalTaskNum<-ceiling(numOfgeneCheck/unitSize)
      #totalTaskNum
      
      loopVector<-{}
      
      for(i in 1:totalTaskNum){
        startNum<-1 + unitSize*(i-1)
        endNum<-startNum + unitSize -1
        if(endNum>numOfgeneCheck){
          endNum<-numOfgeneCheck
        }
        loopVector[[i]]<-c(startNum:endNum)
      }
      
      if(taskNum!=0){
        loopSequence<-loopVector[[taskNum]]
      }else{
        loopSequence<-c(1:numOfgeneCheck)
      }  
      
      
    #######  
    # main p-value calculation section
    #######
      
      detectCores()
      cl<-makeCluster(useCores)
      registerDoParallel(cl)
      
      
      time_start<-proc.time()
      cat(sprintf("Use %s cores\n",useCores))
      cat (sprintf ("Start calculating p-value for %s candidates\n", numOfgeneCheck) )
      
      filePath<-workDir
      fileName<-paste("log_",mutationType,"_triMut_distribution_RT_similarity_",replicationTimingCutOff,"_task_",taskNum,"_.txt",sep="")
      fileName<-file.path(filePath,fileName)
      
      cat(sprintf ("progress can be look up at \n"))
      cat(sprintf("%s\n",fileName))
      
      
      outputDf<-foreach(k=loopSequence, .options.RNG=seedNum, .packages=c("stringr","parallel","plyr","CNCDriver"))  %dorng% {
        #outputDf<-foreach(k=1:numOfgeneCheck, .options.RNG=seedNum, .packages=c("FNN","stringr","parallel","plyr"))  %dorng% {
        #for(k in 1:numOfgeneCheck){ 
        #filePath<-file.path("~/work/Ekta_lab/JasonWong_dataset/compositeFunSeq_result_triMut_match",tumorType)
        filePath<-workDir
        fileName<-paste("log_",mutationType,"_triMut_distribution_RT_similarity_",replicationTimingCutOff,"_task_",taskNum,"_.txt",sep="")
        fileName<-file.path(filePath,fileName)
        
        #cat(sprintf("%s/%s\t",k,numOfgeneCheck),file=fileName,append=TRUE)        
        cat(sprintf("%s/%s\ttype:%s\tgene:%s\n",k,numOfgeneCheck,mutationType,ncdsMutationCheckList[k]),file=fileName,append=TRUE)
        
        #cat(sprintf("%s/%s\t",k,numOfgeneCheck))        
        #cat(sprintf("type:%s\tgene:%s\n",mutationType,ncdsMutationCheckList[k]))
        
        numOfAlterationPos<-compositeScoreDF[rownames(compositeScoreDF) %in% ncdsMutationCheckList[k],]$uniqueVariantPos
        numOfAlteration<-sum(geneDFunique[geneDFunique$geneSymbol %in% ncdsMutationCheckList[k],]$occurence)
        numOfPatient<-geneDFpatient[[ncdsMutationCheckList[k]]]
        compositeScore<-compositeScoreDF[rownames(compositeScoreDF) %in% ncdsMutationCheckList[k],]$compositeScoreScaled
        
        gname<-ncdsMutationCheckList[k]
        #gname<-"TERT"
        #a1<-geneDFunique[geneDFunique$geneSymbol==gname,]
        a1<-variantTriMutCategoryParsed[variantTriMutCategoryParsed$geneSymbol==gname,]
        #reSamplePosSize<-nrow(a1)
        #compositeScore[k]<-sum(a1$compositeScore)
        
        
        #####
        #  still need to modifiy for panCancer triMutIndex complexity
        #####
        #variantTable<-{}
        #variantTable<-table(a1$triMutIndex)
        #variantTriMutCategory<-data.frame(names(variantTable),as.numeric(variantTable),stringsAsFactors = FALSE)
        
        #variantTriMutCategory<-makeTriMutCountsByPosition(a1$categoryCounts)
        #colnames(variantTriMutCategory)<-c("categoryName","variantCounts")
        
        #variantTriMutCategory<-parsePosVariantByTriMutContext(a1$categoryCounts)
        #variantTriMutCategory$compositeScore<-a1$compositeScore
        #variantTriMutCategory$posIndex<-a1$posIndex
        #variantTriMutCategory$geneSymbol<-a1$geneSymbol
        variantTriMutCategory<-a1
        
        #####
        
        #a2<-geneDFunique[geneDFunique$geneSymbol!=gname,]
        #a2<-a2[!is.na(a2$signalValue),]
        
        a2<-variantTriMutCategoryParsed[variantTriMutCategoryParsed$geneSymbol!=gname,]
        a2<-a2[!is.na(a2$signalValue),]
        
        #####
        featureDF<-elementReplicationTimingDF[elementReplicationTimingDF$elementName==gname,]
        #a3<-a2[a2$signalValue<=featureDF$signalUpper & a2$signalValue>=featureDF$signalLower,]
        #geneDFuniqueSelected<-a3
        
        backgroundPosVariant<-a2
        #####
        
        ## this step need reduce calculation to gain speedup
        #backgroundPosVariant<-parsePosVariantByTriMutContext(geneDFuniqueSelected$categoryCounts,useCores=1)
        
        #backgroundPosVariant$compositeScore<-geneDFuniqueSelected$compositeScore
        #backgroundPosVariant$posIndex<-geneDFuniqueSelected$posIndex
        #backgroundPosVariant$geneSymbol<-geneDFuniqueSelected$geneSymbol
        
        #backgroundPosVariant<-a3
        
        #####
        
        #filePath<-file.path("~/work/Ekta_lab/JasonWong_dataset/compositeFunSeq_result_triMut_match",tumorType)
        filePath<-workDir
        debugFileName<-paste("debug_",mutationType,"_triMut_distribution_RT_similarity_",replicationTimingCutOff,"_task_",taskNum,"_.txt",sep="")
        debugFileName<-file.path(filePath,debugFileName)
        
        ## plotting correlation for debug
        if(FALSE){ 
          filePath<-file.path("~/work/Ekta_lab/JasonWong_dataset/checkCorrelation",tumorType)
          fileName<-paste(tumorType,"_correlation_RT_",mutationType,"_",gname,".png",sep="")
          fileName<-file.path(filePath,fileName)
          png(fileName,width=800,height=600)
          plot(geneDFunique$signalValue,geneDFunique$compositeScore,pch=19,xlab="",ylab="",cex=0.6)
          title(main=paste(tumorType," (n=",sampleSize ,"), ",mutationType,sep=""),
                xlab="replication timing [ late to early ]",
                ylab="positional CompositeScore")
          mtext(paste(gname,"num pos:",numOfAlterationPos,sep=" "))
          points(geneDFuniqueSelected$signalValue,geneDFuniqueSelected$compositeScore,col="red",pch=19,cex=0.6)
          #points(a22$signalValue,a22$compositeScore,col="blue",pch=19,cex=0.6)
          points(a1$signalValue,a1$compositeScore,col="green",pch=19,cex=0.6)
          
          dev.off()
        }
        
        #replicationTimingCutOff<-0.1
        #system.time(
        result<-calculatePvalueWithMutCategoryDistributionMatched7(variantTriMutCategory,backgroundPosVariant,mutationDistMatrix,featureDF,reSampleNum,replaceFlag = FALSE,replicationTimingCutOff,debugFileName,debugMode)
        #)
        
        if(result$pValue<=reRunPvalueCutOff){
          
          reSampleNumExpanded<-reSampleNum*10
          result<-calculatePvalueWithMutCategoryDistributionMatched7(variantTriMutCategory,backgroundPosVariant,mutationDistMatrix,featureDF,reSampleNumExpanded,replaceFlag = FALSE,replicationTimingCutOff,debugFileName,debugMode)
          
        }
        
        
        pValue<-result$pValue 
        numOfAboveCScore<-result$numOfAboveCScore
        reSampleSize<-result$reSampleSize
        
        tmpResult<-data.frame(gname,numOfAlterationPos,numOfAlteration,numOfPatient,compositeScore,numOfAboveCScore,reSampleSize,pValue,stringsAsFactors = FALSE)
        #colnames(tmpResult)<-c("geneSymbol","numOfAlterationPos","numOfAlteration","numOfPatient","compositeDriverScore","numOfAboveCDscore","reSampleNum","pValue")
        colnames(tmpResult)<-c("elementPos","numOfAlterationPos","numOfAlteration","numOfPatient","compositeDriverScore","numOfAboveCDscore","reSampleNum","pValue")
        
        
        #filePath<-file.path("~/work/Ekta_lab/JasonWong_dataset/compositeFunSeq_result_triMut_match",tumorType)
        filePath<-workDir
        fileName<-paste("tmpOutput_",mutationType,"_triMut_distribution_RT_similarity_",replicationTimingCutOff,"_task_",taskNum,"_.txt",sep="")
        fileName<-file.path(filePath,fileName)
        
        write.table(tmpResult,file=fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)        
        
        
        
        return(tmpResult)
      }
          
      stopCluster(cl)
      time_elapesed<-proc.time()-time_start
      cat (sprintf ("Time for calculating p-value of %s candidates: %.2f sec\n", numOfgeneCheck, time_elapesed[3]) )
      
      outputDf<-rbind.fill(outputDf)
          
      
      #####
      
      #outputDf<-data.frame(ncdsMutationCheckList[1:numOfgeneCheck],numOfAlterationPos,numOfAlteration,numOfPatient,compositeScore,numOfAboveCScore,reSampleSize,pValue)
      outputDf<-outputDf[order(outputDf$pValue),]
      outputDf$qValue<-p.adjust(outputDf$pValue,method = "BH")
      #colnames(outputDf)<-c("geneSymbol","numOfAlterationPos","numOfAlteration","numOfPatient","compositeDriverScore","numOfAboveCDscore","reSampleNum","pValue","qValue")
      colnames(outputDf)<-c("elementPos","numOfAlterationPos","numOfAlteration","numOfPatient","compositeDriverScore","numOfAboveCDscore","reSampleNum","pValue","qValue")
      
      
      fileName<-paste(tumorType,"_outputDf_",mutationType,"_",groupName,"_with_RT_correction_similarity_",replicationTimingCutOff,"_triMut_match_distribution_task_",taskNum,"_iter_",reSampleNum,"_.txt",sep="")
      fileName<-file.path(workDir,fileName)  
      write.table(outputDf,file=fileName,sep="\t",quote=FALSE,row.names =FALSE,col.names = TRUE)
      
    
    }
    
    return(outputDf)

}
