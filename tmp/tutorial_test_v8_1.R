library(CNCDriver)
source('~/work/Ekta_lab/codeRepository/CNCDriver/tmp/makeQQplot_v2.R')
#####
# global parameters setup
#####

dataContextDir<-"/path/to/dataContext"
annotatedInput<-"/path/to/Output.vcf"
outputDir<-"/path/to/output"
functionalImpactScoreCDS<-"FunSeq2"
functionalImpactScoreNonCoding<-"FunSeq2"

#tumorType<-"GBM"
#cellType<-"SKNSH"

#tumorType<-"CLL"
#cellType<-"K562"

#tumorType<-"BRCA"
#cellType<-"MCF7"

#tumorType<-"BLCA"
#cellType<-"average"


replicationTimingCutOff<-0.2
filterOutBlacklistMutations<-TRUE

seedNum<-42
reSampleIterations<-10000
reRunPvalueCutOff<-0.1


useCores<-4
debugMode<-TRUE

taskNum<-0
unitSize<-100

inputFileDir<-"/Users/mil20411/work/Ekta_lab/cncdriver_analysis_cluster_Mar_2018/compositeDriver_input"
outputFileDir<-"/Users/mil20411/work/Ekta_lab/Eric_et_al_Cell_Systems_2019_example/cncdriver_analysis_cluster_example_July_2020"
triNucleotideDistributionFile<-"~/work/Ekta_lab/compositeDriver_data/triNucleotidedistribution/compositeDriver_mutationDistributionMatrix_all_kmer_3_counts_Mar_2018_v0.txt"
mutationBlacklistFile<-"~/work/Ekta_lab/compositeDriver_data/genome_blacklist/wgEncodeDacMapabilityConsensusExcludable.bed"
replicationTimingGenomeBinnedFile<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows.txt"


codingRegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.cds.bed"
proteinDomainFile<-"~/work/Ekta_lab/compositeDriver_data/cds_data/domainData/pfam28.9606.processed_w_eValue_cutOff.txt"
proteinLengthFile<-"~/work/Ekta_lab/compositeDriver_data/ensembl_hg37/hg37_ensembl_txCDSwithProteinLength.txt"
replicationTimingElementBinnedFileCDS<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_cds.txt"


promoterRegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.promoter.bed"
replicationTimingElementBinnedFilePromoter<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_promoter.txt"


lincRNARegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.lincRNA.bed"
replicationTimingElementBinnedFileLincRNA<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_lincrna.txt"

enhancerRegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/drm_FunSeq2.gene.bed"
replicationTimingElementBinnedFileEnhancer<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_enhancer_May17_2018.txt"
#replicationTimingElementBinnedFileEnhancer<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_enhancer.txt"

###########

filePath<-file.path("~/work/Ekta_lab/cncdriver_analysis_Mar_2018/sample_summary")
fileName<-"Eric_sampleDist_Mar_2018.txt"
#fileName<-"Eric_sampleDist_OCT_2016_pancancer.txt"
fileName<-file.path(filePath,fileName)
tumorSampleSize<-read.table(fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE)

filePath<-file.path("~/work/Ekta_lab/cncdriver_analysis_Mar_2018/sample_summary")
fileName<-"tumorType_w_replicationTime_cellLine.txt"
#fileName<-"tumorType_w_replicationTime_cellLine_pancancer.txt"
fileName<-file.path(filePath,fileName)
cellLineMatch<-read.table(fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE)
colnames(cellLineMatch)<-c("histologyName","cellLineRTtype")

datConfig<-merge(tumorSampleSize,cellLineMatch,by="histologyName")
datConfig<-datConfig[order(datConfig$histologySize),]


minPoints<-2
dRadius<-50

datConfig<-datConfig[!(datConfig$histologyName %in% "Melanoma"),]


tumorTypeList<-datConfig$histologyName
cellTypeList<-datConfig$cellLineRTtype




for(i in 1:length(tumorTypeList)){
#for(i in 1:1){
  
  tumorType<-tumorTypeList[i]
  cellType<-cellTypeList[i]   

  #tumorType<-"GBM"
  #cellType<-"SKNSH"


#####

#preProcessVCF(annotatedInput,functionalImpactScoreCDS,functionalImpactScoreNonCoding,outputDir,tumorType,useCores)

#####

#mutationType<-"CDS"
#cdsOutputDf<-getCDSpvalue(outputDir,tumorType,mutationType,
#                          reSampleIter=reSampleIter,
#                          seedNum=seedNum,debugMode=debugMode)

  
  
if(FALSE){
  
mutationType<-"cds_cluster"
elementKeyWord<-"CDS"

#minPoints<-2
#dRadius<-50

cdsOutputDf<-getCDSPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                       codingRegionBedFile,elementKeyWord,
                                       proteinDomainFile,proteinLengthFile,
                                       minPoints,dRadius,
                                       triNucleotideDistributionFile,
                                       filterOutBlacklistMutations,mutationBlacklistFile,
                                       replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileCDS,
                                       tumorType,mutationType,cellType,replicationTimingCutOff,
                                       seedNum,reSampleIterations,reRunPvalueCutOff,
                                       useCores,taskNum,unitSize,debugMode)

}

if(FALSE){
  
mutationType<-"promoter_cluster"
elementKeyWord<-"Promoter"

promoterOutputDf<-getPromoterPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                    promoterRegionBedFile,elementKeyWord,
                                    minPoints,dRadius,
                                    triNucleotideDistributionFile,
                                    filterOutBlacklistMutations,mutationBlacklistFile,
                                    replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFilePromoter,
                                    tumorType,mutationType,cellType,replicationTimingCutOff,
                                    seedNum,reSampleIterations,reRunPvalueCutOff,
                                    useCores,taskNum,unitSize,debugMode)

#makeQQplot(promoterOutputDf)

}

  
if(TRUE){
    
mutationType<-"lincRNA_cluster"
elementKeyWord<-"lincRNA"

lincRNAOutputDf<-getLincRNAPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                  lincRNARegionBedFile,elementKeyWord,
                                  minPoints,dRadius,
                                  triNucleotideDistributionFile,
                                  filterOutBlacklistMutations,mutationBlacklistFile,
                                  replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileLincRNA,
                                  tumorType,mutationType,cellType,replicationTimingCutOff,
                                  seedNum,reSampleIterations,reRunPvalueCutOff,
                                  useCores,taskNum,unitSize,debugMode)


#makeQQplot(lincRNAOutputDf)


}
  
#####

#####
if(FALSE){
mutationType<-"promoter"
elementKeyWord<-"Promoter"

promoterOutputDf<-getPromoterPvalue(inputFileDir,outputFileDir,
                                    promoterRegionBedFile,elementKeyWord,
                                    triNucleotideDistributionFile,
                                    filterOutBlacklistMutations,mutationBlacklistFile,
                                    replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFilePromoter,
                                    tumorType,mutationType,cellType,replicationTimingCutOff,
                                    seedNum,reSampleIterations,reRunPvalueCutOff,
                                    useCores,taskNum,unitSize,debugMode)








#####

mutationType<-"lincRNA"
elementKeyWord<-"lincRNA"

lincRNAOutputDf<-getLincRNAPvalue(inputFileDir,outputFileDir,
                                  lincRNARegionBedFile,elementKeyWord,
                                  triNucleotideDistributionFile,
                                  filterOutBlacklistMutations,mutationBlacklistFile,
                                  replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileLincRNA,
                                  tumorType,mutationType,cellType,replicationTimingCutOff,
                                  seedNum,reSampleIterations,reRunPvalueCutOff,
                                  useCores,taskNum,unitSize,debugMode)

}
#####

if(FALSE){
mutationType<-"enhancerUnit"
elementKeyWord<-"Distal"

enhancerUnitOutputDf<-getEnhancerUnitPvalue(inputFileDir,outputFileDir,
                                            enhancerRegionBedFile,elementKeyWord,
                                            triNucleotideDistributionFile,
                                            filterOutBlacklistMutations,mutationBlacklistFile,
                                            replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileEnhancer,
                                            tumorType,mutationType,cellType,replicationTimingCutOff,
                                            seedNum,reSampleIterations,reRunPvalueCutOff,
                                            useCores,taskNum,unitSize,debugMode)


makeQQplot(enhancerUnitOutputDf)

}
  
#####

if(FALSE){  
  
  mutationType<-"enhancerUnit_cluster"
  elementKeyWord<-"Distal"
  
  enhancerUnitOutputDf<-getEnhancerUnitPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                              enhancerRegionBedFile,elementKeyWord,
                                              minPoints,dRadius,
                                              triNucleotideDistributionFile,
                                              filterOutBlacklistMutations,mutationBlacklistFile,
                                              replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileEnhancer,
                                              tumorType,mutationType,cellType,replicationTimingCutOff,
                                              seedNum,reSampleIterations,reRunPvalueCutOff,
                                              useCores,taskNum,unitSize,debugMode)
  
  
  #makeQQplot(enhancerUnitOutputDf)


}



}

