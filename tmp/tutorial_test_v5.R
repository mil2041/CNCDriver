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


useCores<-10
debugMode<-FALSE

taskNum<-0
unitSize<-100

inputFileDir<-"/home/mil2041/work/Ekta_lab/TCGA_LUAD_WGS_AWG/cncdriver_analysis_NYGC_April_2018/compositeDriver_input"
outputFileDir<-"/home/mil2041/work/Ekta_lab/TCGA_LUAD_WGS_AWG/cncdriver_analysis_NYGC_April_2018"
triNucleotideDistributionFile<-"~/work/Ekta_lab/compositeDriver_data/triNucleotidedistribution/compositeDriver_mutationDistributionMatrix_all_kmer_3_counts_LUAD_WGS_NYGC_Apr18_2018_v0.txt"
mutationBlacklistFile<-"~/work/Ekta_lab/compositeDriver_data/genome_blacklist/wgEncodeDacMapabilityConsensusExcludable.bed"
replicationTimingGenomeBinnedFile<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows.txt"


codingRegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.cds.bed"
proteinDomainFile<-"~/work/Ekta_lab/compositeDriver_data/cds_data/domainData/pfam28.9606.processed_w_eValue_cutOff.txt"
proteinLengthFile<-"~/work/Ekta_lab/compositeDriver_data/ensembl_hg37/hg37_ensembl_txCDSwithProteinLength.txt"
replicationTimingElementBinnedFileCDS<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_cds.txt"


promoterRegionBedFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS_AWG/FunSeq2/data_context_PCAWG_checked/gencode/gencode.v19.promoterCore.bed"
replicationTimingElementBinnedFilePromoter<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_promoterCore_PCAWG_OCT_2016.txt"


utr5primeRegionBedFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS_AWG/FunSeq2/data_context_PCAWG_checked/gencode/gencode.v19.5utr.PCAWG_OCT_2016.bed"
replicationTimingElementBinnedFileUTR5prime<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_5UTRs_PCAWG_OCT_2016.txt"


utr3primeRegionBedFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS_AWG/FunSeq2/data_context_PCAWG_checked/gencode/gencode.v19.3utr.PCAWG_OCT_2016.bed"
replicationTimingElementBinnedFileUTR3prime<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_3UTRs_PCAWG_OCT_2016.txt"


lincRNARegionBedFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS_AWG/FunSeq2/data_context_PCAWG_checked/gencode/gencode.v19.lincRNA.PCAWG_OCT_2016.bed"
replicationTimingElementBinnedFileLincRNA<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_lincRNA_ncrna_PCAWG_OCT_2016.txt"

enhancerRegionBedFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS_AWG/FunSeq2/data_context_PCAWG_checked/drm_FunSeq2.gene.bed"
replicationTimingElementBinnedFileEnhancer<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_enhancer.txt"

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


tumorTypeList<-datConfig$histologyName
cellTypeList<-datConfig$cellLineRTtype


#for(i in 1:length(tumorTypeList)){
#for(i in 20:20){
  
  #tumorType<-tumorTypeList[i]
  #cellType<-cellTypeList[i]   

  tumorType<-"LUAD"
  cellType<-"IMR90"


#####

#preProcessVCF(annotatedInput,functionalImpactScoreCDS,functionalImpactScoreNonCoding,outputDir,tumorType,useCores)

#####

  mutationType<-"cds_original"
  elementKeyWord<-"CDS"
  
  #minPoints<-2
  #dRadius<-50
  
  cdsOutputDf<-getCDSPvalue(inputFileDir,outputFileDir,
                                         codingRegionBedFile,elementKeyWord,
                                         proteinDomainFile,proteinLengthFile,
                                         triNucleotideDistributionFile,
                                         filterOutBlacklistMutations,mutationBlacklistFile,
                                         replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileCDS,
                                         tumorType,mutationType,cellType,replicationTimingCutOff,
                                         seedNum,reSampleIterations,reRunPvalueCutOff,
                                         useCores,taskNum,unitSize,debugMode)

  makeQQplot(cdsOutputDf)
  
  
if(FALSE){  
  
  #mutationType<-"cds_cluster_pt_2_distance_50"
  mutationType<-"cds_cluster_pt_5_distance_50"
  mutationType<-"cds_cluster_pt_2_distance_200"
  elementKeyWord<-"CDS"
  
  minPoints<-2
  dRadius<-200
  
  cdsOutputDf<-getCDSPvalueWithPreFilter(inputFileDir,outputFileDir,
                                         codingRegionBedFile,elementKeyWord,
                                         proteinDomainFile,proteinLengthFile,
                                         minPoints,dRadius,
                                         triNucleotideDistributionFile,
                                         filterOutBlacklistMutations,mutationBlacklistFile,
                                         replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileCDS,
                                         tumorType,mutationType,cellType,replicationTimingCutOff,
                                         seedNum,reSampleIterations,reRunPvalueCutOff,
                                         useCores,taskNum,unitSize,debugMode)
  
  makeQQplot(cdsOutputDf)

}
  
############  
  
  minPoints<-2
  dRadius<-50
  
  mutationType<-"promoter_cluster_pt_2_distance_50"
  elementKeyWord<-"PromoterCore"
  
  promoterOutputDf<-getPromoterPvalueWithPreFilter(inputFileDir,outputFileDir,
                                      promoterRegionBedFile,elementKeyWord,
                                      minPoints,dRadius,
                                      triNucleotideDistributionFile,
                                      filterOutBlacklistMutations,mutationBlacklistFile,
                                      replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFilePromoter,
                                      tumorType,mutationType,cellType,replicationTimingCutOff,
                                      seedNum,reSampleIterations,reRunPvalueCutOff,
                                      useCores,taskNum,unitSize,debugMode)
  
  makeQQplot(promoterOutputDf)


  
#####

  
  minPoints<-2
  dRadius<-50
  
  mutationType<-"utr5prime_cluster_pt_2_distance_50"
  elementKeyWord<-"UTR"
  
  
  utr5primeOutputDf<-get5utrPvalueWithPreFilter(inputFileDir,outputFileDir,
                                                    utr5primeRegionBedFile,elementKeyWord,
                                                     minPoints,dRadius,
                                                     triNucleotideDistributionFile,
                                                     filterOutBlacklistMutations,mutationBlacklistFile,
                                                     replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileUTR5prime,
                                                     tumorType,mutationType,cellType,replicationTimingCutOff,
                                                     seedNum,reSampleIterations,reRunPvalueCutOff,
                                                     useCores,taskNum,unitSize,debugMode)
  
  makeQQplot(utr5primeOutputDf)
  
  
  
  ############
  
  minPoints<-3
  dRadius<-50
  
  mutationType<-"utr3prime_cluster_pt_3_distance_50"
  elementKeyWord<-"UTR"
  
  
  utr3primeOutputDf<-get3utrPvalueWithPreFilter(inputFileDir,outputFileDir,
                                                    utr3primeRegionBedFile,elementKeyWord,
                                                    minPoints,dRadius,
                                                    triNucleotideDistributionFile,
                                                    filterOutBlacklistMutations,mutationBlacklistFile,
                                                    replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileUTR3prime,
                                                    tumorType,mutationType,cellType,replicationTimingCutOff,
                                                    seedNum,reSampleIterations,reRunPvalueCutOff,
                                                    useCores,taskNum,unitSize,debugMode)
  
  makeQQplot(utr3primeOutputDf)
  
  
  ############
  
  mutationType<-"utr3prime_original"
  elementKeyWord<-"UTR"
  
  
  utr3primeOutputDf<-get3utrPvalue(inputFileDir,outputFileDir,
                                      utr3primeRegionBedFile,elementKeyWord,
                                      triNucleotideDistributionFile,
                                      filterOutBlacklistMutations,mutationBlacklistFile,
                                      replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileUTR3prime,
                                      tumorType,mutationType,cellType,replicationTimingCutOff,
                                      seedNum,reSampleIterations,reRunPvalueCutOff,
                                      useCores,taskNum,unitSize,debugMode)
  
  makeQQplot(utr3primeOutputDf)
  
#####

  
mutationType<-"promoter_original"
elementKeyWord<-"PromoterCore"

# TODO: figure out reason
# Warning message:
#  In getPromoterPvalue(inputFileDir, outputFileDir, promoterRegionBedFile,  :
#                         NAs introduced by coercion
                       
promoterOutputDf<-getPromoterPvalue(inputFileDir,outputFileDir,
                                    promoterRegionBedFile,elementKeyWord,
                                    triNucleotideDistributionFile,
                                    filterOutBlacklistMutations,mutationBlacklistFile,
                                    replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFilePromoter,
                                    tumorType,mutationType,cellType,replicationTimingCutOff,
                                    seedNum,reSampleIterations,reRunPvalueCutOff,
                                    useCores,taskNum,unitSize,debugMode)

makeQQplot(promoterOutputDf)






#####

mutationType<-"lincRNA_original"
elementKeyWord<-"lincRNA"

lincRNAOutputDf<-getLincRNAPvalue(inputFileDir,outputFileDir,
                                  lincRNARegionBedFile,elementKeyWord,
                                  triNucleotideDistributionFile,
                                  filterOutBlacklistMutations,mutationBlacklistFile,
                                  replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileLincRNA,
                                  tumorType,mutationType,cellType,replicationTimingCutOff,
                                  seedNum,reSampleIterations,reRunPvalueCutOff,
                                  useCores,taskNum,unitSize,debugMode)

makeQQplot(lincRNAOutputDf)


#####

minPoints<-2
dRadius<-50

mutationType<-"lincRNA_cluster_pt_2_distance_50"
elementKeyWord<-"lincRNA"

lincRNAOutputDf<-getLincRNAPvalueWithPreFilter(inputFileDir,outputFileDir,
                                  lincRNARegionBedFile,elementKeyWord,
                                  minPoints,dRadius,
                                  triNucleotideDistributionFile,
                                  filterOutBlacklistMutations,mutationBlacklistFile,
                                  replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileLincRNA,
                                  tumorType,mutationType,cellType,replicationTimingCutOff,
                                  seedNum,reSampleIterations,reRunPvalueCutOff,
                                  useCores,taskNum,unitSize,debugMode)

makeQQplot(lincRNAOutputDf)


#####

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

#####



