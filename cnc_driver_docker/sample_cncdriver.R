suppressMessages(library(CNCDriver))
print("CNCDriver loaded Successfully!")

#####
# change the path for the parsed Funseq results as input for CNCDriver
#####

#dataContextDir<-"/Users/nalomran/docker-image-funseq2-latest/data_context_funseq2"

# input files
inputFileDir<-"/home/parsed_FunSeq2_result_example/"

# output files
outputFileDir<-"/home/output/"

# required pre-processed additional information
cncdriver_data_Dir<-"/home/CNCDriver_data/"

###

triNucleotideDistributionFile<-file.path(cncdriver_data_Dir,"triNucleotidedistribution/compositeDriver_mutationDistributionMatrix_all_kmer_3_counts_Mar_2018_v0.txt")
mutationBlacklistFile<-file.path(cncdriver_data_Dir,"genome_blacklist/wgEncodeDacMapabilityConsensusExcludable.bed")
replicationTimingGenomeBinnedFile<-file.path(cncdriver_data_Dir,"replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows.txt")

###

codingRegionBedFile<-file.path(cncdriver_data_Dir,"FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.cds.bed")
proteinDomainFile<-file.path(cncdriver_data_Dir,"cds_data/domainData/pfam28.9606.processed_w_eValue_cutOff.txt")
proteinLengthFile<-file.path(cncdriver_data_Dir,"ensembl_hg37/hg37_ensembl_txCDSwithProteinLength.txt")
replicationTimingElementBinnedFileCDS<-file.path(cncdriver_data_Dir,"replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_cds.txt")

###

promoterRegionBedFile<-file.path(cncdriver_data_Dir,"FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.promoter.bed")
replicationTimingElementBinnedFilePromoter<-file.path(cncdriver_data_Dir,"replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_promoter.txt")

###

lincRNARegionBedFile<-file.path(cncdriver_data_Dir,"FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.lincRNA.bed")
replicationTimingElementBinnedFileLincRNA<-file.path(cncdriver_data_Dir,"replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_lincrna.txt")

###

enhancerRegionBedFile<-file.path(cncdriver_data_Dir,"FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/drm_FunSeq2.gene.bed")
replicationTimingElementBinnedFileEnhancer<-file.path(cncdriver_data_Dir,"replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_enhancer_May17_2018.txt")


###########
# load cancer information (name, cohort size, and replication time cell line)
###########

filePath<-"/home/"
fileName<-"cancer_type_info.txt"
fileName<-file.path(filePath,fileName)

cancer_info<-read.table(fileName,sep="\t",header=TRUE)

# only run GBM cancer in this example
cancer_info<-cancer_info[cancer_info$histologyName %in% "GBM",]

tumorTypeList<-cancer_info$histologyName
cellTypeList<-cancer_info$cellLineRTtype


#####
# Global parameters setup section
#####
# Default parameter settings
#####

replicationTimingCutOff<-0.2
filterOutBlacklistMutations<-TRUE

seedNum<-42
reSampleIterations<-10000
reRunPvalueCutOff<-0.1

minPoints<-2
dRadius<-50

useCores<-4
debugMode<-FALSE

taskNum<-0
unitSize<-100



######

for(i in 1:length(tumorTypeList)){
    
  tumorType<-tumorTypeList[i]
  cellType<-cellTypeList[i]

if(TRUE){
  
    mutationType<-"cds_cluster"
    elementKeyWord<-"CDS"
    cdsOutputDf<-getCDSPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                           codingRegionBedFile,elementKeyWord,
                                           proteinDomainFile,proteinLengthFile,
                                           minPoints,dRadius,
                                           triNucleotideDistributionFile,
                                           filterOutBlacklistMutations,
                                           mutationBlacklistFile,
                                           replicationTimingGenomeBinnedFile,
                                           replicationTimingElementBinnedFileCDS,
                                           tumorType,mutationType,cellType,
                                           replicationTimingCutOff,
                                           seedNum,reSampleIterations,
                                           reRunPvalueCutOff,
                                           useCores,taskNum,unitSize,debugMode)
    
    makeQQplot(cdsOutputDf)

}
}


