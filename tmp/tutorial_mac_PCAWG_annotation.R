library(CNCDriver)

#####
# global parameters setup
#####

#dataContextDir<-"/path/to/dataContext"
#annotatedInput<-"/path/to/Output.vcf"
#outputDir<-"/path/to/output"
#functionalImpactScoreCDS<-"FunSeq2"
#functionalImpactScoreNonCoding<-"FunSeq2"



#tumorType<-options[1]
#cellType<-options[2]

#tumorType<-"GBM"
#cellType<-"SKNSH"

#tumorType<-"Medulloblastoma"
#cellType<-"SKNSH"

#tumorType<-"CLL"
#cellType<-"K562"

tumorType<-"LUAD"
cellType<-"IMR90"


replicationTimingCutOff<-0.2
filterOutBlacklistMutations<-TRUE

seedNum<-42
reSampleIterations<-10000
reRunPvalueCutOff<-0.1


useCores<-6
debugMode<-FALSE

taskNum<-0
unitSize<-100

inputFileDir<-"~/work/Ekta_lab/TCGA_LUAD_WGS/cncdriver_analysis_Mar_2018/compositeDriver_input"
outputFileDir<-"~/work/Ekta_lab/TCGA_LUAD_WGS/cncdriver_analysis_Mar_2018"
triNucleotideDistributionFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS/cncdriver_data_context_LUAD/triNucleotidedistribution/compositeDriver_mutationDistributionMatrix_all_kmer_3_counts_LUAD_WGS_Mar24_2018_v0.txt"
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


#enhancerRegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/drm_FunSeq2.gene.bed"
enhancerRegionBedFile<-"~/work/Ekta_lab/TCGA_LUAD_WGS/data_context_PCAWG_checked/enhancers.PCAWG_OCT_2016.bed"
replicationTimingElementBinnedFileEnhancer<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_enhancer.txt"

#####

#preProcessVCF(annotatedInput,functionalImpactScoreCDS,functionalImpactScoreNonCoding,outputDir,tumorType,useCores)

#####

mutationType<-"cds"
elementKeyWord<-"CDS"

cdsOutputDf<-getCDSPvalue(inputFileDir,outputFileDir,
                          codingRegionBedFile,elementKeyWord,
                          proteinDomainFile,proteinLengthFile,
                          triNucleotideDistributionFile,
                          filterOutBlacklistMutations,mutationBlacklistFile,
                          replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFileCDS,
                          tumorType,mutationType,cellType,replicationTimingCutOff,
                          seedNum,reSampleIterations,reRunPvalueCutOff,
                          useCores,taskNum,unitSize,debugMode)


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



