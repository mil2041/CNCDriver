library(CNCDriver)

#####
# global parameters setup
#####

dataContextDir<-"/path/to/dataContext"
annotatedInput<-"/path/to/Output.vcf"
outputDir<-"/path/to/output"
functionalImpactScoreCDS<-"FunSeq2"
functionalImpactScoreNonCoding<-"FunSeq2"

tumorType<-"GBM"
cellType<-"SKNSH"
replicationTimingCutOff<-0.2
filterOutBlacklistMutations<-TRUE

seedNum<-42
reSampleIterations<-10000
reRunPvalueCutOff<-0.1


useCores<-10
debugMode<-TRUE

taskNum<-0
unitSize<-100

inputFileDir<-"/home/mil2041/work/Ekta_lab/cncdriver_test/compositeDriver_input"
outputFileDir<-"/home/mil2041/work/Ekta_lab/cncdriver_test"
promoterRegionBedFile<-"~/work/Ekta_lab/FunSeq2_compositeDriver_SEP_2017/data_context_SEP_2017/gencode/gencode.v19.promoter.bed"
triNucleotideDistributionFile<-"~/work/Ekta_lab/compositeDriver_data/triNucleotidedistribution/compositeDriver_mutationDistributionMatrix_all_kmer_3_counts_Mar_2018_v0.txt"
mutationBlacklistFile<-"~/work/Ekta_lab/compositeDriver_data/genome_blacklist/wgEncodeDacMapabilityConsensusExcludable.bed"
replicationTimingGenomeBinnedFile<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows.txt"
replicationTimingElementBinnedFile<-"~/work/Ekta_lab/compositeDriver_data/replication_timing/UW/processed/UW_RepliSeq_wavelet_1Mb_windows_promoter.txt"


#####

preProcessVCF(annotatedInput,functionalImpactScoreCDS,functionalImpactScoreNonCoding,outputDir,tumorType,useCores)

#####

mutationType<-"CDS"
cdsOutputDf<-getCDSpvalue(outputDir,tumorType,mutationType,
                          reSampleIter=reSampleIter,
                          seedNum=seedNum,debugMode=debugMode)

#####

#####

mutationType<-"promoter"
elementKeyWord<-"Promoter"

promoterOutputDf<-getPromoterPvalue(inputFileDir,outputFileDir,
                                    promoterRegionBedFile,elementKeyWord,
                                    triNucleotideDistributionFile,
                                    filterOutBlacklistMutations,mutationBlacklistFile,
                                    replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFile,
                                    tumorType,mutationType,cellType,replicationTimingCutOff,
                                    seedNum,reSampleIterations,reRunPvalueCutOff,
                                    useCores,taskNum,unitSize,debugMode)




#####

mutationType<-"lincRNA"
elementKeyWord<-"lincRNA"

lincRNAOutputDf<-getLincRNAPvalue(inputFileDir,outputFileDir,
                                  promoterRegionBedFile,elementKeyWord,
                                  triNucleotideDistributionFile,
                                  filterOutBlacklistMutations,mutationBlacklistFile,
                                  replicationTimingGenomeBinnedFile,replicationTimingElementBinnedFile,
                                  tumorType,mutationType,cellType,replicationTimingCutOff,
                                  seedNum,reSampleIterations,reRunPvalueCutOff,
                                  useCores,taskNum,unitSize,debugMode)


#####

fileName<-"drm.gene.bed"
fileName<-file.path(dataContextDir,fileName)
enhancerGeneInteractionFileName<-fileName

mutationType<-"enhancerUnit"
enhancerUnitOutputDf<-getEnhancerUnitPvalue(outputDir,tumorType,mutationType,
                                            reSampleIter=reSampleIter,
                                            seedNum=seedNum,
                                            enhancerGeneInteractionFileName,
                                            useCores=useCores,
                                            debugMode=debugMode)

#####
