# CNCDriver

[![Travis-CI Build Status](https://travis-ci.org/mil2041/CNCDriver.svg?branch=master)](https://travis-ci.org/mil2041/CNCDriver)
[![Last commit ](https://img.shields.io/github/last-commit/mil2041/CNCDriver.svg)]()
[![Repo Size ](https://img.shields.io/github/repo-size/mil2041/CNCDriver.svg)]()
[![Repo Size ](https://img.shields.io/github/release/mil2041/CNCDriver.svg)]()
[![License ](https://img.shields.io/github/license/mil2041/CNCDriver.svg?style=flat-square)]()

### References

* Eric Minwei Liu, Alexander Martinez Fundichely, Bianca Jay Diaz, Boaz Aronson, Tawny Cuykendall, Matthew MacKay, Priyanka Dhingra, Elissa WP Wong, Ping Chi, Effie Apostolou, Neville E Sanjana, Ekta Khurana [“Identification of Cancer Drivers at CTCF Insulators in 1,962 Whole Genomes”](https://www.ncbi.nlm.nih.gov/pubmed/31078526), Cell Systems, 8,4460455. e8 (2019)

### Overview

CNCDriver combined mutation recurrence and functional impact to identify coding and non-coding cancer drivers

### Version notes

* CNCDriver (version 0.3.3) add tutorial, Oct-10-2019
* CNCDriver (version 0.3.2) bugs fix
* CNCDriver (version 0.3.1) variable refacoring 
* CNCDriver (version 0.3) supports SNV coding drivers, promoter, enhancer, lincRNA and CTCF/cohesin insulator   




## Installation
User will need to install devtools in R for running CNCDriver package

``` r
library("remotes")
remotes::install_github("khuranalab/CNCDriver", ref="master", build_vignette=TRUE)
```

## Usage

``` r

library(CNCDriver)

#####
# global parameters setup
#####

  funseq2OutFile<-"path/to/funseq2/annotatedfile"

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


#####

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

#####

    mutationType<-"promoter_cluster"
    elementKeyWord<-"Promoter"
    
    promoterOutputDf<-getPromoterPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                        promoterRegionBedFile,elementKeyWord,
                                        minPoints,dRadius,
                                        triNucleotideDistributionFile,
                                        filterOutBlacklistMutations,
                                        mutationBlacklistFile,
                                        replicationTimingGenomeBinnedFile,
                                        replicationTimingElementBinnedFilePromoter,
                                        tumorType,mutationType,cellType,
                                        replicationTimingCutOff,
                                        seedNum,reSampleIterations,
                                        reRunPvalueCutOff,
                                        useCores,taskNum,unitSize,debugMode)
    
    makeQQplot(promoterOutputDf)


#####
    
    mutationType<-"lincRNA_cluster"
    elementKeyWord<-"lincRNA"
    
    lincRNAOutputDf<-getLincRNAPvalueWithPreFilter2(inputFileDir,outputFileDir,
                                      lincRNARegionBedFile,elementKeyWord,
                                      minPoints,dRadius,
                                      triNucleotideDistributionFile,
                                      filterOutBlacklistMutations,
                                      mutationBlacklistFile,
                                      replicationTimingGenomeBinnedFile,
                                      replicationTimingElementBinnedFileLincRNA,
                                      tumorType,mutationType,cellType,
                                      replicationTimingCutOff,
                                      seedNum,reSampleIterations,
                                      reRunPvalueCutOff,
                                      useCores,taskNum,unitSize,debugMode)
    
    
    makeQQplot(lincRNAOutputDf)

  

#####

    mutationType<-"enhancerUnit_cluster"
    elementKeyWord<-"Distal"
  
    enhancerUnitOutputDf<-getEnhancerUnitPvalueWithPreFilter2(inputFileDir,
                                                            outputFileDir,
                                              enhancerRegionBedFile,
                                              elementKeyWord,
                                              minPoints,dRadius,
                                              triNucleotideDistributionFile,
                                              filterOutBlacklistMutations,
                                              mutationBlacklistFile,
                                              replicationTimingGenomeBinnedFile,
                                              replicationTimingElementBinnedFileEnhancer,
                                              tumorType,mutationType,cellType,
                                              replicationTimingCutOff,
                                              seedNum,reSampleIterations,
                                              reRunPvalueCutOff,
                                              useCores,taskNum,unitSize,
                                              debugMode)
  
  
     makeQQplot(enhancerUnitOutputDf)
     
#####


```
## Contacts
For any questions, comments and suggestions, please email

* ekk2003 [at] med.cornell.edu 
* mil2041 [at] med.cornell.edu

Copyright © 2016-2019 Ekta Khurana Lab, WCMC

### License 
This project is licensed under the 
[![License ](https://img.shields.io/github/license/mil2041/CNCDriver.svg?style=flat-square)]()


