# CNCDriver

[![Travis-CI Build Status](https://travis-ci.org/mil2041/CNCDriver.svg?branch=master)](https://travis-ci.org/mil2041/CNCDriver)
[![Last commit ](https://img.shields.io/github/last-commit/mil2041/CNCDriver.svg)]()
[![Repo Size ](https://img.shields.io/github/repo-size/mil2041/CNCDriver.svg)]()
[![Repo Size ](https://img.shields.io/github/release/mil2041/CNCDriver.svg)]()
[![License ](https://img.shields.io/github/license/mil2041/CNCDriver.svg?style=flat-square)]()

### Overview

CNCDriver combined mutation recurrence and functional impact to identify coding and non-coding cancer drivers

### Version notes

* CNCDriver (version 0.3.2) bugs fix
* CNCDriver (version 0.3.1) variable refacoring 
* CNCDriver (version 0.3) supports SNV coding drivers, promoter, enhancer, lincRNA and CTCF/cohesin insulator   




## Installation
User will need to install devtools in R for running CNCDriver package

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("khuranalab/CNCDriver")
```

## Usage

``` r

library(CNCDriver)

#####
# global parameters setup
#####

  dataContextDir<-"/path/to/dataContext"
  annotatedInput<-"/path/to/Output.vcf"
  outputDir<-"/path/to/output"
  tumorType<-"GBM"
  seedNum<-42
  functionalImpactScoreCDS<-"FunSeq2"
  functionalImpactScoreNonCoding<-"FunSeq2"
  reSampleIter<-1000
  useCores<-6
  debugMode<-FALSE

#####
  
  preProcessVCF(annotatedInput,functionalImpactScoreCDS,functionalImpactScoreNonCoding,outputDir,tumorType,useCores)

#####

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

#####


```
## Contacts
For any questions, comments and suggestions, please email

* ekk2003 [at] med.cornell.edu 
* mil2041 [at] med.cornell.edu

Copyright © 2016-2018 Ekta Khurana Lab, WCMC

### License 
This project is licensed under the 
[![License ](https://img.shields.io/github/license/mil2041/CNCDriver.svg?style=flat-square)]()


