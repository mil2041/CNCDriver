#' Add trinucleotide distribution annotations
#'
#' @param dat reducedFunseqOutput data frame
#'
#'
#' @return triNucleotideDat data frame
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
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom SomaticSignatures mutationContext
#' @importFrom SomaticSignatures ucsc
#' @importFrom VariantAnnotation VRanges
#' @importFrom GenomicRanges GRanges





#library(plyr)
#library(GenomicRanges)
#library(data.table)
#library(SomaticSignatures)
#library("BSgenome.Hsapiens.UCSC.hg19")
#library("TxDb.Hsapiens.UCSC.hg19.knownGene")

addTriNucleotideDistribution<-function(dat){

  #dat<-reducedFunseqOutputNCDS
  
  txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  chromSize<-seqlengths(txdb)[1:24]

  #tmpDF<-reducedFunseqOutputWhole
  tmpDF<-dat
  #tmpDF<-reducedFunseqOutputCDS


  sca_gr = GRanges(
    seqnames = Rle(tmpDF$chr),
    ranges = IRanges(tmpDF$posStart,tmpDF$posEnd),
    ref = tmpDF$ref,
    alt = tmpDF$alt,
    sampleID=tmpDF$sampleID,
    study=tmpDF$tumorType,
    seqlengths = chromSize)

  sca_vr = VRanges(
    seqnames = seqnames(sca_gr),
    ranges = ranges(sca_gr),
    ref = mcols(sca_gr)$ref,
    alt = mcols(sca_gr)$alt,
    sampleNames = mcols(sca_gr)$sampleID,
    seqinfo = seqinfo(sca_gr),
    study = mcols(sca_gr)$study)

  # converting the seqnames notation to ’UCSC’ 
  sca_vr = ucsc(sca_vr)

  #sca_vr

  #sort(table(sca_vr$study), decreasing = TRUE)

  sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19, unify = TRUE)
  #head(sca_motifs)

  tmpDat<-mcols(sca_motifs)
  alteration<-unlist(strsplit(toString(tmpDat$alteration),"\\,"))
  alteration<-gsub(" ","",alteration)
  context<-unlist(strsplit(toString(tmpDat$context),"\\,"))
  context<-gsub(" ","",context)
  tmpDat<-data.frame(alteration,context,stringsAsFactors = FALSE)
  tmpDat$index<-paste(tmpDat$alteration,tmpDat$context,sep="@")

  return(tmpDat)

}


