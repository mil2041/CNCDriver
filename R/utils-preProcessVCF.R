#' Parse Funseq2 VCF file
#'
#' @param funseq2OutputFile Where is annotated funseq2 VCF file
#' @param outputDir path to store processed Rd file
#' @param tumorType study name
#' @param useCores number of cores to use for parallelization
#' @param debugMode TRUE or FALSE
#'
#' @return flag good or bad
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
#' @import VariantAnnotation
#' @importFrom utils write.table
#' @importFrom parallel mclapply
preProcessVCF<-function(funseq2OutputFile,outputDir,tumorType,useCores=6,debugMode=FALSE){

  #filePath<-file.path("~/work/Ekta_lab/PCAWG_AUG_2016/FunSeq2_annotation",chromosome)
  #fileName<-paste("chr",chrChar,".vcf.gz",sep="")
  #fileName<-"Output.vcf"
  #fileName<-file.path(filePath,fileName)
  #workDir<-"~/work/Ekta_lab/JasonWong_dataset/input_reduced"
  #workDir<-file.path("~/work/Ekta_lab/PCAWG_AUG_2016/FunSeq2_annotation",chromosome)

  cat(sprintf("Start parsing VCF file\n"))

  vcf<-readVcf(funseq2OutputFile,"hg19")
  #vcfTmp<-read.table(paste(filePath,fileName,sep="/"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  #vcf

  vcfDF<-info(vcf)

  #head(vcf)
  #Number Type   Description
  #SAMPLE  .      String Sample id
  #CDS     .      String Coding Variants or not
  #VA      .      String Coding Variant Annotation
  #HUB     .      String Network Hubs, PPI (protein protein interaction network), REG (regulatory network), PHOS (phosp...
  #GNEG    .      String Gene Under Negative Selection
  #GERP    .      String Gerp Score
  #NCENC   .      String NonCoding ENCODE Annotation
  #HOT     .      String Highly Occupied Target Region
  #MOTIFBR .      String Motif Breaking
  #MOTIFG  .      String Motif Gain
  #SEN     .      String In Sensitive Region
  #USEN    .      String In Ultra-Sensitive Region
  #UCONS   .      String In Ultra-Conserved Region
  #GENE    .      String Target Gene (For coding - directly affected genes ; For non-coding - promoter, enhancer regula...
  #CDSS    .      String Coding Score
  #NCDS    .      String NonCoding Score
  #RECUR   .      String Recurrent elements / variants
  #DBRECUR .      String Recurrence database

  keepVector<-as.logical(c(1,1,1,0,0,0,1,0,1,1,0,0,0,1,1,1,1,1))
  reducedVCF<-vcfDF[,keepVector]
  values(rowRanges(vcf))<-cbind(values(rowRanges(vcf)),reducedVCF)

  chr<-as.character(seqnames(vcf))
  posStart<-start(vcf)
  posEnd<-end(vcf)
  sampleID<-mcols(vcf)$SAMPLE
  sampleID[lapply(sampleID,length)==0]<-NA
  sampleID<-unlist(sampleID)

  ref<-as.character(ref(vcf))
  alt<-as.character(unlist(alt(vcf)))

  CDS<-mcols(vcf)$CDS
  CDS[lapply(CDS,length)==0]<-NA
  CDS<-unlist(CDS)

  VA<-mcols(vcf)$VA
  VA[lapply(VA,length)==0]<-NA
  VA <- mclapply(VA, function(x) {
    tmp <- paste(x, collapse=",")
    return(tmp)
  },mc.cores=useCores)

  VA<-unlist(VA)

  GENE<-mcols(vcf)$GENE
  GENE[lapply(GENE,length)==0]<-NA

  GENE <- mclapply(GENE, function(x) {
    tmp <- paste(x, collapse=",")
    return(tmp)
  },mc.cores=useCores)

  GENE<-unlist(GENE)

  CDSS<-mcols(vcf)$CDSS
  CDSS[lapply(CDSS,length)==0]<-NA
  CDSS<-unlist(CDSS)

  NCDS<-mcols(vcf)$NCDS
  NCDS[lapply(NCDS,length)==0]<-NA
  NCDS<-unlist(NCDS)

  NCENC<-mcols(vcf)$NCENC
  NCENC[lapply(NCENC,length)==0]<-NA


  #NCENC <- lapply(NCENC, function(x) {
  #  tmp <- paste(x, collapse=",")
  #  return(tmp)
  #})

  NCENC <- mclapply(NCENC, function(x) {
    tmp <- paste(x, collapse=",")
    return(tmp)
  },mc.cores=useCores)

  NCENC<-unlist(NCENC)


  MOTIFBR<-mcols(vcf)$MOTIFBR
  MOTIFBR[lapply(MOTIFBR,length)==0]<-NA
  MOTIFBR <- mclapply(MOTIFBR, function(x) {
    tmp <- paste(x, collapse=",")
    return(tmp)
  },mc.cores=useCores)

  MOTIFBR<-unlist(MOTIFBR)

  MOTIFG<-mcols(vcf)$MOTIFG
  MOTIFG[lapply(MOTIFG,length)==0]<-NA
  MOTIFG <- mclapply(MOTIFG, function(x) {
    tmp <- paste(x, collapse=",")
    return(tmp)
  },mc.cores=useCores)

  MOTIFG<-unlist(MOTIFG)

  RECUR<-mcols(vcf)$RECUR
  RECUR[lapply(RECUR,length)!=0]<-"Yes"
  RECUR[lapply(RECUR,length)==0]<-"No"
  RECUR<-unlist(RECUR)

  DBRECUR<-mcols(vcf)$DBRECUR
  DBRECUR[lapply(DBRECUR,length)!=0]<-"Yes"
  DBRECUR[lapply(DBRECUR,length)==0]<-"No"
  DBRECUR<-unlist(DBRECUR)

#####

  reducedDf<-data.frame(chr,posStart,posEnd,ref,alt,sampleID,CDS,VA,GENE,CDSS,NCDS,NCENC,MOTIFBR,MOTIFG,RECUR,DBRECUR,stringsAsFactors = FALSE)
  reducedDf<-reducedDf[order(reducedDf$chr,reducedDf$posStart,reducedDf$posEnd),]

  #reducedDfCombined<-merge(reducedDf,sampleTable,by="sampleID")
  reducedDf$tumorType<-rep(tumorType,nrow(reducedDf))

  cat(sprintf("Finish parsing VCF file\n"))

  outDir<-file.path(outputDir,tumorType,"input")

  if( !file.exists(paste(outDir,sep="/")) ){
    dir.create(paste(outDir,sep=""),recursive=TRUE)
  }


  if(debugMode){

  #workDir<-file.path("~/work/Ekta_lab/CDS_debug",tumorType)
  #fileName<-paste(chr[1],"_output_reduced.txt.gz",sep="")
  fileName<-paste(tumorType,"_output_reduced.txt.gz",sep="")
  outDir<-file.path(outputDir,tumorType,"input")

  fileName<-file.path(outDir,fileName)

  gz1 <- gzfile(fileName, "w")
  #write.table(reducedDf,file=paste(workDir,"/","chr",chr[1],"_output_reduced.txt",sep=""),quote=FALSE,row.names =FALSE,col.names = TRUE)
  write.table(reducedDf,gz1,sep="\t",quote=FALSE,row.names =FALSE,col.names = TRUE)
  close(gz1)

  }

  if(TRUE){
    reducedFunseqOutput<-reducedDf

    reducedFunseqOutputCDS<-reducedFunseqOutput[reducedFunseqOutput$CDS=="Yes",]
    reducedFunseqOutputNCDS<-reducedFunseqOutput[reducedFunseqOutput$CDS=="No",]

    fileName<-paste("reducedFunseqOutputCDS_",tumorType,".Rd",sep="")
    fileName<-file.path(outDir,fileName)
    save(reducedFunseqOutputCDS,file=fileName)

    fileName<-paste("reducedFunseqOutputNCDS_",tumorType,".Rd",sep="")
    fileName<-file.path(outDir,fileName)
    save(reducedFunseqOutputNCDS,file=fileName)
  }

  flag<-"Finish"

  return(flag)

}
