#' Parse NCENC annotations from FunSeq2
#'
#' @param dat reducedFunseqOutput data frame
#' @param useCores default is one
#'
#' @return parsedNCENC data frame
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
#' @importFrom parallel mclapply
#' @importFrom plyr rbind.fill

parseNCENC<-function(dat,useCores=1){
  
  TFPname<-{}
  TFMname<-{}
  TFMfullName<-{}
  #a1<-strsplit(aa,",")
  
  a0<-dat$NCENC
  #posIndex<-paste(aa$chr,paste(aa$posStart,aa$posEnd,sep="-"),sep=":")
  
  #a1<-str_extract_all(a0,"([^,][A-Z0-9,_]+#[A-Za-z0-9_.:-]+#[0-9]+#[0-9]+#[+-]+#[0-9]+#[0-9.]+#[0-9.]+)")
  #a1<-unlist(a1)
  
  tmp<-mclapply(1:length(a0), function(x){
  #tmp<-lapply(1:length(a0), function(x){ 
    cat(sprintf("gene:%s / %s\n",x,length(a0)))  
    
    t3 <- strsplit(a0[[x]], "\\,")
    
    if(!is.na(t3[[1]][1])){
      
      tmpVector<-t3[[1]]
      tmpTFMvector<-tmpVector[grepl("TFM",tmpVector)]
      tmpTFPvector<-tmpVector[grepl("TFP",tmpVector)]
      
      #####
      if(length(tmpTFMvector)!=0){
        t4<-strsplit(tmpTFMvector,"\\|")
        t46<-data.frame(do.call("rbind",t4),stringsAsFactors=FALSE)
        
        TFMfullName<-paste(unique(t46[,2]),collapse=",")
        
        t47<-strsplit(t46[,2],"\\_")
        t48<-data.frame(do.call("rbind",t47),stringsAsFactors=FALSE)
        TFMname<-paste(unique(t48[,1]),collapse=",")
      
      }else{
        TFMname<-NA
        TFMfullName<-NA
      }
      #####
      if(length(tmpTFPvector)!=0){
        t5<-strsplit(tmpTFPvector,"\\(")
        t56<-data.frame(do.call("rbind",t5),stringsAsFactors=FALSE)
        t57<-strsplit(t56[,2],"\\|")
        t58<-data.frame(do.call("rbind",t57),stringsAsFactors=FALSE)
        TFPname<-paste(unique(t58[,1]),collapse=",")
      }else{
        TFPname<-NA
      }
      
      
      
    }else{
      
      TFPname<-NA
      TFMname<-NA
      TFMfullName<-NA
      
    }
    
    tmp<-data.frame(TFPname,TFMname,TFMfullName,
                    stringsAsFactors = FALSE)
    
    
    return(tmp)
  #})  
  },mc.cores=useCores)
  
  #motifBreakParsed<-list(fullName=motifBreakFullName,familyName=motifBreakFamilyName,strand=motifBreakStrand)
  parsedNCENC<-rbind.fill(tmp)
  return(parsedNCENC)  
}