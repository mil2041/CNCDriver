#' Parse MOTIFBR annotations from FunSeq2
#'
#' @param dat reducedFunseqOutput data frame
#' @param useCores default is one
#'
#' @return motifBreakParsed data frame
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
#' 
parseMOTIFBR<-function(dat,useCores=1){
  
  motifBreakMotifFullName<-{}
  motifBreakMotifTFName<-{}
  motifBreakStrand<-{}
  motifBreakParsed<-{}
  #a1<-strsplit(dat,",")
  
  a0<-dat$MOTIFBR
  posIndex<-paste(dat$chr,paste(dat$posStart,dat$posEnd,sep="-"),sep=":")
  
  a1<-str_extract_all(a0,"([^,][A-Z0-9,_]+#[A-Za-z0-9_.:-]+#[0-9]+#[0-9]+#[+-]+#[0-9]+#[0-9.]+#[0-9.]+)")
  #a1<-unlist(a1)
  
  tmp<-mclapply(1:length(a1), function(x){
    #cat(sprintf("gene:%s / %s\n",i,length(a1)))  
    
    t3 <- strsplit(a1[[x]], "\\#")
    
    if(!is.na(t3[[1]][1])){
      cc<-data.frame(do.call("rbind",t3),stringsAsFactors=FALSE)
      colnames(cc)<-c("ChIPseqPeak","motifFullName","motifStart","motifEnd","strand","motifPos","altProb","refProb")
      
      cc$diffScore<-as.numeric(cc$altProb)-as.numeric(cc$refProb)
      
      cc<-cc[order(cc$diffScore),]
      
      # parse motif name
      t4<-strsplit(cc$motifFullName,"\\_")
      
      cc2<-data.frame(do.call(rbind,t4),stringsAsFactors = FALSE)
      colnames(cc2)<-c("motifNameRoot","motifNameExtend")
      
      # merge motif name
      motifBreakChIPpeakTFName<-paste(as.character(unique(unlist(strsplit(cc$ChIPseqPeak,",")))),collapse=",")
      motifBreakMotifTFName<-paste(as.character(unique(cc2$motifNameRoot)),collapse=",")
      motifBreakMotifFullName<-paste(as.character(cc$motifFullName),collapse=",")
      motifBreakStrand<-paste(as.character(cc$strand),collapse=",")
      motifBreakScoreAlt<-paste(as.character(cc$altProb),collapse=",")
      motifBreakScoreRef<-paste(as.character(cc$refProb),collapse=",")
      motifBreakScoreDiff<-paste(cc$diffScore,collapse=",")
      
      motifBreakMotifFullDetail<-paste(paste(cc$motifFullName,"[",paste(posIndex[x],cc$diffScore,sep="#"),"]",sep=""),collapse=",")
      
      maxDiffScore<-min(cc$diffScore)
      selectedIdx<-(cc$diffScore %in% maxDiffScore)
      
      cc3<-cc[selectedIdx,]
      t5<-strsplit(cc3$motifFullName,"\\_")
      cc4<-data.frame(do.call(rbind,t5),stringsAsFactors = FALSE)
      colnames(cc4)<-c("motifNameRoot","motifNameExtend")
      
      motifBreakMotifFullDetailMax<-paste(paste(cc3$motifFullName,"[",paste(posIndex[x],cc3$diffScore,sep="#"),"]",sep=""),collapse=",")
      motifBreakMotifTFNameMax<-paste(as.character(unique(cc4$motifNameRoot)),collapse=",")
            
    }else{
      
      motifBreakChIPpeakTFName<-NA
      motifBreakMotifFullName<-NA
      motifBreakMotifTFName<-NA
      motifBreakStrand<-NA
      motifBreakScoreAlt<-NA
      motifBreakScoreRef<-NA
      motifBreakScoreDiff<-NA
      motifBreakMotifFullDetail<-NA
      motifBreakMotifFullDetailMax<-NA
      motifBreakMotifTFNameMax<-NA
      
    }
    
    tmp<-data.frame(motifBreakChIPpeakTFName,motifBreakMotifTFName,
                    motifBreakMotifFullName,motifBreakStrand,
                    motifBreakScoreAlt,motifBreakScoreRef,
                    motifBreakMotifFullDetail,motifBreakMotifFullDetailMax,
                    motifBreakMotifTFNameMax,
                    stringsAsFactors = FALSE)
    
    
    return(tmp)
    
  },mc.cores=useCores)
  
  #motifBreakParsed<-list(fullName=motifBreakFullName,familyName=motifBreakFamilyName,strand=motifBreakStrand)
  motifBreakParsed<-rbind.fill(tmp)
  return(motifBreakParsed)  
}


