#' Parse MOTIFBR annotations from FunSeq2
#'
#' @param dat reducedFunseqOutput data frame
#' @param useCores default is one
#'
#' @return motifGainParsed data frame
#'
#' @examples
#' # TCF12_disc4[chr19:49990691-49990691#2.23] 
#' # motifmodel # mutation position # motif start # motif end # motif direction # mutation position on the motif model # delta score 
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


parseMOTIFG<-function(dat,useCores=1){
  motifGainMotifFullName<-{}
  motifGainMotifTFName<-{}
  motifGainStrand<-{}
  motifGainParsed<-{}
  #a1<-strsplit(dat,",")
  
  a0<-dat$MOTIFG
  posIndex<-paste(dat$chr,paste(dat$posStart,dat$posEnd,sep="-"),sep=":")
  
  
  a0[a0=="NA"]<-NA
  a1<-strsplit(a0,"\\,")
  #a1<-unlist(a1)
  
  tmp<-mclapply(1:length(a1), function(x){
    #cat(sprintf("gene:%s / %s\n",i,length(a1))) 
    
    if(!is.na(a1[[x]][1])){
      
      t3 <- strsplit(a1[[x]], "\\#")
      cc<-data.frame(do.call("rbind",t3),stringsAsFactors=FALSE)
      colnames(cc)<-c("motifFullName","motifStart","motifEnd","strand","motifPos","altScore","refScore")
      
      cc$diffScore<-as.numeric(cc$altScore)-as.numeric(cc$refScore)
      
      # gain, large score is better 
      cc<-cc[order(-cc$diffScore),]
      
      # parse motif TF name
      t4<-strsplit(cc$motifFullName,"\\_")
      
      cc2<-data.frame(do.call(rbind,t4),stringsAsFactors = FALSE)
      colnames(cc2)<-c("motifNameRoot","motifNameExtend")
      
      # merge motif name
      motifGainMotifTFName<-paste(as.character(unique(cc2$motifNameRoot)),collapse=",")
      motifGainMotifFullName<-paste(as.character(cc$motifFullName),collapse=",")
      motifGainStrand<-paste(as.character(cc$strand),collapse=",")
      motifGainScoreAlt<-paste(as.character(cc$altScore),collapse=",")
      motifGainScoreRef<-paste(as.character(cc$refScore),collapse=",")
      
      motifGainMotifFullDetail<-paste(paste(cc$motifFullName,"[",paste(posIndex[x],cc$motifStart,cc$motifEnd,cc$strand,cc$motifPos,cc$diffScore,sep="#"),"]",sep=""),collapse=",")
      
      maxDiffScore<-max(cc$diffScore)
      selectedIdx<-(cc$diffScore %in% maxDiffScore)
      
      cc3<-cc[selectedIdx,]
      t5<-strsplit(cc3$motifFullName,"\\_")
      cc4<-data.frame(do.call(rbind,t5),stringsAsFactors = FALSE)
      colnames(cc4)<-c("motifNameRoot","motifNameExtend")
      
      #motifGainMotifFullDetailMax<-paste(paste(cc3$motifFullName,"[",paste(posIndex[x],cc$diffScore,sep="#"),"]",sep=""),collapse=",")
      #motifGainMotifTFNameMax<-paste(as.character(unique(cc4$motifNameRoot)),collapse=",")
      
      motifGainMotifFullDetailMax<-paste(paste(cc3$motifFullName,"[",paste(posIndex[x],cc3$motifStart,cc3$motifEnd,cc3$strand,cc3$motifPos,cc3$diffScore,sep="#"),"]",sep=""),collapse=",")
      motifGainMotifTFNameMax<-paste(as.character(unique(cc3$motifFullName)),collapse=",")
      
      
    }else{
      
      motifGainMotifFullName<-NA
      motifGainMotifTFName<-NA
      motifGainStrand<-NA
      motifGainScoreAlt<-NA
      motifGainScoreRef<-NA
      motifGainScoreDiff<-NA
      motifGainMotifFullDetail<-NA
      motifGainMotifFullDetailMax<-NA
      motifGainMotifTFNameMax<-NA
      
      
    }
    
    tmp<-data.frame(motifGainMotifTFName,motifGainMotifFullName,
                    motifGainStrand,motifGainScoreAlt,motifGainScoreRef,
                    motifGainMotifFullDetail,motifGainMotifFullDetailMax,
                    motifGainMotifTFNameMax,
                    stringsAsFactors = FALSE)
    return(tmp)
    
  },mc.cores=useCores)
  
  #motifGainParsed<-list(fullName=motifGainFullName,familyName=motifGainFamilyName,strand=motifGainStrand)
  motifGainParsed<-rbind.fill(tmp)
  
  return(motifGainParsed)  
}
