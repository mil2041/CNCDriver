#' Make QQ-plot
#'
#' @param outputDf Output from CNCDriver p-value calculations
#' @return ggplot2 graph object
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
#' @import ggplot2
#' @import ggrepel
#' 


makeQQplot<-function(outputDf){

    outputDf<-outputDf[order(outputDf$pValue),]
    outputDf$padj<-p.adjust(outputDf$pValue,method="BH")
    
    
    threshold<-0.1
    hitNum<-nrow(outputDf[outputDf$padj<threshold,])
    cat(sprintf("hit number: %s\n",hitNum))
    
    #####
    
    pvector<-outputDf$pValue
    #names(pvector)<-outputDf$geneSymbol
    names(pvector)<-outputDf$elementPos
    observed = -log10(sort(pvector,decreasing=F))
    expected = -log10( 1:length(observed)/length(observed) )
    #axisMax<-max(max(observed),max(expected))
    
    idx<-seq(1,length(observed))
    geneSymbol<-names(observed)
    padj<-outputDf$padj
    df<-data.frame(geneSymbol,expected,observed,idx,padj,stringsAsFactors = FALSE)
    
    #df$group<-rep(tumorType,nrow(df))
    #dd[[tumorType]]<-df
    #dd<-rbind.fill(dd)
    
    
    #filePath<-workDir
    #fileName="df_debug.txt"
    #fileName<-file.path(filePath,fileName)
    #write.table(df,fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)
    
    
    #df2<-head(df,30)
    
    #idx<-seq(1,length(observed))
    #geneSymbol<-names(observed)
    #padj<-outputDf$padj
    
    #label<-rep(paste("nearest ",distCutOff*100," %",sep=""),nrow(outputDf))
    #label<-rep(paste("nearest ",categoryRankCutOff*100," %",sep=""),nrow(outputDf))
    #label<-rep(distCutOff,nrow(outputDf))
    
    
    #df[[i]]<-data.frame(geneSymbol,expected,observed,idx,padj,label,stringsAsFactors = FALSE)
    
    
    #}
    
    #df.combined<-rbind.fill(df)
    #x<-as.factor(df.combined$label)
    #x<-factor(x,levels(x)[c(1:9)])
    #print(levels(x))
    #df.combined$triMutMatch<-x
    
    #df2<-head(df,30)
    
    threshold<-0.1
    #axisMax<-max(max(observed),max(expected))+0.05
    axisMax<-7+0.05
    
    ####
    df$group<-rep("notSignificant",nrow(df))
    
    if(sum(df$padj<=threshold)>0){
    df[df$padj<=threshold,]$group<-"significant"
    }
    ####
    
    
    #pdf(paste(workDir,"/",tumorType,"_",mutationType,"_qqplot_iter_",reSampleNum,"_w_label.pdf",sep=""))
    
    #zp1 <- ggplot(df,
    #              aes(x = expected, y = observed, group=label, colour=triMutMatch))
    
    zp1 <- ggplot(df,
                  aes(x = expected, y = observed))
    
    
    zp1 <- zp1 + labs(x = expression(Expected~~-log[10](italic(p))), y = expression(Observed~~-log[10](italic(p))))
    zp1 <- zp1 + xlim(0,axisMax) + ylim(0,axisMax)
    
    zp1 <- zp1 + scale_color_manual(values=c("black","red"))
    #zp1 <- zp1 + ggtitle(paste("QQ plot for ",tumorType," [ iterations: ",reSampleNum," ]\n",mutationType," ",groupName,sep=""))
    #zp1 <- zp1 + ggtitle(paste(tumorType," [ iterations: ",reSampleNum," ]\n",mutationType,", ","q-value= ",threshold,sep=""))
    
    #zp1 <- zp1 + labs(title=paste(tumorType,", [ q-value = ",threshold," ]",sep=""),
    #                  subtitle=paste(mutationType,", replication timing: ",replicationTimingCutOff,sep=""))
    
    zp1 <- zp1 + labs(title=paste(tumorType,sep=""))
    #                  subtitle=paste(mutationType,sep=""))
    
    zp1 <- zp1 + geom_abline(intercept = 0, slope = 1, colour="grey",linetype="solid")
    #zp1 <- zp1 + geom_abline(intercept = 0, slope = 1, colour="red",linetype="dotted")
    
    zp1 <- zp1 + geom_line(size = 0.5, alpha=0.7)
    zp1 <- zp1 + geom_point(aes(color=factor(group)),size=1.8,shape=20)
    #zp1 <- zp1 + geom_point(data=subset(df,padj>threshold),
    #                        aes(color="black"),size=1.5,shape=20,alpha=0.7)
    #zp1 <- zp1 + geom_point(shape = 20, size = 1.5, alpha=0.7)
    
    
    
    #zp1 <- zp1 + scale_color_gradient(low="red",high="grey")
    #zp1 <- zp1 + scale_color_manual(values=brewer.pal(n=9,name="RdBu"))
    #zp1 <- zp1 + geom_abline(intercept = 0, slope = 1, colour="red",linetype="dotted")
    #zp1 <- zp1 + geom_abline(intercept = 0, slope = 1, colour="grey",linetype="solid")
    #zp1 <- zp1 + geom_abline(intercept = 0, slope = 1, colour="red",linetype="dotted")
    
    
    #zp1 <- zp1 + facet_wrap(~group)
    
    #zp1 <- zp1 + geom_text(aes(label=rownames(df)),hjust=0, vjust=0)
    #zp1 <- zp1 + geom_text_repel(aes(label=ifelse(idx<=10,as.character(rownames(df)),'')),hjust=0, vjust=0,angle=(-45))
    
    if(hitNum<50){
    zp1 <- zp1 + geom_text_repel(data=subset(df,padj<threshold),
                                 aes(label=geneSymbol),
                                 fontface="plain",
                                 segment.color="grey",
                                 size = 5,
                                 box.padding = unit(0.35, "lines"),
                                 point.padding = unit(0.3, "lines"),
                                 force=10)
    }
    
    zp1 <- zp1 + theme_bw()
    
    zp1<- zp1 + theme(axis.text = element_text(size = 10),
                      #axis.title = element_text(size = 7, face="bold"),
                      axis.title.x = element_text(size=7),
                      axis.title.y = element_text(size=7),
                      #axis.title.x = element_blank(),
                      #axis.title.y = element_blank(),
                      #panel.border = element_rect(linetype = "solid", colour = "black"),
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line=element_line(colour="black"),
                      plot.title = element_text(lineheight=1.5, face="bold",size=10,hjust=0.5),
                      plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                      legend.position="none")
                      #legend.position = c(0.65,0.2),
                      #legend.background = element_rect(color = NA, 
                      #                                 fill = "transparent", size = 1, linetype = "solid"),
                      #legend.text = element_text(size = 7, colour = "black", angle = 0))
    
    print(zp1)


####

}
