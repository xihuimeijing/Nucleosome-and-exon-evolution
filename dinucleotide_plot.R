#!/usr/local/bin/Rscript --slave --vanilla 

args <- commandArgs(TRUE)
if(length(args) < 2 || args[1] == "-h"){
  cat("No arguments supplied. Usage:\n")
  helpstr = "dinucleostide_plot.R AT_freq.tsv GC_freq.tsv dinucleotide.pdf
  ## AT_freq   frequency of AA/AT/TA/TT dinucleotide with header line to be'pos AA TT AT TA',5 columns
  ## GC_freq   frequency of GG/GC/CG/CC dinucleotide the same format with AT_freq
  ## Output    One file for output figures \n"
  
  cat (helpstr)
  quit(status=1) 
}
ATfile=args[1]
GCfile=args[2]

library("ggplot2")
library("reshape2")
AT<-read.delim(file=ATfile,header=T)
AT<-AT[order(AT[,1]),]
AT_mean<-rowMeans(AT[,-1])
GC<-read.delim(file=GCfile,header=T)
GC<-GC[order(GC[,1]),]
GC_mean<-rowMeans(GC[,-1])

mean<-as.data.frame(cbind(AT[,1],AT_mean,GC_mean))
colnames(mean)<-c("pos","AA/AT/TA/TT","GG/GC/CG/CC")
mean<-melt(mean,id.vars="pos")
pdf(file=args[3])
ggplot(mean,aes(x=pos,y=value))+geom_line(aes(col=variable))+labs(x="Position relative to nucleosome dyad",
       y="Dinucleotide frequency",col="")+theme(panel.grid.minor.x=element_line(size=0.5,
      color="white"))+scale_x_continuous(minor_breaks = seq(-1000, 1000, 10),breaks=c(-73,0,73),labels=c("-73","0","73"))

sum<-as.data.frame(cbind(AT,GC)[,-6])
sum<-melt(sum,id.vars="pos")
ggplot(sum,aes(x=pos,y=value))+geom_line(aes(col=variable))+labs(x="Position relative to nucleosome dyad",
     y="Dinucleotide frequency",col="")+theme(panel.grid.minor.x=element_line(size=0.5,
     color="white"))+scale_x_continuous(minor_breaks = seq(-1000, 1000, 10),breaks=c(-73,0,73),labels=c("-73","0","73"))
dev.off()