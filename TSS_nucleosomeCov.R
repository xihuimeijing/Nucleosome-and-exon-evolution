#!/usr/local/bin/Rscript --slave --vanilla 

args <- commandArgs(TRUE)
if(length(args) < 1 || args[1] == "-h"){
  cat("No arguments supplied. Usage:\n")
  helpstr = "TSS_nucleosomeCov.R inFile
  ## inFile  Three columns[pos count frequency]
   \n"
  cat (helpstr)
  quit(status=1) 
}

data<-read.delim(file=args[1],header=F)

mean=''
step<-seq(1, nrow(data), 20)
for(i in 1:(length(step)-1)){
  mean[i]<-mean(data[step[i]:step[i+1],3])
}
result<-cbind(data[step+10,1],mean)
pdf(file="TSS_nucleosomeCov.pdf")
plot(result[,1],result[,2],type="n",xlab="TSS Relative Position",ylab="Nucleosome dyad fraction",lab=c(20,5,5),las=3)
lines(spline(result[,1],result[,2]))
dev.off()