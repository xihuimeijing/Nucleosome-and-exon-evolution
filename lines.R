#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script plot lines in a single plot.\n',file=stderr())
  cat('Usage:lines.R -i=input.tsv -x="Position" -y="Signal value" -o=lines.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe tab delimted input data to plot lines(The first column is the x values and the following comumns are y values.Can be read from STDIN.)\n',file=stderr())
  cat('\t-x1\tINT\tThe xlim start(Optional)\n',file=stderr())
  cat('\t-x2\tINT\tThe xlim end(Optional)\n',file=stderr())
  cat('\t-y1\tINT\tThe ylim start(Optional)\n',file=stderr())
  cat('\t-y2\tINT\tThe ylim end(Optional)\n',file=stderr())
  cat('\t-x\tSTRING\tThe x title[Class]\n',file=stderr())
  cat('\t-y\tSTRING\tThe y title[Value]\n',file=stderr())
  cat('\t-c\tSTRING\tThe colors for each line separated by comma(Optional).\n',file=stderr())
  cat('\t-l\tSTRING\tThe legend names for each line separated by comma(Optional).\n',file=stderr())
  cat('\t-m\tSTRING\tThe title for the plot(Optional).\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput file name[lines.pdf]\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

xlab="class"
ylab="value"
out="lines.pdf"

if(length(args)==0 || args[1]=="-h"){
  usage()
}else{
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the input file -i')
      }else{
        inFile=arg.split[2]
      }
    }else if(grepl('^-x=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      xlab=arg.split[2]
    }else if(grepl('^-y=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      ylab=arg.split[2]
    }else if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      out=arg.split[2]
    }
  }
}

if(exists("inFile")){
	data<-read.delim(file=inFile,comment.char = "#",header=F)
}else{
	data<-read.delim(file('stdin'), header = F)
}
colNum = ncol(data)
x1=min(data[,1])
x2=max(data[,1])
y1=min(data[,-1])
y2=max(data[,-1])
colors<-colors()[1:(colNum-1)]
legends<-as.vector(mapply(paste,"sample",sep='',c(1:(colNum-1))))
main=""
for(i in 1:length(args)){
  arg=args[i]
  if(grepl('^-c=',arg)){
    arg.split = strsplit(arg, '=', fixed = T)[[1]]
    colors=strsplit(arg.split[2],',',fixed = T)[[1]]
  }else if(grepl('^-x1',arg)){
    arg.split=strsplit(arg, '=', fixed = T)[[1]]
    x1=arg.split[2]
  }else if(grepl('^-x2',arg)){
    arg.split=strsplit(arg, '=', fixed = T)[[1]]
    x2=arg.split[2]
  }else if(grepl('^-y1',arg)){
    arg.split=strsplit(arg, '=', fixed = T)[[1]]
    y1=as.numeric(arg.split[2])
  }else if(grepl('^-y2',arg)){
    arg.split=strsplit(arg, '=', fixed = T)[[1]]
    y2=as.numeric(arg.split[2])
  }else if(grepl('^-l',arg)){
    arg.split = strsplit(arg, '=', fixed = T)[[1]]
    legends=strsplit(arg.split[2],',',fixed = T)[[1]]
  }else if(grepl('^-m',arg)){
    arg.split = strsplit(arg, '=', fixed = T)[[1]]
    main=strsplit(arg.split[2],',',fixed = T)[[1]]
  }
}

pdf(file=out)
plot(data[,1],data[,2],xlab=xlab,ylab=ylab,type="l",xlim=c(x1,x2),ylim=c(y1,y2),col=colors[1],main=main)
if(ncol(data)>2){
	for(i in 3:colNum){
  		lines(data[,1],data[,i],col=colors[i-1])
	}
	legend("topright",legends,col=colors,lty=1)
}
dev.off()
