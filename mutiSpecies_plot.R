setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/GC_evolution/hg19_regions")
nc<-read.delim(file="H.Hancestor.intergenic.NCgt2.GC",header=F)
nf<-read.delim(file="H.Hancestor.intergenic.NClt0.5.GC",header=F)
par(mfrow=c(2,2))
plot(nc[,2],nc[,3],xlim=c(0,1),ylim=c(0,1),xlab="human GC",ylab="human ancestor GC",main="Intergenic RTMS NC>2")
abline(a=0,b=1,col="red")
plot(nf[,2],nf[,3],xlim=c(0,1),ylim=c(0,1),xlab="human GC",ylab="human ancestor GC",main="Intergenic RTMS NC<0.5")
abline(a=0,b=1,col="red")
boxplot(nc[,2],nc[,3],nf[,2],nf[,3],col=c("darkgreen","cadetblue1"),boxwex=0.3,at=sort(c(1:2-0.2,1:2+0.2)),ylab="GC content",names=c("NC>2\nhg19","NC>2\nancestor","NC<0.5\nhg19","NC<0.5\nancestor"),main="Intergenic 100bp window")
plot(ecdf(nc[,2]-nc[,3]),verticals = T, do.points=F, col="red", xlab="GC difference(hg19-ancestor)",main="Intergenic 100bp window")
lines(ecdf(nf[,2]-nf[,3]),verticals = T, do.points=F,col="blue")
legend("topleft",c("NC>2","NC<0.5"),col=c("red","blue"),lty=1)

#intron
nc<-read.delim(file="H.Hancestor.intron.NCgt2.GC",header=F)
nf<-read.delim(file="H.Hancestor.intron.NClt0.5.GC",header=F)
plot(nc[,2],nc[,3],xlim=c(0,1),ylim=c(0,1),xlab="human GC",ylab="human ancestor GC",main="Intron RTMS NC>2")
abline(a=0,b=1,col="red")
plot(nf[,2],nf[,3],xlim=c(0,1),ylim=c(0,1),xlab="human GC",ylab="human ancestor GC",main="Intron RTMS NC<0.5")
abline(a=0,b=1,col="red")
plot(ecdf(nc[,2]-nc[,6]),verticals = T, do.points=F, col="red", xlab="GC difference(hg19-susScr3)",main="Intron 100bp window")
lines(ecdf(nf[,2]-nf[,6]),verticals = T, do.points=F,col="blue")
legend("topleft",c("NC>2","NC<0.5"),col=c("red","blue"),lty=1)



gt1<-read.delim(file="HRTMS.intron.NCgt1.GC",header=F)
total<-read.delim(file="HRTMS.intron.GC",header = F)
plot(ecdf(total[,2]-total[,6]),verticals = T, do.points=F, col="gray")
lines(ecdf(gt1[,2]-gt1[,6]),verticals = T, do.points=F,col="blue")

setwd("./hg19_rheMac2/")
par(mfrow=c(1,2))
calss1<-read.delim(file="intergenic.class1.anestor.hg19.GC.tsv",header=F)
calss2<-read.delim(file="intergenic.class2.anestor.hg19.GC.tsv",header=F)
calss3<-read.delim(file="intergenic.class3.anestor.hg19.GC.tsv",header=F)
calss4<-read.delim(file="intergenic.class4.anestor.hg19.GC.tsv",header=F)
boxplot(calss1$V6,calss2$V6,calss3$V6,calss4$V6,names=c("class1","class2","class3","class4"),ylab="GC content",main="Intergenic")
calss1<-read.delim(file="intron.class1.anestor.hg19.GC.tsv",header=F)
calss2<-read.delim(file="intron.class2.anestor.hg19.GC.tsv",header=F)
calss3<-read.delim(file="intron.class3.anestor.hg19.GC.tsv",header=F)
calss4<-read.delim(file="intron.class4.anestor.hg19.GC.tsv",header=F)
#Control CpG
setwd("./controlCpG/")
class1<-read.delim(file="intergenic.class1.ancestor.hg19.CpG.tsv",header = F)
class4<-read.delim(file="intergenic.class4.ancestor.hg19.CpG.tsv",header=F)
par(mfrow=c(1,2))
boxplot(class1$V5,class4$V5,outline = F,names=c("Class1","Class4"),ylab="CpG content",border = c("red","blue"),main="All regions",ylim=c(0,0.15))
data<-read.delim(file="intergenic.class1.class4.comparible.ancestor.list",header=F)
boxplot(data$V5,data$V11,outline = F,names=c("Class1","Class4"),ylab="CpG content",border = c("red","blue"),main="Control CpG",ylim=c(0,0.15))

#Control GC
setwd("./ControlGC/")
data<-read.delim(file="intergenic.class4.class1.comparible.ancestor.list",header=F)
boxplot(data$V11,data$V5,data$V12,data$V6,col=c("gray","white"),at=sort(c(1:2-0.2,1:2+0.2)),outline=F,boxwex=0.3,names = c("A_class1","A_class2","H_class1","H_class2"),ylab="GC content")
points(c(0.8,1.2,1.8,2.2),c(0.4514,0.4509,0.4509,0.4497),col="black",pch=19)
subset<-data[data$V5!=data$V6 | data$V11!=data$V12,]
boxplot(subset$V11,subset$V5,subset$V12,subset$V6,col=c("gray","white"),at=sort(c(1:2-0.2,1:2+0.2)),outline=F,boxwex=0.3,names = c("A_class1","A_class2","H_class1","H_class2"),ylab="GC content")
ggplot(data,aes(x=V11, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,500))+labs(x="Ancestor GC content(Class1)",y="Ancestor GC content(Class2)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits = c(0,1))
ggplot(data,aes(x=V12, y=V6))+stat_bin_2d(bins=100)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,500),na.value="red")+labs(x="Human GC content(Class1)",y="Human GC content(Class2)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits = c(0,1))+geom_abline(slope = 1,intercept = 0, col="black")
ancestor=subset$V11/subset$V5
human=subset$V12/subset$V6
data<-cbind(ancestor,human)
colnames(data)<-c("ancestor","human")
data<-melt(as.data.frame(data))
ggplot(data,aes(x=variable, y=value))+geom_violin(trim=F,fill="gray",color=NA)+labs(x="",y="GC content ratio")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  stat_summary(fun.y=mean, geom="point", size=2,color="white")+ylim(0.6,2)

data<-read.delim(file="intergenic.class4.class1.comparible.ancestor.newList",header=F)
data<-data[,c(5,11)]
colnames(data)<-c("class2","class1")
data<-melt(data)
data$variable<-as.factor(data$variable)
ggplot(data,aes(x=variable,y=value,fill=variable))+geom_violin(trim = F,fill=NA)+geom_boxplot(width=0.2)+scale_fill_grey(start=0.5,end=1.0)+labs(y="Ancestor GC content",x="")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+ylim(0,1.2)
human<-read.delim(file="H-R-newAnces-hg19.txt",header=T)
human<-ggplot(human,aes(x=mutation,y=rate,fill=Class))+geom_bar(stat="identity",color="black",width=0.6,position=position_dodge())+scale_fill_grey(start=0.5,end=1.0)+geom_errorbar(aes(ymin=rate-sd,ymax=rate+sd),width=0.1,position=position_dodge(0.6))+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
monkey<-read.delim(file="H-R-newAnces-rheMac2.txt",header=T)
monkey<-ggplot(monkey,aes(x=mutation,y=rate,fill=Class))+geom_bar(stat="identity",color="black",width=0.6,position=position_dodge())+scale_fill_grey(start=0.5,end=1.0)+geom_errorbar(aes(ymin=rate-sd,ymax=rate+sd),width=0.1,position=position_dodge(0.6))+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
grid.arrange(human,monkey,ncol=1,nrow=2) 

ancestor<-read.delim(file="NCdiff.adjacent.3ss.ancestor.tsv",header=F)
human<-read.delim(file="NCdiff.adjacent.3ss.hg19.tsv",header=F)
ancestor<-ancestor[,-c(1:6)]
human<-human[,-c(1:6)]
change<-human-ancestor
pos=0
for(i in 1:nrow(change)){
  pos[i]<-which.max(change[i,])
}

diff<-apply(as.matrix(human)-as.matrix(ancestor),1,FUN=max)
class1<-read.delim(file="class1.adjacent.3ss.hg19.ancestor.tsv",header=F)
diff_1<-apply(class1[,c(4:181)]-class1[,c(185:362)],1,FUN=max)
class4<-read.delim(file="class4.adjacent.3ss.hg19.ancestor.tsv",header=F)
diff_4<-apply(class4[,c(4:181)]-class4[,c(185:362)],1,FUN=max)



setwd("../position/")
data<-read.delim(file="n147.ancestor.hg19.profile.txt",header = T)
names(data)<-c("pos","Ancestor","human")
plot(data$pos,data$Ancestor,type="l",col="black")
lines(data$pos,data$human,col="red")
plot(data$pos,data$human-data$Ancestor,type="l")
data<-read.delim(file="speciesMean.n147.sdLt20.ancestor.hg19.GCprofile.txt",header=T)

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/NRL/")
data<-read.delim(file="NRL.summary.tsv",header = T)
ymin<-data$NRL-data$SE
ymax<-data$NRL+data$SE
ggplot(data, aes(x = species, y = NRL, fill = tissue))+geom_bar(position=position_dodge(0.9), 
    stat="identity",width=0.8)+geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .2,position=position_dodge(0.9))+coord_cartesian(ylim=c(180,200))+scale_y_continuous(minor_breaks=seq(1:200))
H1C<-read.delim(file="../genExp/H1C_FPKM/H1C.FPKM.sum",header = F)
data<-cbind(data[order(data$species,data$tissue),],H1C[order(H1C$V1,H1C$V2),])
data<-data[order(data$V3),]
rownames(data)<-paste(data[,1],data[,2],sep="-")
data[,7]<-log10(data[,7])
pheatmap(data[,c(3,7)],cluster_rows = F,cluster_cols = F,scale = "column",show_rownames = T,border_color =NA)
sd=0
sd[1]<-sd(data[data$tissue=="brain",3])
sd[2]<-sd(data[data$tissue=="heart",3])
sd[3]<-sd(data[data$tissue=="kidney",3])
sd[4]<-sd(data[data$tissue=="liver",3])
sd[5]<-sd(data[data$tissue=="muscle",3])
sd[1]<-sd(data[data$species=="human",3])
sd[2]<-sd(data[data$species=="monkey",3])
sd[3]<-sd(data[data$species=="treeShrew",3])
sd[4]<-sd(data[data$species=="mouse",3])
sd[5]<-sd(data[data$species=="pig",3])

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/occupancy")
data<-read.delim(file="5species.homo.n100.pearsonCor.tsv",header = 1,row.names = 1)
rownames(data)<-gsub("'","",rownames(data))
colnames(data)<-rownames(data)
pheatmap(data,cluster_rows=T, cluster_cols=T)
pheatmap(data,cluster_rows=T, cluster_cols=T, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation")
data<-read.delim(file="5species.homo.n147.pearsonCor.tsv",header=1,row.names = 1)
annotation<-as.data.frame(strsplit(rownames(data),"-"))
rownames(annotation)<-c("species","tissue")
colnames(annotation)<-colnames(data)
annotation<-as.data.frame(t(annotation))
myColors=rainbow(10)
ann_colors=list(tissue=c(brain=myColors[1],heart=myColors[2],kidney=myColors[3],liver=myColors[9],muscle=myColors[10])
                ,species=c(human=myColors[4],monkey=myColors[6],tree=myColors[7],mouse=myColors[8],pig="black"))
pheatmap(data,cluster_rows=T, cluster_cols=T, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",show_rownames = F,show_colnames = F
         ,annotation_row = annotation,annotation_col = annotation, annotation_colors = ann_colors,border_color = NA)
pca<-read.delim(file="5species.homo.n100.PCA.tsv",header = T)
percent<-round(100*pca[c(1,2),15]/sum(pca$Eigenvalue),digits = 2)
data<-t(pca[c(1,2),c(2:14)])
data<-cbind(t(as.data.frame(strsplit(rownames(data),'.',fixed = T))),data)
colnames(data)<-c("species","tissue","PC1","PC2")
data<-as.data.frame(data)
ggplot(data, aes(PC1, PC2, color=species, shape=tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percent[1],"% variance")) +
  ylab(paste0("PC2: ",percent[2],"% variance")) + 
  scale_x_continuous(breaks =seq(-0.1,0.5,0.1),labels=seq(-0.1,0.5,0.1))+
  scale_y_continuous(breaks =seq(-0.3,0.5,0.6),labels=seq(-0.3,0.5,0.6))
##human and monkey
data<-read.delim(file="human.monkey.n147.tsv",header=F)
dcols<-densCols(data$V1,data$v2,colramp = colorRampPalette(c('blue', 'yellow', 'red')))
plot(data$V1,data$V2,col=dcols,pch=19,xlab="human NC(147bp window)",ylab="monkey NC")
brain_d<-data$V1-data$V3
brain_d[brain_d<(-5)]=-5
brain_d[brain_d>5]=5
h<-hist(brain_d,plot=F,breaks=100)
muscle_d<-data$V2-data$V4
muscle_d[muscle_d<(-5)]=-5
muscle_d[muscle_d>5]=5
h<-hist(muscle_d,plot=F,breaks=100)
par(mfrow=c(1,2))
h$counts <- h$counts / sum(h$counts)#correlation: brain0.42,muscle0.29
plot(h, freq=TRUE, ylab="Frequency",xlab="Muscle NC(human-monkey)",main="",xlim=c(-5,5),ylim=c(0,0.07))
#Whole genome position
setwd("D:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/position")
data<-read.delim(file="homo.n147.pearsonCor.tsv",header = T,row.names = 1)
subData<-data[,-c(8,9,15)]
subData<-subData[-c(8,9,15),]
annotation<-as.data.frame(strsplit(rownames(data),"_"))
rownames(annotation)<-c("species","tissue")
colnames(annotation)<-colnames(data)
annotation<-as.data.frame(t(annotation))
myColors=rainbow(10)
ann_colors=list(tissue=c(brain=myColors[1],heart=myColors[2],kidney=myColors[3],liver=myColors[9],muscle=myColors[10])
                ,species=c(human=myColors[4],monkey=myColors[6],tree=myColors[7],mouse=myColors[8],pig="black"))
pheatmap(subData,cluster_rows=T, cluster_cols=T, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",show_rownames = F,show_colnames = F
         ,annotation_row = annotation,annotation_col = annotation, annotation_colors = ann_colors, border_color = NA)
annotation<-as.data.frame(strsplit(rownames(data),"_"),stringsAsFactors = F)
annotation[8,2]<-"liver"
annotation[9,2]<-"liver"
rownames(annotation)<-c("species","tissue")
colnames(annotation)<-colnames(data)
annotation<-as.data.frame(t(annotation))
myColors=rainbow(10)
ann_colors=list(tissue=c(brain=myColors[1],heart=myColors[2],kidney=myColors[3],liver=myColors[9],muscle=myColors[10],sperm="brown")
                ,species=c(human=myColors[4],monkey=myColors[6],tree=myColors[7],mouse=myColors[8],pig="black"))
pheatmap(data,cluster_rows=T, cluster_cols=T, clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",show_rownames = F,show_colnames = F
         ,annotation_row = annotation,annotation_col = annotation, annotation_colors = ann_colors)

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/occupancy/homo_exons/")
data<-read.delim(file="meanSummary.tsv",header=F)
names=c("human-brain","human-muscle","mouse-brain","mouse-heart","mouse-kidney","mouse-liver","mouse-muscle","monkey-brain","monkey-heart","monkey-kidney","monkey-liver","monkey-muscle","pig-brain","pig-heart","pig-kidney","pig-liver","pig-muscle","tree-brain","tree-heart","tree-kidney","tree-muscle")
colnames(data)<-c("exon",names)
data<-data[-7516,]
cor<-cor(data[,-1])
pheatmap(cor,cluster_rows=T, cluster_cols=T)
cor<-cor(data[,-1],method="spearman")
data<-read.delim(file="upNCratio.sum.tsv",header=F)
data<-data[,-1]
tmp<-as.vector(as.matrix(data))#min:-11.208;max: 10.63
max(tmp[is.finite(tmp)])
min(tmp[is.finite(tmp)])
tmp[grepl("-Inf",tmp)]<-(-15)
tmp[grepl("Inf",tmp)]<-15
data<-matrix(tmp,ncol=ncol(data),byrow = F)
colnames(data)<-names
cor<-cor(data,use="pairwise.complete.obs")
cor<-cor(data,use="pairwise.complete.obs",method="spearman")

data<-read.delim(file="downNCratio.sum.tsv",header=F)
data<-data[,-1]
#Exon Usage
data<-read.delim(file="exonUsage.sum.uniq.tsv",header=T,row.names = 1)
data[data<0]=NA
cor<-cor(as.matrix(data),use="pairwise.complete.obs")
pheatmap(cor,cluster_rows = T,cluster_cols = T,show_rownames = T,show_colnames = T)
#In the same tissues
data<-read.delim(file="muscle.upNCratio.usage.tsv",header=F)
data[data$V8<0,8]=NA
data[data$V9<0,9]=NA
data[data$V10<0,10]=NA
data[data$V11<0,11]=NA
data[data$V12<0,12]=NA
cor=0
for(i in 1:nrow(data)){
  cor[i]=cor(as.vector(as.numeric(data[i,c(2:6)])),as.vector(as.numeric(data[i,c(8:12)])),use="pairwise.complete.obs")
}
plot(ecdf(cor),col="blue",do.points=F,verticals=TRUE,xlim=c(-1,1),main="",ylab="Cumlative fraction",xlab="Pearson correlation of psi and exon/upIntron logNC")
abline(a=0.5,b=0.5,col="black",lty=1)
ks.test()
cor_s=0
for(i in 1:nrow(data)){
  cor_s[i]=cor(as.vector(as.numeric(data[i,c(2:6)])),as.vector(as.numeric(data[i,c(8:12)])),use="pairwise.complete.obs",method="spearman")
}
#Switch score
usage=0
nc=0
num=0
for(i in 1:nrow(data)){
  max<-max(data[i,c(8:12)],na.rm = T)
  min<-min(data[i,c(8:12)],na.rm = T)
  max_i<-which(data[i,c(8:12)] == max(data[i,c(8:12)], na.rm = TRUE))
  min_i<-which(data[i,c(8:12)] == min(data[i,c(8:12)], na.rm = TRUE))
  if(length(max_i)==1 & length(min_i)==1){
    num=num+1
    usage[num]=max-min
    nc[num]=data[i,1+max_i]-data[i,1+min_i]
  }
}
#No significant result!
setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/occupancy/homo_genes/")
data<-read.delim(file="TSSprofile.sum.txt",header=F)
names2<-c("monkey-liver","tree-heart","tree-kidney","tree-muscle","mouse-brain","mouse-heart","mouse-kidney","mouse-liver","mouse-muscle","pig-brain","pig-heart","pig-kidney","pig-liver")
colnames(data)<-c("pos",names2)
cor<-cor(data[,-1])
pheatmap(cor,cluster_rows=T, cluster_cols=T, display_numbers=T, number_color = "black")
cor<-cor(data[data$pos>=(-500) & data$pos<500,-1])
setwd("./version2")
pos=''
for(i in 1:13){
  pos[i]<-which.max(data[c(1500:1700),i])
}
minus=''
for(i in 1:13){
  minus[i]<-which.max(data[c(1000:1500),i])
}
#NFR & plusNuc-NFR
NFR<-read.delim(file="NFR.sum.tsv",header=F)
colnames(NFR)<-c("exon",names)
pheatmap(cor(as.matrix(NFR[,-1]),use="pairwise.complete.obs"),cluster_rows=T, cluster_cols=T, display_numbers=T, number_color = "black")
plus<-read.delim(file="plusNucleosome.sum.tsv",header=F)
colnames(plus)<-c("exon",names)
data<-plus[,-1]-NFR[,-1]
colnames(data)<-names
pheatmap(cor(data,use="pairwise.complete.obs"),cluster_rows=T, cluster_cols=T, display_numbers=T, number_color = "black")
#expression
setwd("./expression/")
data<-read.delim(file="FPKM.summary.tsv",header=T)
cor<-cor(data[,c(3:ncol(data))])
pheatmap(cor,cluster_rows = T,cluster_cols = T,show_rownames = T,show_colnames = T)

#position
data<-read.delim(file="summary.distance.filter.tsv",header=F)
colnames(data)<-c("chr","start","end","monkey-liver","mouse-brain","mouse-heart","mouse-kidney","mouse-liver","mouse-muscle","pig-brain","pig-heart","pig-kidney","pig-liver","tree-heart","tree-kidney","tree-muscle")
cor<-cor(data[,-c(1,2,3)],use="pairwise.complete.obs")
pheatmap(cor,cluster_rows=T, cluster_cols=T, display_numbers=T, number_color = "black")

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/exonAge")
brain<-read.delim(file="brain.colMeanNor.tsv",header=F)
brain<-read.delim(file="brain.logNC.tsv",header=F)
boxplot(brain[brain$V2=="H----",3],brain[brain$V2=="H----",4],brain[brain$V2=="H----",5],brain[brain$V2=="H----",6],brain[brain$V2=="H----",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
boxplot(brain[brain$V2=="HR---",3],brain[brain$V2=="HR---",4],brain[brain$V2=="HR---",5],brain[brain$V2=="HR---",6],brain[brain$V2=="HR---",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
boxplot(brain[brain$V2=="HRT--",3],brain[brain$V2=="HRT--",4],brain[brain$V2=="HRT--",5],brain[brain$V2=="HRT--",6],brain[brain$V2=="HRT--",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
boxplot(brain[brain$V2=="HRTM-",3],brain[brain$V2=="HRTM-",4],brain[brain$V2=="HRTM-",5],brain[brain$V2=="HRTM-",6],brain[brain$V2=="HRTM-",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
boxplot(brain[brain$V2=="HRTMS",3],brain[brain$V2=="HRTMS",4],brain[brain$V2=="HRTMS",5],brain[brain$V2=="HRTMS",6],brain[brain$V2=="HRTMS",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")

wilcox.test(brain[brain$V2=="H----",3],brain[brain$V2=="H----",4],alternative = "g")
wilcox.test(brain[brain$V2=="H----",3],brain[brain$V2=="H----",5],alternative = "g")
wilcox.test(brain[brain$V2=="H----",3],brain[brain$V2=="H----",6],alternative = "g")
wilcox.test(brain[brain$V2=="H----",3],brain[brain$V2=="H----",7],alternative = "g")
muscle<-read.delim(file="muscle.logNC.tsv",header=F)
boxplot(muscle[muscle$V2=="H----",3],muscle[muscle$V2=="H----",4],muscle[muscle$V2=="H----",5],muscle[muscle$V2=="H----",6],muscle[muscle$V2=="H----",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
boxplot(muscle[muscle$V2=="HRTM-",3],muscle[muscle$V2=="HRTM-",4],muscle[muscle$V2=="HRTM-",5],muscle[muscle$V2=="HRTM-",6],muscle[muscle$V2=="HRTM-",7],outline = F, notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
wilcox.test(muscle[muscle$V2=="H----",3],muscle[muscle$V2=="H----",4],alternative = "g")
wilcox.test(muscle[muscle$V2=="H----",3],muscle[muscle$V2=="H----",5],alternative = "g")
wilcox.test(muscle[muscle$V2=="H----",3],muscle[muscle$V2=="H----",6],alternative = "g")
wilcox.test(muscle[muscle$V2=="H----",3],muscle[muscle$V2=="H----",7],alternative = "g")

##splice score
ss<-read.delim(file="H----.5ss.sum.tsv",header=F)
boxplot(ss[,-1],outline = F,notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="5' splice score")
ss<-read.delim(file="H----.3ss.sum.tsv",header=F)
boxplot(ss[,-1],outline = F,notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="3' splice score")

##GC ratio
data<-read.delim(file="H----.up.logGC.sum.tsv",header=F)
boxplot(data[,-1],outline = F, notch = T, names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron GC)")

## Normlazated by GC ratio
NC<-read.delim(file="brain.logNC.tsv",header=F)
GC<-read.delim(file="H----.up.logGC.sum.tsv",header=F)
norm=NC[NC$V2=="H----",c(3:7)]-GC[,c(2:6)]
boxplot(GC[,-1],outline = F, names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron GC)")
wilcox.test(GC$V4,GC$V5,alternative = "g")
#Tissue mean of logNC
data<-read.delim(file="5species.logNC.tissueMean.tsv",header=F)
vioplot(data[data$V2=="H----",3],data[data$V2=="H----",4],data[data$V2=="H----",5],data[data$V2=="H----",6],data[data$V2=="H----",7],col="gray",c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron NC)")
vioplot(data[data$V2=="HRTMS",3],data[data$V2=="HRTMS",4],data[data$V2=="HRTMS",5],data[data$V2=="HRTMS",6],data[data$V2=="HRTMS",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
vioplot(data[data$V2=="HRTM-",3],data[data$V2=="HRTM-",4],data[data$V2=="HRTM-",5],data[data$V2=="HRTM-",6],data[data$V2=="HRTM-",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
vioplot(data[data$V2=="HRT--",3],data[data$V2=="HRT--",4],data[data$V2=="HRT--",5],data[data$V2=="HRT--",6],data[data$V2=="HRT--",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
vioplot(data[data$V2=="HR---",3],data[data$V2=="HR---",4],data[data$V2=="HR---",5],data[data$V2=="HR---",6],data[data$V2=="HR---",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
wilcox.test(data[data$V2=="HR---",5],data[data$V2=="HR---",6],alternative = "g")

#NC prediction
setwd("./NC_prediction/")
data<-read.delim(file="5species.uplogNC.tsv",header=F)
boxplot(data[,-1],notch = T,names = c("human","monkey","treeShrew","mouse","pig"),ylab="Predicted log2(Exon/upIntron NC)")
text(c(1.5,2.5,3.5,4.5),c(2,2,2,2),labels = c(0.13,2.352e-15,6.854e-05,1.5e-4))

setwd("./list_v3/")
gc<-read.delim(file="summary.uplogGC.tsv",header=F)
nc<-read.delim(file="summary.uplogNC.tsv",header = F)
vioplot(nc[nc$V2=="H----",3],nc[nc$V2=="H----",4],nc[nc$V2=="H----",5],nc[nc$V2=="H----",6],nc[nc$V2=="H----",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
vioplot(nc[nc$V2=="HR---",3],nc[nc$V2=="HR---",4],nc[nc$V2=="HR---",5],nc[nc$V2=="HR---",6],nc[nc$V2=="HR---",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
text(c(1.5,2.5,3.5,4.5),c(5,5,5,5),c("NS","NS","NS","p=1.169e-06"))
vioplot(nc[nc$V2=="HRT--",3],nc[nc$V2=="HRT--",4],nc[nc$V2=="HRT--",5],nc[nc$V2=="HRT--",6],nc[nc$V2=="HRT--",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
text(c(1.5,2.5,3.5,4.5),c(5,5,5,5),c("NS","0.006","0.0007","1.888e-06"))
vioplot(nc[nc$V2=="HRTM-",3],nc[nc$V2=="HRTM-",4],nc[nc$V2=="HRTM-",5],nc[nc$V2=="HRTM-",6],nc[nc$V2=="HRTM-",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))
text(c(1.5,2.5,3.5,4.5),c(5,5,5,5),c("NS","0.03","0.001","p<2.2e-16"))
vioplot(nc[nc$V2=="HRTMS",3],nc[nc$V2=="HRTMS",4],nc[nc$V2=="HRTMS",5],nc[nc$V2=="HRTMS",6],nc[nc$V2=="HRTMS",7],col="gray",names=c("human","monkey","treeShrew","mouse","pig"))

boxplot(gc[gc$V2=="H----",3],gc[gc$V2=="H----",4],gc[gc$V2=="H----",5],gc[gc$V2=="H----",6],gc[gc$V2=="H----",7],ylab="log2(Exon/upIntron GC)",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
boxplot(gc[gc$V2=="HR---",3],gc[gc$V2=="HR---",4],gc[gc$V2=="HR---",5],gc[gc$V2=="HR---",6],gc[gc$V2=="HR---",7],ylab="log2(Exon/upIntron GC)",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
text(c(1.5,2.5,3.5,4.5),c(5,5,5,5),c("NS","0.005","0.002","0.03"))
boxplot(gc[gc$V2=="HRT--",3],gc[gc$V2=="HRT--",4],gc[gc$V2=="HRT--",5],gc[gc$V2=="HRT--",6],gc[gc$V2=="HRT--",7],ylab="log2(Exon/upIntron GC)",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
boxplot(gc[gc$V2=="HRTM-",3],gc[gc$V2=="HRTM-",4],gc[gc$V2=="HRTM-",5],gc[gc$V2=="HRTM-",6],gc[gc$V2=="HRTM-",7],ylab="log2(Exon/upIntron GC)",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
boxplot(gc[gc$V2=="HRTMS",3],gc[gc$V2=="HRTMS",4],gc[gc$V2=="HRTMS",5],gc[gc$V2=="HRTMS",6],gc[gc$V2=="HRTMS",7],ylab="log2(Exon/upIntron GC)",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)

wilcox.test(gc[gc$V2=="HR---",3],gc[gc$V2=="HR---",7],alternative = "g")

ss5<-read.delim(file="summary.5ss.tsv",header=F)
ss3<-read.delim(file="summary.3ss.tsv",header=F)
ss<-cbind(ss5[,c(1,2)],(ss5[,c(3:7)]+ss3[,c(3:7)])/2)
boxplot(ss[ss$V2=="H----",3],ss[ss$V2=="H----",4],ss[ss$V2=="H----",5],ss[ss$V2=="H----",6],ss[ss$V2=="H----",7],ylab="Splice score",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
text(c(1.5,2.5,3.5,4.5),c(10,10,10,10),c("p<2.2e-16","p<2.2e-16","p<2.2e-16","p<2.2e-16"))
boxplot(ss[ss$V2=="HR---",3],ss[ss$V2=="HR---",4],ss[ss$V2=="HR---",5],ss[ss$V2=="HR---",6],ss[ss$V2=="HR---",7],ylab="Splice score",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
wilcox.test(ss[ss$V2=="HR---",3],ss[ss$V2=="HR---",4],alternative = "g")
text(c(1.5,2.5,3.5,4.5),c(10,10,10,10),c("p=0.25","p=2.4e-12","p<2.2e-16","p=2.343e-12"))
boxplot(ss[ss$V2=="HRT--",3],ss[ss$V2=="HRT--",4],ss[ss$V2=="HRT--",5],ss[ss$V2=="HRT--",6],ss[ss$V2=="HRT--",7],ylab="Splice score",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
text(c(1.5,2.5,3.5,4.5),c(10,10,10,10),c("NS","NS","p=3.962e-05","p=0.06"))
boxplot(ss[ss$V2=="HRTM-",3],ss[ss$V2=="HRTM-",4],ss[ss$V2=="HRTM-",5],ss[ss$V2=="HRTM-",6],ss[ss$V2=="HRTM-",7],ylab="Splice score",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
text(c(1.5,2.5,3.5,4.5),c(10,10,10,10),c("NS","NS","NS","P=1.427e-06"))
boxplot(ss[ss$V2=="HRTMS",3],ss[ss$V2=="HRTMS",4],ss[ss$V2=="HRTMS",5],ss[ss$V2=="HRTMS",6],ss[ss$V2=="HRTMS",7],ylab="Splice score",names=c("human","monkey","treeShrew","mouse","pig"),outline=F)
text(c(1.5,2.5,3.5,4.5),c(10,10,10,10),c("NS","NS","NS","P=1.427e-06"))
ancestor<-read.delim(file="ss_motif/H----.hg19.ancestorNew.ssScore.tsv",header=F)
boxplot(ancestor$V8,ss3[ss3$V2=="H----",3],ss3[ss3$V2=="H----",4],names=c("ancestor","human","monkey"),ylab="3'ss splice score")
#Read count
x<-matrix(c(237,330,106,461,206,361,157,410),nrow=2)
barplot(x,names.arg =c("monkey","treeShrew","mouse","pig"),ylab="Exon number",ylim=c(0,600))
legend("topright",c("With","Without"),col=grey.colors(2),border="white",pch=15)
#NC prediction(cis directed NC)
data<-read.delim(file="H----.5species.predict.uplogNC.tsv",header=F)
vioplot(data$V2,data$V3,data$V4,data$V5,data$V6,names = c("human","monkey","treeShrew","mouse","pig"),col="gray")
wilcox.test(data$V2,data$V3,alternative = "g",paired = T)

setwd("./final_list/")
##splice score
ss<-read.delim(file="H----.summary.5ss.tsv",header=F)
boxplot(ss[,-c(1,2)],outline = F,names = c("human","monkey","treeShrew","mouse","pig"),ylab="5' splice score")
ss<-read.delim(file="H----.summary.3ss.tsv",header=F)
boxplot(ss[,-c(1,2)],outline = F,names = c("human","monkey","treeShrew","mouse","pig"),ylab="3' splice score")
##GC ratio
data<-read.delim(file="H----.summary.uplogGC.tsv",header=F)
boxplot(data[,-c(1,2)],outline = F, names = c("human","monkey","treeShrew","mouse","pig"),ylab="log2(Exon/upIntron GC)")
wilcox.test(data$V3,data$V4,alternative = "g")#0.002
wilcox.test(data$V3,data$V5,alternative = "g")#p<2.2e-16
wilcox.test(data$V3,data$V6,alternative = "g")#p<2.2e-16
wilcox.test(data$V3,data$V7,alternative = "g")#p=1.708e-11
##NC ratio
data<-read.delim(file="H----.summary.uplogNC.tsv",header=F)
data<-data[,c(3:7)]
colnames(data)<-c("human","monkey","treeShrew","mouse","pig")
data<-melt(data)
ggplot(data,aes(x=variable, y=value, fill=variable))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-5,5)+labs(y="log2(Exon/upIntron NC)",x="")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(data[data$variable=="human",2],data[data$variable=="treeShrew",2],alternative = "g",paired = T)
#Pvalues(one tail, paired): 0.1459,0.00565,2.339e-08,<2.2e-16
###Final version
data<-read.delim(file="brain.5species.uplogNC.final.tsv",header=F,row.names = 1)
data<-read.delim(file="NCratio/brain.5tissues.final.downLogNC.tsv",header=F,row.names = 1)
boxplot(value~variable,data=data,outline=F,ylim=c(-4,5),ylab="Exon/up-intron NC log ratio",par(bty='n'))
#Pvalues(one tail, paired): 0.13,0.02,0.004,0.01; down: 0.17,0.0005547,3.926e-06,0.03
data<-read.delim(file="H----.summary.uplogGC.final.tsv",header=F)
boxplot(data[,-c(1,2)],outline = F, names = c("human","monkey","treeShrew","mouse","pig"),ylab="Exon/up-intron GC log ratio",ylim=c(-0.6,1.2),par(bty='n'))
#Pvalues(one tail, paired): 0.01,3.304e-11,6.836e-16,1.144e-08; down:  7.7e-10,<2.2e-16, <2.2e-16,<2.2e-16
ss<-read.delim(file="H----.summary.3ss.final.tsv",header=F)
boxplot(ss[,-c(1,2)],outline = F,names = c("human","monkey","treeShrew","mouse","pig"),ylab="3' splice score",ylim=c(-30,20),par(bty='n'))
#Pvalues(one tail, paired): 4.627e-11,< 2.2e-16,< 2.2e-16,< 2.2e-16
###H---- in muscle
up<-read.delim(file="muscle.5tissues.H.upLogNC.tsv",header=F)
down<-read.delim(file="muscle.5tissues.H.downLogNC.tsv",header=F)
par(bty='n',mfrow=c(1,2))
boxplot(up[,-1],outline=F,ylim=c(-5,5),ylab="Exon/up-intron NC log ratio",names = c("human","monkey","treeShrew","mouse","pig"))
boxplot(down[,-1],outline=F,ylim=c(-5,5),ylab="Exon/down-intron NC log ratio",names = c("human","monkey","treeShrew","mouse","pig"))
#p values-up:ns,ns,0.0003866, 0.04656
#p values-down: ns,ns,0.001616,0.09677
####HR--- exons
#p values: uplogGC(NS,0.001,0.003,7.841e-06);downlogGC(ns,ns,0.03,0.06)
###All exons in diff ages
setwd("D:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/exonAge/final_list/NCratio")
h<-read.delim(file="brain.5species.uplogNC.final.tsv",header=F)
hr<-read.delim(file="brain.5tissues.HR.upLogNC.tsv",header=F)
others<-read.delim(file="brain.5tissues.HRT.upLogNC.tsv",header=F)
boxplot(h$V2,hr$V3,others[others$V2=="HRTM-",3],others[others$V2=="HRTMS",3],outline = F,names=c("H","HR","HRTM","HRTMS"),ylim=c(-5,6),ylab="Exon/up-intron logNC ratio")
abline(h=0,col="red")
boxplot(others[others$V2=="HRTM-",3],others[others$V2=="HRTM-",4],others[others$V2=="HRTM-",5],others[others$V2=="HRTM-",6],others[others$V2=="HRTM-",7],names = c("human","monkey","treeShrew","mouse","pig"),ylab="Exon/up-intron logNC ratio",outline = F,ylim=c(-6,8))
boxplot(others[others$V2=="HRTMS",3],others[others$V2=="HRTMS",4],others[others$V2=="HRTMS",5],others[others$V2=="HRTMS",6],others[others$V2=="HRTMS",7],names = c("human","monkey","treeShrew","mouse","pig"),ylab="Exon/up-intron logNC ratio",ylim=c(-6,8),outline = F)
setwd("../../burge/")
data<-read.delim(file="brain_logNC.tsv",header = F)
boxplot(V8~V7,data,ylab="Exon/up-intron logNC ratio",outline=F,ylim=c(-4,6))
#DAF in rhesus
setwd("./final_list/DAF/")
exon<-read.delim(file="H----.rheMac2.exon.snps.GC-AT.daf.final",header=F)
upExon<-read.delim(file="H----.rheMac2.exonUp150.snps.GC-AT.daf.final",header=F)
qqplot(exon$V8,upExon$V8,xlab="Exon",ylab="upExon",type="l",col="red")
abline(a=0,b=1,col="gray")
exon2<-read.delim(file="H----.rheMac2.exon.snps.AT-GC.daf",header=F)
upExon2<-read.delim(file="H----.rheMac2.exonUp150.snps.AT-GC.daf",header=F)
plot<-qqplot(exon2$V8,upExon2$V8,plot.it = F)
lines(plot$x,plot$y,col="blue")
plot(ecdf(exon$V8),do.points=F,verticals=TRUE,col="red")
lines(ecdf(upExon$V8),do.points=F,verticals=TRUE,col="blue")


#Cross-tissues comparison
setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mouse_tissue_data/comparison/exon_intron/final_AS")
diff<-read.delim(file="AS.differentialGC.psi.logNC.summary.tsv",header=F)
level<-read.delim(file="AS.leveledGC.psi.logNC.summary.tsv",header=F)
cor_d=''
cor_l=''
for(i in 1:nrow(diff)){cor_d[i]=cor(as.vector(as.numeric(diff[i,c(7:11)])),as.vector(as.numeric(diff[i,c(12:16)])))}
cor_d<-na.omit(as.numeric(cor_d))
for(i in 1:nrow(level)){cor_l[i]=cor(as.vector(as.numeric(level[i,c(7:11)])),as.vector(as.numeric(level[i,c(12:16)])))}
cor_l<-na.omit(as.numeric(cor_l))
plot(ecdf(cor_d),col="blue",do.points=F,verticals=TRUE,xlim=c(-1,1),main="",ylab="Cumlative fraction",xlab="Pearson correlation of psi and exon/upIntron logNC")
lines(ecdf(cor_l),col="red",do.points=F,verticals=TRUE)
abline(a=0.5,b=0.5,col="black",lty=1)
legend("topleft",c("Differential GC exons","Leveled GC exons","Uniform distribution"),lty=1,col=c("blue","red","black"))
ks.test(cor_d,cor_l)#p=0.01
boxplot(cor_d,cor_l,border =c("blue","red"),names=c("Differential GC exons","Leveled GC exons"),ylab="Pearson correlation of psi and exon/upIntron logNC")
wilcox.test(cor_d,cor_l,alternative = "g")#p=0.049
setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/monkey_tissue_data/comparison/exon_intron")

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/evaluation/")
data<-read.delim(file="insertSize.all.tsv",header=T,check.names = F)
freq<-cbind(data$size,data[,2:22]/t(replicate(151,colSums(data[,2:22]))))
colnames(freq)[1]<-"size"
data<-melt(freq,id.vars = 1)
ggplot(data,aes(x=size,y=value,color=variable))+geom_line(lwd=1)+labs(x="Fragment size(bp)",y="Frequency")
data<-read.delim(file="TSS.profile.all.tsv",header=T,check.names = F)
data<-melt(data,id.vars = 1)
ggplot(data,aes(x=size,y=value,color=variable))+geom_line(lwd=1)+labs(x="TSS relative position(bp)",y="Normalized nucleosome occupancy")
data<-read.delim(file="AT_freq.tsv",header = T,check.names = F)
data<-melt(data,id.vars = 1)
ggplot(data,aes(x=size,y=value,color=variable))+geom_line(lwd=1)+labs(x="Dyad relative position(bp)",y="Dinucleotide frequency")+ylim(0.025,0.125)+
  theme(panel.grid.minor.x=element_line(size=0.5,color="white"))+scale_x_continuous(minor_breaks = seq(-1000, 1000, 10),breaks=c(-73,0,73),labels=c("-73","0","73"))
data<-read.delim(file="GC_freq.tsv",header=T,check.names = F)
data<-melt(data,id.vars = 1)
ggplot(data,aes(x=size,y=value,color=variable))+geom_line(lwd=1)+labs(x="Dyad relative position(bp)",y="Dinucleotide frequency")+
  theme(panel.grid.minor.x=element_line(size=0.5,color="white"))+scale_x_continuous(minor_breaks = seq(-1000, 1000, 10),breaks=c(-73,0,73),labels=c("-73","0","73"))

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/mutiSpecies/Exon_intron/")
data<-read.delim(file="exon.meta.agg.tsv",comment.char = "#",header=T,check.names = F)
normF<-as.matrix(data[,-1])
normF<-normF/t(replicate(nrow(normF),normF[1,]))
normF<-cbind(data[,1],normF)
data<-melt(data,id.vars = 1)
ggplot(data,aes(x=size,y=value,color=variable))+geom_line(lwd=1)+labs(x="Exon relative position(bp)",y="Normalized nucleosome occupancy")+
  scale_x_continuous(limits = c(-500,700),breaks = seq(-500,700,by=100))
colnames(normF)[1]<-"size"
data<-melt(as.data.frame(normF),id.vars = 1)
normMean<-as.matrix(data[,-1])
colM<-colMeans(normMean)
normMean<-normMean/t(replicate(nrow(normMean),colM))
normMean<-cbind(data[,1],normMean)
normMean<-melt(as.data.frame(normMean),id.vars = 1)
ggplot(normMean,aes(x=V1,y=value,color=variable))+geom_line(lwd=1)+labs(x="Exon relative position(bp)",y="Normalized nucleosome occupancy")+
  scale_x_continuous(limits = c(-500,700),breaks = seq(-500,700,by=100))

setwd("./intra-sample/")
human<-read.delim(file="human.leveledGC.usage.uplogNC.tsv",header=F)
human<-human[,7:10]
colnames(human)[3:4]<-c("human-brain","human-muscle")
monkey<-read.delim(file="monkey.leveledGC.usage.uplogNC.tsv",header=F)
monkey<-monkey[,7:ncol(monkey)]
colnames(monkey)[6:10]<-c("monkey-brain","monkey-heart","monkey-kidney","monkey-liver","monkey-muscle")
tree<-read.delim(file="treeshrew.leveledGC.usage.uplogNC.tsv",header=F)
tree<-tree[,7:ncol(tree)]
colnames(tree)[5:8]<-c("treeshrew-brain","treeshrew-heart","treeshrew-kidney","treeshrew-muscle")
pig<-read.delim(file="pig.leveledGC.usage.uplogNC.tsv",header=F)
pig<-pig[,7:ncol(pig)]
colnames(pig)[6:10]<-c("pig-brain","pig-heart","pig-kidney","pig-liver","pig-muscle")
mouse<-read.delim(file="mouse.leveledGC.usage.uplogNC.tsv",header=F)
mouse<-mouse[,7:ncol(mouse)]
colnames(mouse)[6:10]<-c("mouse-brain","mouse-heart","mouse-kidney","mouse-liver","mouse-muscle")
data=matrix(ncol=4,nrow=0)
colnames(data)<-c("species","tissue","usage","log2NC")
myGetdata<-function(a,b){
  tmp_data<-matrix(ncol=4,nrow=0)
  colnames(tmp_data)<-c("species","tissue","usage","log2NC")
  for(i in 1:b){
    subset<-a[,c(i,i+b)]
    colnames(subset)[1]<-"usage"
    s_t<-strsplit(colnames(subset)[2],"-")[[1]]
    subset<-subset[subset$usage>=0,]
    subset[subset$usage>0.5,1]="H" #avoid "L" defined as "H"
    subset[subset$usage<0.5,1]="L"
    subset<-cbind(matrix(s_t,nrow=1,byrow = T),subset)
    colnames(subset)<-c("species","tissue","usage","log2NC")
    tmp_data<-rbind(tmp_data,subset)
  }
  return(tmp_data)
  
}
data<-rbind(data,myGetdata(human,2))
data<-rbind(data,myGetdata(monkey,5))
data<-rbind(data,myGetdata(tree,4))
data<-rbind(data,myGetdata(pig,5))
data<-rbind(data,myGetdata(mouse,5))
data<-data[data$usage=="H" | data$usage=="L",]
data$species<-factor(data$species,levels=c("human","monkey","treeshrew","pig","mouse"))
data$tissue<-as.factor(data$tissue)
data$usage<-as.factor(data$usage)
ggplot(data,aes(x=tissue,y=log2NC))+geom_boxplot(aes(fill=usage),position=position_dodge(.9),outlier.shape = NA)+facet_grid(.~species, scales="free", space="free")+ylim(-10,10)

species<-unique(data$species)
tissue<-unique(data$tissue)
pvalue=0
tag=0
for(i in 1:5){
  for(j in 1:5){
    tag=tag+1
    try({pvalue[tag]=wilcox.test(data[data$species==species[i] & data$tissue==tissue[j] & data$usage=="H",4],data[data$species==species[i] & data$tissue==tissue[j] & data$usage=="L",4],alternative = "g")$p.value})
  }
}
#[1] 1.981117e-02 1.442101e-02           NA           NA           NA
#[6] 1.249119e-05 1.946916e-01 5.326317e-01 8.045727e-03 5.081492e-01
#[11] 1.148467e-01 1.163676e-01 8.786367e-01 2.831382e-01           NA
#[16] 9.795253e-27 6.294362e-09 9.049568e-08 1.082686e-36 5.198643e-18
#[21] 5.089161e-03 3.974735e-01 3.168473e-04 3.923903e-02 3.023383e-01
setwd("./logGC_logNC/")
files=list.files(pattern = ".tsv$")
data<-matrix(ncol=3,nrow=0)
colnames(data)<-c("GCdiff","logNC","sample")
for(i in 1:length(files)){
  test<-read.delim(file=files[i],header=F)
  name=sub(".GCdiff.logNC.tsv","",files[i])
  test<-cbind(aggregate(V3~V4,data=test[is.finite(test$V3),],FUN="mean"),rep(name,7))
  colnames(test)<-c("GCdiff","logNC","sample")
  data<-rbind(data,test)
}
data$GCdiff<-factor(data$GCdiff,levels = c("<=-20","-20--10","-10-0","0-10","10-20","20-30",">=30"))
ggplot(data, aes(x=GCdiff, y=logNC, color=sample))+geom_point(shape=2)+geom_line(aes(group=sample),lwd=1)+labs(x="Exon/up-intron GC difference",y="log2(Exon/up-intron NC)")

setwd("E:/百度云同步盘/my PC D/projects/exon_evolution/summaryFigures/v2")
data<-read.delim(file="胶图quantity.forR.txt",header=T,check.names = F)
size1<-data[,c(1,2,3)]
size2<-na.omit(data[,c(4,5)])
size3<-na.omit(data[,c(6,7,8,9)])
size4<-na.omit(data[,c(10,11)])
size5<-na.omit(data[,12:16])
size6<-na.omit(data[,c(17,18)])
size7<-na.omit(data[,c(19,20)])
size8<-na.omit(data[,c(21,22)])
size9<-na.omit(data[,c(23,24,25)])
size10<-na.omit(data[,26:31])
colors<-rainbow(21)
myLine<-function(a,b){
  a<-as.matrix(a)
  colNum<-ncol(a)
  if(colNum>2){
    a[,2:colNum]<-a[,2:colNum]/matrix(colSums(a[,2:colNum]),nrow=nrow(a),ncol=colNum-1,byrow = TRUE)
    for(i in 2:colNum){
      lines(a[,1],a[,i],col=b[i-1],type="l",lwd=2)
    }
  }else{
    a[,2]<-a[,2]/sum(a[,2])
    lines(a[,1],a[,2],col=b[1],type="l",lwd=2)
  }
}
plot(data[,1],data[,2],type="n",xlim=c(0,1000),ylim=c(0,0.04),xlab="Size of migrating DNA population(bp)",ylab="Relative fluorescence ratio")
myLine(size1,colors[c(1,2)])
myLine(size2,colors[3])
myLine(size3,colors[4:6])
myLine(size4,colors[7])
myLine(size5,colors[8:11])
myLine(size6,colors[11])
myLine(size7,colors[13])
myLine(size8,colors[14])
myLine(size9,colors[15:16])
myLine(size10,colors[17:21])
legend<-colnames(data)[c(2,3,5,7:9,11,13:16,18,20,22,24,25,27:31)]
legend("topright",legend=legend,col=colors,lty=1,lwd=2)
