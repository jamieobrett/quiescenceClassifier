library(pvclust); library(ape)

#Read in data
input<-read.table("data.txt",header=T,as.is=T,sep="\t",quote="") #data.txt: first column = "Class", second column = "Color", third column = "Cells", fourth column = "Label", rest = gene names and expression values
#For subscripted numbers, flank with "z" and then "Z", such as Az10Z
data<-input[,5:dim(input)[2]]; rownames(data)=input$Cell; dcol<-input$Color; dlabels<-as.vector(input$Label); names(dcol)<-dlabels;

#Cluster data and obtain confidence levels
datascaled<-scale(data) #mean-centered, UV-scaled genes
result<-pvclust(t(datascaled),method.hclust="ward",method.dist="correlation",nboot=10000)
hr<-result$hclust
hr$labels<-dlabels
hrd1 <- as.dendrogram(hr,0.1)

colEdge<-function(n) {
  a<-attributes(n)
  mycol<-"mediumorchid2"
  attr(n,"edgePar") <- c(a$edgePar,list(col=mycol,lwd=thick))
  n
}
dendr1=hrd1
dendr1[[1]][[1]]=dendrapply(dendr1[[1]][[1]],colEdge)

#Make normal dendrogram
#Set up plot
png(filename="dendrogram.png",width=19,height=12,res=1200,units="in")
par(mar=c(0,2.5,0,0),oma=c(10,0,0,0))
plot(dendr1,type="triangle",edgePar=list(lwd=thick),axes=F,leaflab="none",edge.root=F) #or type="rectangle"
axis(side=2,at=seq(0,max(hr$height),max(hr$height)/ticks),lwd=3,font=2,labels=round(seq(0,max(hr$height),max(hr$height)/ticks)),line=0)
#Make labels with subscripts for compound concentrations
mye<-c()
for (i in 1:length(dlabels)) {
  s<-unlist(strsplit(labels(dendr1)[i],"z|Z"))

  if(length(s)==1) {
    e<-bquote(expression(.(s[1])))
  }
  else if(length(s)==2) {
    e<-bquote(expression(.(s[1])[.(s[2])]))
  }
  else if(length(s)==4) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]))
  }
  else if(length(s)==6) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]*.(s[5])[.(s[6])]))
  }
  else if(length(s)==8) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]*.(s[5])[.(s[6])]*.(s[7])[.(s[8])]))
  }
  else if(length(s)>8) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]*.(s[5])[.(s[6])]*.(s[7])[.(s[8])]*.(s[9])))
  }
  mye<-c(mye,e)
}

for(i in 1:length(dlabels)) {
  mtext(eval(mye[[i]]),side=1,line=0,outer=F,col=dcol[labels(dendr1)[i]],font=2,las=2,at=i)
}
dev.off()

#Make fan dendrogram
#Make labels with subscripts for compound concentrations
phr<-as.phylo(hr)
mylabels1<-c()
for (i in 1:length(dlabels)) {
  s<-unlist(strsplit(dlabels[i],"z|Z"))
  if(length(s)==1) {
    e<-bquote(expression(.(s[1])))
  }
  else if(length(s)==2) {
    e<-bquote(expression(.(s[1])[.(s[2])]))
  }
  else if(length(s)==4) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]))
  }
  else if(length(s)==6) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]*.(s[5])[.(s[6])]))
  }
  else if(length(s)==8) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]*.(s[5])[.(s[6])]*.(s[7])[.(s[8])]))
  }
  else if(length(s)>8) {
    e<-bquote(expression(.(s[1])[.(s[2])]*.(s[3])[.(s[4])]*.(s[5])[.(s[6])]*.(s[7])[.(s[8])]*.(s[9])))
  }
  mylabels1<-c(mylabels1,eval(e))
}
phr$tip.label<-mylabels1
#Make plot
png(filename="fan.png",width=8,height=8,res=1200,units="in")
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(phr,type="fan",no.margin=T,cex=0.8,label.offset=0.1,edge.color=c("black",rep("mediumorchid2",15),rep("black",188)),tip.color=dcol,font=2,edge.width=1.5)
dev.off()
