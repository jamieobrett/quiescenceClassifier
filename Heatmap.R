#Custom heatmap script

#Read in data
anndata=read.table("anno.txt",header=T); #anno.txt: first column = condition name, column "A" = conc of A, column "B" = conc of B, column "C" = conc of C, column "D" = conc of D, column "QM" = Y or N
input<-read.table("data.txt",header=T,as.is=T,sep="\t",quote="") #data.txt: first column = "Class", second column = "Color", third column = "Cells", fourth column = "Label", rest = gene names and expression values
roworder<-read.table("rorder.txt") #rorder.txt: for each gene, which column to put it in the heatmap
colorder<-read.table("corder.txt") #corder.txt: for each cell, which row to put it in the heatmap
data<-input[,5:dim(input)[2]]; rownames(data)=input$Cell; dcol<-input$Color; dlabels<-as.vector(input$Label); names(dcol)<-dlabels; datascaled<-scale(data); forhm<-datascaled #mean-centered, UV-scaled genes

#Some parameters
colors<-c("#00FF00","#00F900","#00F300","#00ED00","#00E800","#00E200","#00DC00","#00D700","#00D100","#00CB00","#00C500","#00C000","#00BA00","#00B400","#00AF00","#00A700","#009B00","#008F00","#008300","#007700","#006B00","#005F00","#005300","#004700","#003B00","#002F00","#002300","#001700","#000B00","#000000","#0B0000","#170000","#230000","#2F0000","#3B0000","#470000","#530000","#5F0000","#6B0000","#770000","#830000","#8F0000","#9B0000","#A70000","#AF0000","#B40000","#BA0000","#C00000","#C50000","#CB0000","#D10000","#D70000","#DC0000","#E20000","#E80000","#ED0000","#F30000","#F90000","#FF0000")
lowcutoff<- -5; highcutoff<- 5; #Set extreme outlier values to these cutoffs for coloring purposes
forhm<-replace(forhm,forhm < lowcutoff,lowcutoff); forhm<-replace(forhm,forhm > highcutoff, highcutoff)

#Set up names
rown=c();
for (i in 1:length(roworder)) {
  rown<-c(rown,dlabels[roworder[i] ])
}
coln=c()
for (i in 1:length(colorder)) {
  myname<-names(data)[colorder[i] ]
  myname<-unlist(strsplit(myname,"_"))[1] #remove _a or _b from probesets
  coln<-c(coln,myname)
}

#Arrange and normalize data for graphing
R<-length(roworder); C<-length(colorder); A<-dim(anndata)[2]
m<-matrix(data=NA,nrow=R,ncol=C); rownames(m)<-rown; colnames(m)<-coln
an<-matrix(data=NA,nrow=R,ncol=A); rownames(an)<-rown
for (cell in 1:R) {
  for (gene in 1:C) {
    m[cell,gene]<-forhm[roworder[cell],colorder[gene]]
  }
  for (anno in 1:A) {
    an[which(roworder==cell),anno]<-anndata[cell,anno]
  }
}
datamax<-max(forhm); datamin<-min(forhm); datarange<-datamax-datamin
normalized<-(m-datamin)/datarange

#Heatmap
png(filename="hm.png",width=14,height=10,res=600,units="in")
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,C+9),c(-14,R),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
for(myr in 1:R) {
  for(myc in 1:C) {
    colindex=round(normalized[myr,myc]*(length(colors)-1))+1
    rect(myc-1,R-myr,myc,R-myr+1,col=colors[colindex],border=NA) #[1,1] starts at the top left
  }
}

text(seq(0.5,length(coln)-0.5,1),-0.15,labels=coln,cex=1,adj=c(0,0.5),srt=270)

#Supplment Conc
mygreys<-rev(gray.colors(100,end=0.95,start=0,gamma=0.5));
for(myr in 1:R) {
  for(mya in 1:(A-1)) {
    colindex=round(an[myr,mya]/max(an[,mya])*(length(mygreys)-1))+1
    rect(C+mya,R-myr,C+1+mya,R-myr+1,col=mygreys[colindex],border=NA) #[1,1] starts at the top left
  }
}

#QM label
for(myr in 1:R) {
  if(an[myr,A]==2) {
    qcol<-"blue"
  } else {
    qcol<-"white"
  }
  rect(C+A,R-myr,C+1+A,R-myr+1,col=qcol,border=NA) #[1,1] starts at the top left
}

#Days after isolation label
for(myr in 1:R) {
  if(rown[myr]=="0DPI") {
    days<-0
  } else if (rown[myr]=="Early_QM") {
    days<-1.25
  } else {
    days<-2.5
  }
  points(C+2+days+A,R-myr+0.5,col="black",pch=20)
}

text(seq(C+1.5,C+4.5,1),-0.22,labels=c("A","B","C","D"),cex=1,adj=c(0.5,1))
text(C+5.3,-0.47,labels=c("QM"),cex=1,adj=c(0,0.5),srt=315)
text(c(C+2+A,C+3.2+A,C+4.7+A),-0.22,labels=c("0","1.5","2.5"),cex=0.9,adj=c(0.5,1))
text(C+A+3.6,-2.6,labels=c("Days"),cex=1,adj=c(0.5,1))
text(C+A+3.6,-4.5,labels=c("In Vitro"),cex=1,adj=c(0.5,1))
rect(C+1,R,C+1+A,0,border="black")
rect(C+1.5+A,R,C+5+A,0,border="black")

#Highlight the QM cluster
for (myr in 1:R) {
  if(rown[myr]=="0DPI" | rown[myr]=="Early_QM" | rown[myr]=="Az100ZBz100ZCz100ZDz100ZQM") {
    rect(C,R-myr,C+0.6,R-myr+1,col="mediumorchid2",border=NA)
  }
}
dev.off()

#Legend
png(filename="legend.png",width=3,height=1,res=600,units="in")
par(mar=c(1.2,0.7,0,0.7),oma=c(0,0,0,0),mgp=c(0,0.3,0))
image(matrix(1:length(colors),length(colors),1),0.1,col=colors,ylab="",axes=F,xlim=c(0,1)); axis(1,at=(0:6 / 6),labels=sprintf("%.1f",seq(from=min(forhm),to=max(forhm),length.out=7)),cex.axis=0.8,font=2)
box()
dev.off()

png(filename="legend2.png",width=3,height=1,res=600,units="in")
par(mar=c(1.2,0.7,0,0.7),oma=c(0,0,0,0),mgp=c(0,0.3,0))
image(x=0:length(mygreys),y=0.1,matrix(1:length(mygreys),length(mygreys),1),col=mygreys,ylab="",axes=F,xlim=c(0,100),xlab=""); axis(1,at=c(0,10,20,50,100),labels=c(0,10,20,50,100),cex.axis=0.8,font=2)
box()
dev.off()
