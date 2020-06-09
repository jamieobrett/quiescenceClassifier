library(TeachingDemos); library(MetabolAnalyze); library(ade4); library(wordcloud); library(flexclust)

#Read in data
input<-read.table("data.txt"); col<-read.table("colors.txt",header=T,as.is=T); #data.txt: cells are rows, genes are columns; colors.txt: "Class" column and "Color" column

#Set up colors
outlinecolor<-col$Color
rg<-col2rgb(outlinecolor)/255
fillcolor<-c()
for (i in 1:dim(rg)[2]) {
  fillcolor<-c(fillcolor,rgb(rg[1,i],rg[2,i],rg[3,i],alpha=0.1))
}
ellipsecolor=unique(col$Color);

#Calculate components
data <- t(na.omit(t(input)))
pcavar <- prcomp(data,center=T,scale=T) #centered genes, UV-scaled genes

#Plot cells
#Set up axes
ticks=6
xmin=min(pcavar$x[,'PC1']); xmax=max(pcavar$x[,'PC1']); ymin=min(pcavar$x[,'PC2']); ymax=max(pcavar$x[,'PC2'])
xran=ceiling((ceiling(xmax)-floor(xmin))/ticks)*ticks
yran=ceiling((ceiling(ymax)-floor(ymin))/ticks)*ticks
deltaxmin=ceiling((xran-(xmax-xmin))/2)
underxmin=trunc(xmin-deltaxmin)
overxmax=trunc(xmin+xran)
deltaymin=ceiling((yran-(ymax-ymin))/2)
underymin=trunc(ymin-deltaymin)
overymax=trunc(ymin+yran)
xincr=xran/ticks
yincr=yran/ticks
#Plot
png(filename="pca_cells.png",width=6,height=6,res=1200,units="in")
par(mar=c(4.5,5,0.1,0.1)) #bottom, left, top, right
plot.new()
par(mgp=c(3,2,0)) #axis.title.position,axis.label.position,axis.line.position
plot.window(xlim=c(underxmin,overxmax),ylim=c(underymin,overymax))
axis(1,at=seq(underxmin,overxmax,xincr),pos=underymin,lwd=2,cex.axis=2)
axis(2,at=seq(underymin,overymax,yincr),pos=underxmin,lwd=2,cex.axis=2)
s.class(pcavar$x,fac=factor(col$Class,levels=unique(col$Class)),add.plot=TRUE,col=ellipsecolor,axesell=F,label="",cpoint=0,cellipse=1,cstar=0) #standard deviation ellipse (centered on mean, orthogonal axes with the first along the great dispersion, each of length SD)
my.symbols(pcavar$x[,"PC1"],pcavar$x[,"PC2"],ms.filled.polygon,n=3,r=0.08,adj=seq(from=0,to=1,by=0.33),fg=outlinecolor,bg=fillcolor)
xtitle<-paste("PC1 Score (",round(100*summary(pcavar)$importance["Proportion of Variance",][1]),"% of Variance)",sep="")
ytitle<-paste("PC2 Score (",round(100*summary(pcavar)$importance["Proportion of Variance",][2]),"% of Variance)",sep="")
title(xlab=xtitle,ylab=ytitle,font.lab=2,cex.lab=2)
box(which="inner",col="white",lwd=6) #to remove the border placed by s.class
#par(oma=c(1.5,0,0,0)); mtext(format(Sys.time(),"%Y-%m-%d %H:%M"),cex=0.75,line=0,side=1,adj=0,outer=TRUE)
dev.off()

#Plot genes
#K-Medians
mydata<-scale(pcavar$rotation[,c(1,2)],center=T,scale=F)
fit<-kcca(mydata,k=2,family=kccaFamily("kmedians"),control=list(initcent="kmeanspp"),save.data=T); clusters(fit);
kcol<-vector(length=length(clusters(fit)))
mycolors<-c("blue","red")
rg<-col2rgb(mycolors)/255
kcol[clusters(fit)==1]=mycolors[1]
kcol[clusters(fit)==2]=mycolors[2]
ellcolor<-c()
for (i in 1:dim(rg)[2]) {
  ellcolor<-c(ellcolor,rgb(rg[1,i],rg[2,i],rg[3,i],alpha=0.4))
}
#Set up axes
ticks=4
xmin=min(pcavar$rotation[,'PC1']); xmax=max(pcavar$rotation[,'PC1']); ymin=min(pcavar$rotation[,'PC2']); ymax=max(pcavar$rotation[,'PC2'])
xran=(xmax-xmin)*1.2
yran=(ymax-ymin)*1.2
deltaxmin=(xran-(xmax-xmin))/2
underxmin=xmin-deltaxmin
overxmax=xmin+xran
deltaymin=(yran-(ymax-ymin))/2
underymin=ymin-deltaymin
overymax=ymin+yran
xincr=xran/ticks
yincr=yran/ticks
#Plot
png(filename="pca_genes.png",width=6,height=6,res=1200,units="in")
par(mar=c(4.5,5,0.1,0.1)) #bottom, left, top, right
plot.new()
par(mgp=c(3,2,0)) #axis.title.position,axis.label.position,axis.line.position
plot.window(xlim=c(underxmin,overxmax),ylim=c(underymin,overymax))
axis(1,at=round(seq(underxmin,overxmax,xincr),digits=2),pos=round(underymin,digits=2),lwd=2,cex.axis=1.5)
axis(2,at=round(seq(underymin,overymax,yincr),digits=2),pos=round(underxmin,digits=2),lwd=2,cex.axis=1.5)
s.class(pcavar$rotation,fac=as.factor(clusters(fit)),add.plot=TRUE,col=ellcolor,axesell=F,label="",cpoint=0,cellipse=1,cstar=0) #standard deviation ellipse (centered on mean, orthogonal axes with the first along the great dispersion, each of length SD)
textplot(pcavar$rotation[,1],pcavar$rotation[,2],row.names(pcavar$rotation),show.lines=T,new=F,col=kcol,font=2)
xtitle<-"PC1 Loading"
ytitle<-"PC2 Loading"
title(xlab=xtitle,ylab=ytitle,font.lab=2,cex.lab=2)
box(which="inner",col="white",lwd=6) #to remove the border placed by s.class
dev.off()
