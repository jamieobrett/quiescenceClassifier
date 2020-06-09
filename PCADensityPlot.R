library(MASS); library(RColorBrewer); library(wordcloud); library(TeachingDemos)

#Read in data
input<-read.table("data.txt"); col<-read.table("colors.txt",header=T,as.is=T); #data.txt: cells are rows, genes are columns; colors.txt: "Class" column and "Color" column

#Calculate components
data <- t(na.omit(t(input)))
pcavar <- prcomp(data,center=T,scale=T) #centered genes, UV-scaled genes

#Separate out classes
DPI0 <- pcavar$x[col$Class=="0dpi",c("PC1","PC2")]
DPI3 <- pcavar$x[col$Class=="3.5dpi",c("PC1","PC2")]
QM <- pcavar$x[col$Class=="QM_plastic",c("PC1","PC2")]
GM <- pcavar$x[col$Class=="GM_plastic",c("PC1","PC2")]
AMF <- pcavar$x[col$Class=="QM_AMF",c("PC1","PC2")]

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

#Plot setup
png(filename="pca_density.png",width=6,height=6,res=1200,units="in")
par(mar=c(4.5,5,0.1,0.1)) #bottom, left, top, right
plot.new()
par(mgp=c(3,2,0)) #axis.title.position,axis.label.position,axis.line.position
plot.window(xlim=c(underxmin,overxmax),ylim=c(underymin,overymax))
axis(1,at=seq(underxmin,overxmax,xincr),pos=underymin,lwd=2,cex.axis=2)
axis(2,at=seq(underymin,overymax,yincr),pos=underxmin,lwd=2,cex.axis=2)
xtitle<-paste("PC1 Score (",round(100*summary(pcavar)$importance["Proportion of Variance",][1]),"% of Variance)",sep="")
ytitle<-paste("PC2 Score (",round(100*summary(pcavar)$importance["Proportion of Variance",][2]),"% of Variance)",sep="")
title(xlab=xtitle,ylab=ytitle,font.lab=2,cex.lab=2)

#jamieDensity Function
jamieDensity <- function(pcares,samplecolor,n=101,colres=10,myunderxmin,myoverxmax,myunderymin,myoverymax) {
  mydensity <- kde2d(pcares[,"PC1"],pcares[,"PC2"],n=n,lims=c(myunderxmin,myoverxmax,myunderymin,myoverymax))

  #Colors
  ramp<-colorRampPalette(c("white",samplecolor)); rg<-col2rgb(ramp(colres))/255
  colors<-c()
  for (i in 1:dim(rg)[2]) {
    colors<-c(colors,rgb(rg[1,i],rg[2,i],rg[3,i],alpha=0.6))
  }

  #Draw density
  xinc<-(myoverxmax-myunderxmin)/(n-1)
  yinc<-(myoverymax-myunderymin)/(n-1)
  xstart<-myunderxmin+xinc/2
  ystart<-myunderymin+yinc/2
  zmax<-max(mydensity$z)

  for (x in 0:(length(mydensity$x)-1)) {
    for (y in 0:(length(mydensity$y)-1)) {
      density<-ceiling(colres*mydensity$z[x+1,y+1]/zmax)
      if (density > 1) {
        color<-colors[density]
        rect(xstart+x*xinc,ystart+y*yinc,xstart+(x+1)*xinc,ystart+(y+1)*yinc,col=color,border=NA)
      }
    }
  }

  #Draw contour
  contour(mydensity,add=T,axes=F,nlevels=5,drawlabels=F,col=samplecolor)

}

#Add densities to plot
jamieDensity(DPI3,"red",n=201,colres=10,underxmin,overxmax,underymin,overymax)
jamieDensity(DPI0,"blue",n=201,colres=10,underxmin,overxmax,underymin,overymax)

#jamiePCA Function
jamiePCA <- function(pcares,outlinecolor) {
  rg<-col2rgb(outlinecolor)/255
  fillcolor <- rgb(rg[1],rg[2],rg[3],alpha=0.4)
  my.symbols(pcares[,"PC1"],pcares[,"PC2"],ms.filled.polygon,n=3,r=0.06,adj=seq(from=0,to=1,by=0.33),fg=outlinecolor,bg=fillcolor)
}

#Add cells to plot
jamiePCA(QM,"brown")
jamiePCA(GM,"magenta")
jamiePCA(AMF,"green3")
dev.off()
