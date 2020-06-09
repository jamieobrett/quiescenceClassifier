library(languageR); library(rms); library(party); library(e1071)

#Read in data
data<-read.table("trainingdata.txt",header=T) #trainingdata.txt: first column = cell names, second column = "Class", rest of columns = gene expression values
testdata<-read.table("testdata.txt",header=T) #formatted as trainingdata.txt

#Random forest construction
ngenes<-dim(data)[2]-2
data.controls <- cforest <- unbiased(ntree=8000, mtry=round(sqrt(ngenes))) #cannot overfit, so more trees is better (at least 1000, Metaboanalyst does 5000, recommended max is 8000)
set.seed(1234) #any random number, make sure to do a few times with different seeds
genes<-names(data[,-1])
fmla<-as.formula(paste("Class ~ ",paste(genes,collapse="+")))
class(data$Class); #make sure is "factor"
data.cforest <- cforest(fmla, data=data,controls=data.controls)

#Variable importance calculation
data.cforest.varimp <- varimp(data.cforest,conditional=TRUE,nperm=1,threshold=0.2,OOB=TRUE)
write.table(data.cforest.varimp, file="RFvarimp.csv")
png(filename="rf_varimp.png",width=6,height=7,res=1200,units="in")
par(cex.axis=1.5,font=2,lwd=2,cex.lab=1.5,mgp=c(2.5,1,0))
dotchart2(rev(sort(data.cforest.varimp)),xlab="Variable Importance",pch=16,lty=1,dotsize=1.2,cex.labels=1,lcolor="gray",col=c(rep("red",8),rep("black",17)),axisat=seq(0,0.03,0.01),axislabels=seq(0,0.03,0.01))
box()
dev.off()

#Prediction using random forest
write.csv(table(predict(data.cforest,OOB=FALSE),data$Class),file="All.csv") #OOB=FALSE to use all cells
write.csv(table(predict(data.cforest,OOB=TRUE),data$Class),file="OOB.csv") #OOB=TRUE to get the real OOB error: column categories are real classes, row categories are predicted classes
write.csv(table(Predict(data.cforest,newdata=testdata),testdata$Class),file="output.csv")
party:::prettytree(data.cforest@ensemble[[1]], names(data.cforest@data@get("input")))
#treeresponse(data.cforest,newdata=testdata) #to get probabilities instead of classifications
