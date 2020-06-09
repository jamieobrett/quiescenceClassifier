#Note this takes a very long time unless run on a high-power cluster
library(party)

#Read data
data<-read.table("data.txt",header=T); #data.txt: first column = cell names, second column = "Class", rest of columns = gene expression values

#Random forest construction
ngenes<-dim(data)[2]-2
data.controls <- cforest_unbiased(ntree=8000, mtry=round(sqrt(ngenes)))
set.seed(1234)
genes<-names(data[,-1])
fmla<-as.formula(paste("Class ~ ",paste(genes,collapse="+")))

#Read variable importances
input<-read.table("varimps.txt",header=T) #varimps.txt: first column = gene names, second column = variable importances from conditionalRF.R
data.cforest.varimp<-input[,1]
names(data.cforest.varimp)<-row.names(input)

#Permute
counts<-rep(0,length(data.cforest.varimp));  names(counts)<-names(data.cforest.varimp)
nrep<-1000
for (i in 1:nrep) {
  data$Class <- sample(data$Class)
  cf <- cforest(fmla,data=data,controls=data.controls); vi <- varimp(cf,conditional=FALSE,nperm=1,threshold=0.2,OOB=TRUE)
  counts[vi>=data.cforest.varimp] <- counts[vi>=data.cforest.varimp]+1
}

#Export
write.table(counts,file="counts.txt") #counts.txt: first column = gene name, second column = number of times this gene had a variable importance of at least that in varimps.txt
