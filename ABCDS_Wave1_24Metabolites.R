T<-read.csv("MERGED_UPDATED_ABCDS_10_15_19_CHECK.csv",na.strings="")%>%filter(DiseaseState!="UNKNOWN")
H<-read.csv("ABCDS_24_Data.csv",header=FALSE)[-c(1:6),]
colnames(H)<-unlist(read.csv("ABCDS_24_Data.csv",header=FALSE)[6,])

posed<-t(H)[-c(1),]
colnames(posed)<-t(H)[1,]
final<-cbind(as.data.frame(rownames(posed)),as.data.frame(posed))
colnames(final)[1]<-"GeorgetownID"
rownames(final)<-NULL

R<-dplyr::inner_join(T,final,by="GeorgetownID")
write.csv(R,file="ABCDS_24Metabs_Clinical_9_2_2020.csv")


read<-read.csv("ABCDS_24Metabs_Clinical_9_2_2020.csv",check.names=FALSE)

abunds<-log2(read[is.na(read$ANYE4)==FALSE,-c(1:40)])
pheno<-as.factor(read[is.na(read$ANYE4)==FALSE,]$Consensus_DX)
levels(pheno)<-list("NC/MCI"=c("0-No dementia, no MCI","1-MCI"),"AD"=c("2-DEMENTIA"))
pheno<-cbind(pheno,read[is.na(read$ANYE4)==FALSE,c(8,10,13)])
colnames(pheno)<-c("DiseaseState","Age","Sex","ANYE4")
clean<-empiricalBayesLM(abunds, removedCovariates=pheno[,-c(1)])$adjustedData


cctrl1<- trainControl(method = "repeatedcv", repeats = 10, number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
set.seed(122)
levels(pheno$DiseaseState)<-c("Control","Case")
BoostedLogistic<-train(clean,pheno$DiseaseState, method="rf", trControl = cctrl1, metric = "ROC")

pdf(file="plots_24Features.pdf")
for (i in 1:dim(abunds)[2]){
	vioplot(abunds[,i]~as.factor(read$Consensus_DX),main=colnames(abunds)[i])
}
dev.off()

remove<-as.data.frame(read[,c(8,10)])
K<-empiricalBayesLM(abunds, removedCovariates=remove)$adjustedData

pdf(file="plots_24Features_RESIDUALIZED.pdf")
for (i in 1:dim(K)[2]){
	vioplot(K[,i]~as.factor(read$Consensus_DX),main=colnames(K)[i])
}
dev.off()




