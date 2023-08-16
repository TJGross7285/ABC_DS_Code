#################################################(10/19/19) ABCDS Lipidyzer DE  
library(impute)
library(dplyr)
T<-read.csv("MERGED_UPDATED_ABCDS_10_15_19_CHECK.csv",na.strings="")%>%filter(DiseaseState!="UNKNOWN")

####Import lipidyzer abundances
C<-read.csv("Lipid_Species_Concentration.csv",header=FALSE, na.strings=c(".","0"))
spec<-C[-c(1:5),-c(1:5)]
colnames(spec)<-unlist(C[4,6:dim(C)[2]])
D<-read.csv("Lipid_Class_Concentration.csv",header=FALSE, na.strings=c(".","0"))
class<-D[-c(1:5),-c(1:5)]
colnames(class)<-unlist(D[4,6:dim(D)[2]])
E<-read.csv("Fatty_Acid_Concentration.csv",header=FALSE, na.strings=c(".","0"))
FA<-E[-c(1:5),-c(1:5)]
colnames(FA)<-unlist(E[4,6:238])
F<-read.csv("Total_Fatty_Acid.csv",header=FALSE, na.strings=c(".","0"))
total<-F[-c(1:5),-c(1:5)]
colnames(total)<-unlist(F[4,6:32])

####Check that Georgetown IDs are same across tables 
identical(C[-c(1:5),1],D[-c(1:5),1],E[-c(1:5),1],F[-c(1:5),1])
GeorgetownID<-as.character(C[-c(1:5),1])

####Combine tables and join to metadata; threshold those with excessive NA 
bind<-cbind(as.data.frame(GeorgetownID),as.data.frame(spec),as.data.frame(class),as.data.frame(FA),as.data.frame(total))
index<-caret::nearZeroVar(bind)
bind<-bind[,-index]

combined<-dplyr::inner_join(T,bind,by="GeorgetownID")
na<-apply(combined,2,is.na)
na_index<-apply(na,2,sum)/dim(combined)[1]
combined<-combined[,na_index<.333]
edata<-log2(t(apply(combined[,-c(1:31)],2,as.numeric)))
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data


####Save data object for modeling 
###########save(edata_prime,combined,file="ABCDS_for_fsva.RData")

pheno<-combined[,c(1:31)]
pheno$Consensus_DX<-as.factor(factor(pheno$Consensus_DX))
index<-is.na(pheno$ANYE4)==FALSE
pheno<-pheno[index,]
edata<-edata[,index]

pdf(file="lognorm_ABCDS.pdf") 
for (i in 1:ncol(t(edata))){
	boxplot(t(edata)[,i]~ pheno$Consensus_DX, main=colnames(t(edata))[i], ylab="Log2 Abundance")
}
dev.off()

####Impute missing data with KNN imputation
#!#!#!#!#!#!#edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data
mod<-model.matrix(~as.factor(Consensus_DX)+as.factor(factor(ANYE4)),data=pheno)
mod0<-model.matrix(~as.factor(factor(ANYE4)),data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(as.data.frame(pheno$Consensus_DX),as.data.frame(svobj))
colnames(total)<-c("Main","SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10","SV11","SV12","SV13","SV14")

####Set up accessory objects for DE incorporating pre-post timepoints #SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14
design1<-as.data.frame(model.matrix(~0+Main,data=total))
colnames(design1)[1:3]<-c("Control","MCI","Dementia")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm <- makeContrasts(
	`MCI-Control` = MCI-Control,
	`Dementia-MCI`= Dementia-MCI,
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"

write.csv(T, file="DS_DE_ABCDS_noSV.csv")


library(EnhancedVolcano)
colnames(T)[2]<-"Log2FC"
pdf(file="WAVE1_ABCDS_Lipidyzer_Volcano_10_17_19.pdf")
EnhancedVolcano(T,lab=T$Feature,x="Log2FC",y="P.Value",xlim=c(-.5,.5),
	            ylim=c(0,10),FCcutoff=.2, pLabellingCutoff=.025)
dev.off()


















