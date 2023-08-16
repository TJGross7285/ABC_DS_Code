library(dplyr)
library(impute)
library(sva)
library(limma)

pos<-read.csv("mapstone_abcds_wave_1_metabolomics_POS.csv",header=FALSE,check.names=FALSE, na.strings="0")
neg<-read.csv("mapstone_abcds_wave_1_metabolomics_NEG.csv",header=FALSE,check.names=FALSE, na.strings="0")

pos_meta<-as.data.frame(pos[,1:2])
colnames(pos_meta)<-c("MZ","RT")
pos_meta<-cbind(pos_meta,rep("ESI+",dim(pos_meta)[1]))
colnames(pos_meta)[3]<-"Mode"
pos_meta<-pos_meta[-c(1),]
pos_abunds<-log2(apply(as.data.frame(t(pos[2:dim(pos)[1],3:dim(pos)[2]])),2,as.numeric))
SampleID<-as.data.frame(pos[1,3:dim(pos)[2]])
pos_combine<-cbind(t(SampleID),pos_abunds)
colnames(pos_combine)<-c("SampleID",paste0(seq(1,dim(pos_abunds)[2]),rep(".POS",dim(pos_abunds)[2])))
write.csv(pos_combine,file="ABCDS_POS.csv")
pos_meta<-cbind(colnames(pos_combine)[-c(1)],pos_meta)
colnames(pos_meta)[1]<-"Feature"

neg_meta<-as.data.frame(neg[,1:2])
colnames(neg_meta)<-c("MZ","RT")
neg_meta<-cbind(neg_meta,rep("ESI-",dim(neg_meta)[1]))
colnames(neg_meta)[3]<-"Mode"
neg_meta<-neg_meta[-c(1),]
neg_abunds<-log2(apply(as.data.frame(t(neg[2:dim(neg)[1],3:dim(neg)[2]])),2,as.numeric))
SampleID<-as.data.frame(neg[1,3:dim(neg)[2]])
neg_combine<-cbind(t(SampleID),neg_abunds)
colnames(neg_combine)<-c("SampleID",paste0(seq(1,dim(neg_abunds)[2]),rep(".NEG",dim(neg_abunds)[2])))
write.csv(neg_combine,file="ABCDS_NEG.csv")
neg_meta<-cbind(colnames(neg_combine)[-c(1)],neg_meta)
colnames(neg_meta)[1]<-"Feature"

total_meta<-rbind(pos_meta,neg_meta)
write.csv(total_meta,file="Meta.csv")





library(dplyr)
library(impute)
library(sva)
library(limma)
####DO NOT OVERWRITE!!!!!!!!!!
pos<-read.csv("ABCDS_POS.csv",check.names=FALSE)
neg<-read.csv("ABCDS_NEG.csv",check.names=FALSE)

all.equal(pos$SampleID,neg$SampleID)

combined<-cbind(pos,neg[,-c(1)])
colnames(combined)[1]<-"GeorgetownID"
T<-read.csv("MERGED_UPDATED_ABCDS_10_15_19_CHECK.csv",na.strings="")%>%filter(DiseaseState!="UNKNOWN")


combined<-dplyr::inner_join(T,combined)
combined<-combined[is.na(combined$ANYE4)==FALSE,]

pheno<-combined[,1:38]
pheno$Consensus_DX<-as.factor(factor(pheno$Consensus_DX))
edata<-t(combined[,-c(1:38)])

####Impute missing data with KNN imputation; n.sv=25
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data
mod<-model.matrix(~as.factor(Consensus_DX)+as.factor(ANYE4),data=pheno)
mod0<-model.matrix(~as.factor(ANYE4),data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv




total<-cbind(as.data.frame(pheno$Consensus_DX),as.data.frame(svobj))
colnames(total)<-c("Main","SV1","SV2","SV3","SV4",
					"SV5","SV6","SV7","SV8","SV9",
					"SV10","SV11","SV12","SV13","SV14","SV15",
					"SV16","SV17","SV18","SV19","SV20","SV21",
					"SV22","SV23","SV24","SV25")
####Set up accessory objects for DE incorporating pre-post timepoints
design1<-as.data.frame(model.matrix(~0+Main+SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14+SV15+SV16+SV17+SV18+SV19+SV20+SV21
					+SV22+SV23+SV24+SV25,data=total))
colnames(design1)[1:3]<-c("Control","MCI","Dementia")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm1 <- makeContrasts(
	`MCI-Control` = MCI-Control,
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

fit2<-lmFit(edata,design1)
cm2 <- makeContrasts(
	`Dementia-MCI`= Dementia-MCI,
	levels=design1)
fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)

fit3<-lmFit(edata,design1)
cm3 <- makeContrasts(
	`Dementia-Control`= Dementia-Control,
	levels=design1)
fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"


U<-topTableF(fit2_F,adjust="BH",number=100000)
U<-cbind(rownames(U),U)
colnames(U)[1]<-"Feature"

V<-topTableF(fit3_F,adjust="BH",number=100000)
V<-cbind(rownames(V),V)
colnames(V)[1]<-"Feature"

total_meta<-read.csv("Meta.csv")

final_DE_matrix_T<-dplyr::inner_join(total_meta,T,by="Feature")
final_DE_matrix_U<-dplyr::inner_join(total_meta,U,by="Feature")
final_DE_matrix_V<-dplyr::inner_join(total_meta,V,by="Feature")



R<-WGCNA::empiricalBayesLM(as.data.frame(t(edata)),removedCovariates=svobj)$adjustedData
total<-cbind(pheno[,c(1,5)],as.data.frame(R))%>%dplyr::filter(Consensus_DX!="0-No dementia, no MCI")
index_table<-final_DE_matrix_U%>%filter(adj.P.Val<.05)%>%dplyr::select(Feature)
index<-colnames(total) %in% index_table$Feature
final_MCI_AD_table<-cbind(total[,1:2],total[,index])

write.csv(final_MCI_AD_table,file="SVA_Adjusted_MCI_AD_12_17.csv")










write.csv(final_DE_matrix_T,file="ABCDS_MCI_Control_DE.csv")
write.csv(final_DE_matrix_U, file="ABCDS_AD_MCI_DE.csv")
write.csv(final_DE_matrix_V, file="ABCDS_AD_Control_DE.csv")

write.table(final_DE_matrix_T%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCI.Control"),
			sep="\t",file="DS_WAVE1_POS_MCI_Control.txt",row.names=FALSE)
write.table(final_DE_matrix_T%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCI.Control"),
			sep="\t",file="DS_WAVE1_NEG_MCI_Control.txt",row.names=FALSE)

write.table(final_DE_matrix_U%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Dementia.MCI"),
			sep="\t",file="DS_WAVE1_POS_AD_MCI.txt",row.names=FALSE)
write.table(final_DE_matrix_U%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Dementia.MCI"),
			sep="\t",file="DS_WAVE1_NEG_AD_MCI.txt",row.names=FALSE)

write.table(final_DE_matrix_V%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Dementia.Control"),
			sep="\t",file="DS_WAVE1_POS_AD_Control.txt",row.names=FALSE)
write.table(final_DE_matrix_V%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Dementia.Control"),
			sep="\t",file="DS_WAVE1_NEG_AD_Control.txt",row.names=FALSE)

cd "/Users/TGross/Desktop/ABCDS Pheno" 
source activate py27
mummichog -f DS_WAVE1_POS_MCI_Control.txt -o DS_WAVE1_POS_MCI_Control -m positive -u 7 
mummichog -f DS_WAVE1_NEG_MCI_Control.txt -o DS_WAVE1_NEG_MCI_Control -m negative -u 7 

mummichog -f DS_WAVE1_POS_AD_MCI.txt -o DS_WAVE1_POS_AD_MCI -m positive -u 7 
mummichog -f DS_WAVE1_NEG_AD_MCI.txt -o DS_WAVE1_NEG_AD_MCI -m negative -u 7

mummichog -f DS_WAVE1_POS_AD_Control.txt -o DS_WAVE1_POS_AD_Control -m positive -u 7 
mummichog -f DS_WAVE1_NEG_AD_Control.txt -o DS_WAVE1_NEG_AD_Control -m negative -u 7


mummichog -f PACEMCI_NEG.txt -o test1 -m negative -u 7


mummichog -f PACEMCI_POS.txt -o test2 -m positive -u 7
