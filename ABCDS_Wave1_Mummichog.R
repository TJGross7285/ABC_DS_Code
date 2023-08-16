library(dplyr)
library(impute)
library(sva)
library(limma)
library(EnhancedVolcano)

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

pdf(file="UPDATED_ABCDS_Volcano Plots_12_16.pdf")
with(final_DE_matrix_T, plot(MCI.Control, -log10(adj.P.Val), pch=20, main="DS-MCI vs. CS-NAD",xlab="Log2 Fold Change",ylab="-log10 Adjusted P",xlim=c(-2,2), ylim=c(0,2)))
with(subset(final_DE_matrix_T, adj.P.Val<.05 ), points(MCI.Control, -log10(adj.P.Val), pch=20, col="red"))
with(subset(final_DE_matrix_T, adj.P.Val>.05 ), points(MCI.Control, -log10(adj.P.Val), pch=20, col="darkgreen"))
abline(a=1.301,b=0,lty="dashed",col="firebrick4")



with(final_DE_matrix_U, plot(Dementia.MCI, -log10(adj.P.Val), pch=20, main="DS-AD vs. DS-MCI", xlab="Log2 Fold Change",ylab="-log10 Adjusted P",xlim=c(-2,2), ylim=c(0,2.5)))
with(subset(final_DE_matrix_U, adj.P.Val<.05 ), points(Dementia.MCI, -log10(adj.P.Val), pch=20, col="red"))
with(subset(final_DE_matrix_U, adj.P.Val>.05 ), points(Dementia.MCI, -log10(adj.P.Val), pch=20, col="darkgreen"))
abline(a=1.301,b=0,lty="dashed",col="firebrick4")



with(final_DE_matrix_V, plot(Dementia.Control, -log10(adj.P.Val), pch=20, main="DS-AD vs. CS-NAD", xlab="Log2 Fold Change",ylab="-log10 Adjusted P",xlim=c(-2,2), ylim=c(0,4)))
with(subset(final_DE_matrix_V, adj.P.Val<.05 ), points(Dementia.Control, -log10(adj.P.Val), pch=20, col="red"))
with(subset(final_DE_matrix_V, adj.P.Val>.05 ), points(Dementia.Control, -log10(adj.P.Val), pch=20, col="darkgreen"))
abline(a=1.301,b=0,lty="dashed",col="firebrick4")

dev.off()












pdf(file="ABCDS_UNTARGETED_FDR_Volcano_12_16.pdf")
EnhancedVolcano(final_DE_matrix_T,
	lab=NULL,
	title="DS-MCI vs. CS-NAD",
    x = "MCI.Control",
    y = 'adj.P.Val',
    xlim=c(-2,2),
    ylim=c(0,2),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    pCutoff = 0.05,
    FCcutoff = 0,
    labSize = 4.0,
    colAlpha = 1,
    legend=c('NS','Log2 FC','Adjusted p-value',
      'Adjusted p-value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0)
















EnhancedVolcano(final_DE_matrix_U,
	lab=NULL,
	title="DS-AD vs. DS-MCI",
    x = "Dementia.MCI",
    y = 'adj.P.Val',
    xlim=c(-2,2),
    ylim=c(0,2.5),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    pCutoff = 0.05,
    FCcutoff= 0,
    labSize = 4.0,
    colAlpha = 1,
    legend=c('NS','Log2 FC','Adjusted p-value',
      'Adjusted p-value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0)

EnhancedVolcano(final_DE_matrix_V,
	lab=NULL,
	title="DS-AD vs. CS-NAD",
    x = "Dementia.Control",
    y = 'adj.P.Val',
    xlim=c(-2,2),
    ylim=c(0,4),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    pCutoff = 0.05,
    FCcutoff= 0,
    labSize = 4.0,
    colAlpha = 1,
    legend=c('NS','Log2 FC','Adjusted p-value',
      'Adjusted p-value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0)

dev.off()










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