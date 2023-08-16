library(dplyr)
library(impute)
library(sva)
library(limma)

####monocle 2 used here; use latest stable release from Bioconductor  
library(monocle)

####DO NOT OVERWRITE!!!!!!!!!!
pos<-read.csv("ABCDS_POS.csv",check.names=FALSE)
neg<-read.csv("ABCDS_NEG.csv",check.names=FALSE)

all.equal(pos$SampleID,neg$SampleID)

combined<-cbind(pos,neg[,-c(1)])
colnames(combined)[1]<-"GeorgetownID"
T<-read.csv("MERGED_UPDATED_ABCDS_10_15_19_CHECK.csv",na.strings="")%>%filter(DiseaseState!="UNKNOWN")


combined<-dplyr::inner_join(T,combined)
####combined<-combined[is.na(combined$ANYE4)==FALSE,]

pheno<-combined[,1:38]
pheno$Consensus_DX<-as.factor(factor(pheno$Consensus_DX))
edata<-t(combined[,-c(1:38)])

####Impute missing data with KNN imputation; n.sv=25
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data
mod<-model.matrix(~as.factor(Consensus_DX)+as.factor(ANYE4),data=pheno)
mod0<-model.matrix(~as.factor(ANYE4),data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv

####Set up gene-level metadata 
total_meta<-read.csv("Meta.csv", check.names=FALSE)
rownames(total_meta)<-total_meta$Feature
total_meta<-total_meta[,-c(1)]
colnames(total_meta)[1]<-"gene_short_name"

####Set up sample-level metadata 
SV_pheno<-cbind(pheno,svobj)
colnames(SV_pheno)[39:63]<-c("SV1","SV2","SV3","SV4",
					"SV5","SV6","SV7","SV8","SV9",
					"SV10","SV11","SV12","SV13","SV14","SV15",
					"SV16","SV17","SV18","SV19","SV20","SV21",
					"SV22","SV23","SV24","SV25")

####Initialize metadata objects and total expression set
fd <- new("AnnotatedDataFrame", data = total_meta)
pd <- new("AnnotatedDataFrame", data = pheno)

HSMM <- newCellDataSet(as.matrix(edata),
    phenoData = pd, featureData = fd, expressionFamily= uninormal())


AD_normal<-colnames(read.csv("AD_Control_SVA_Adjusted_Abunds_11_13.csv",check.names=FALSE))[-c(1:7)]
AD_mci<-colnames(read.csv("MCI_AD_SVA_Adjusted_Abunds_11_14.csv", check.names=FALSE))[-c(1:2)]
DE_genes<-c(AD_normal,AD_mci)

subsetted <- row.names(subset(fData(HSMM),
          gene_short_name %in% DE_genes))

####Carry out dimension reduction/sample ordering using SVs as residual medel formula   
HSMM_myo<- reduceDimension(HSMM, reduction_method = 'DDRTree',norm_method="none", pseudo_expr=0,
                            residualModelFormulaStr="~as.numeric(age)+as.factor(Gender)")
HSMM_myo1 <- orderCells(HSMM_myo)

pdf(file="ABCDS_Monocle_Heatmap_FDRsigONLY.pdf")
A<-plot_pseudotime_heatmap(HSMM_myo1[subsetted,],
                num_clusters = 10,
                cores = 60,
                show_rownames = FALSE)
dev.off()

save(A,HSMM_myo1,subsetted,file="ABCDS_Monocle2_12_23.RData")

