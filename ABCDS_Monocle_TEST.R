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



####Carry out dimension reduction/sample ordering using SVs as residual medel formula   
HSMM_myo<- reduceDimension(HSMM, reduction_method = 'DDRTree',norm_method="none", pseudo_expr=0,
                            residualModelFormulaStr="~as.numeric(age)+as.factor(Gender)", max_components=5)
HSMM_myo1 <- orderCells(HSMM_myo)

diffGenes<-differentialGeneTest(HSMM_myo1,
              fullModelFormulaStr = "~sm.ns(Pseudotime)")

sig_gene_names <- row.names(subset(diffGenes, qval < 0.1))

pdf(file="ABCDS_PseudotimeHeatmap_12_9.pdf")
Pseudotime_Heatmap_Data<-plot_pseudotime_heatmap(HSMM_myo1[sig_gene_names,],
                num_clusters = 10,
                cores = 100,
                show_rownames = FALSE)
dev.off()


save(Pseudotime_Heatmap_Data,file="Pseudotime_Heatmap_Data_12_23.RData")






















BEAM_res <- BEAM(HSMM_myo1, branch_point = 1, cores = 100)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf(file="ABCDS_Monocle_Heatmap.pdf")
plot_genes_branched_heatmap(HSMM_myo1[row.names(subset(BEAM_res,
                                          qval < 1e-4)),],
                                          branch_point = 1,
                                          num_clusters = 5,
                                          cores = 100,
                                          use_gene_short_name = TRUE,
                                          show_rownames = TRUE)
dev.off()



####Carry out dimension reduction/sample ordering w/o SVs 
HSMM_myo <- reduceDimension(HSMM, reduction_method = 'DDRTree',
							norm_method="none", pseudo_expr=0,max_components=5)
HSMM_myo2 <- orderCells(HSMM_myo)




####Plot Trajectories w/ and w/o SV adjustment 
pdf(file="Monocle_ABCDS_Untargeted_SVs.pdf")
plot_cell_trajectory(HSMM_myo1, color_by = "Consensus_DX")
plot_cell_trajectory(HSMM_myo1, color_by = "DiseaseState")
dev.off()

pdf(file="Monocle_ABCDS_Untargeted_noSVs.pdf")
plot_cell_trajectory(HSMM_myo2, color_by = "Consensus_DX")
plot_cell_trajectory(HSMM_myo2, color_by = "DiseaseState")
dev.off()


####CODE IN QUESTION 
####CODE RUNS BUT PLOTTING FUNCTION HANGS ON EXECUTION 
####Plot heatmap with respect to branchpoint 1; Cores?; bad plotting arguments? 
pdf(file="ABCDS_Monocle_Heatmap.pdf")
plot_genes_branched_heatmap(HSMM_myo1, 
							branch_point=1,
							num_clusters = 5,
                          	cores = 100,
                          	use_gene_short_name = T,
                          	show_rownames = T)
dev.off()
