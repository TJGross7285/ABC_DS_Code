############################################################(7/17) ADDS Pheno Setup 
library(dplyr)
library(sva)
library(limma)
library(impute)


####Read in NIAD to Georgetown mappings 
K<-read.csv("26046 plasma key_FINAL.csv")[,-c(1)]
colnames(K)[1:2]<-c("SubjectID","GeorgetownID")

####Read in APOE genotype data and dichotomize on E4 allele status 
join<-K%>%filter(GeorgetownID!="NA")
G<-read.csv("niad_apoe_master_6.3.2019.csv",na.strings="")
colnames(G)[1]<-"SubjectID"

HasE4<-rep(0,dim(G)[1])
index<-G$"A1"=="E4"|G$"A2"=="E4"
HasE4[index==TRUE]<-1
G<-cbind(G,HasE4)

join_NIAD<-dplyr::full_join(join,G,by="SubjectID")

J<-read.csv("Final_Pheno.csv")
colnames(J)[1]<-"SubjectID"

####Write final file to spreadsheet; Unclear what is meant by #NULL! code in CONSENSUS_DX field/not in redcap codebook/filter out for now 
final_join_NIAD<-dplyr::inner_join(join_NIAD,J,by="SubjectID")
final_join_NIAD$consensus<-as.factor(final_join_NIAD$consensus)
index<-duplicated(final_join_NIAD$SubjectID)
final_join_NIAD<-final_join_NIAD[index==FALSE,]



NIAD_FOR_MERGE<-final_join_NIAD[,c(1,2,3,14)]
colnames(NIAD_FOR_MERGE)[4]<-"DiseaseState"

T<-read.csv("ABCDS_UPDATED.csv")
T$STUDY<-gsub("NiAD","NIAD",T$STUDY)
NIAD_table<-T%>%filter(STUDY=="NIAD")
colnames(NIAD_table)[1]<-"SubjectID"


NIAD_FOR_MERGE<-dplyr::inner_join(NIAD_FOR_MERGE,NIAD_table)[,c(1,2,13,6)]
colnames(NIAD_FOR_MERGE)[3:4]<-c("HasE4","DiseaseState")

####################################################
####################################################
substrRight <- function(x, n){
  substr(x, nchar(as.character(x))-n+1, nchar(as.character(x)))
}

####Read in NIAD to Georgetown mappings 
K<-read.csv("26046 plasma key_FINAL.csv")[,-c(1)]
colnames(K)[1:2]<-c("SubjectID","GeorgetownID")
K<-cbind(as.data.frame(substrRight(K$SubjectID,4)),K)%>%filter(grepl("^ADDS.+",GeorgetownID))
colnames(K)[1]<-"study_id"
K$study_id<-as.character(K$study_id)



adds_pheno<-read.csv("ADDS_Pheno_8_12.csv",check.names=FALSE)
adds_apoe<-read.csv("ADDS_APOE_8_12.csv",check.names=FALSE,na.strings="")
T<-dplyr::inner_join(adds_apoe,adds_pheno)%>%select("study_id","CONSENSUS_DX_CYCLE1","apoe_genotype")
T$study_id<-as.character(T$study_id)

final_join_ADDS<-dplyr::inner_join(K,T,by="study_id")
index<-duplicated(final_join_ADDS$SubjectID)
final_join_ADDS<-final_join_ADDS[index==FALSE,]
ApoE4<-grepl("4",as.character(final_join_ADDS$apoe_genotype))*1
ApoE4[is.na(final_join_ADDS$apoe_genotype)==TRUE]<-NA

final_join_ADDS<-cbind(final_join_ADDS,ApoE4)
T<-read.csv("ABCDS_UPDATED.csv")
T$STUDY<-gsub("NiAD","NIAD",T$STUDY)
ADDS_table<-T%>%filter(STUDY=="ADDS")
colnames(ADDS_table)[1]<-"study_id"

merge_adds<-dplyr::inner_join(final_join_ADDS,ADDS_table)

ADDS_FOR_MERGE<-final_join_ADDS[,c(2,3,6,4)]
colnames(ADDS_FOR_MERGE)[3:4]<-c("HasE4","DiseaseState")

final_pheno<-rbind(ADDS_FOR_MERGE,NIAD_FOR_MERGE)

write.csv(final_pheno,file="Pheno_8_14.csv")



Hi Anne,

My name is Tom Gross and work with Dr. Mark 
































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
pheno<-combined[,c(1:31)]
pheno$DiseaseState<-as.factor(factor(pheno$DiseaseState))
index<-is.na(pheno$ANYE4)==FALSE
pheno<-pheno[index,]
edata<-edata[,index]

####Impute missing data with KNN imputation
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data
mod<-model.matrix(~as.factor(DiseaseState)+as.factor(factor(ANYE4)),data=pheno)
mod0<-model.matrix(~as.factor(factor(ANYE4)),data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(as.data.frame(pheno$DiseaseState),as.data.frame(svobj))
colnames(total)<-c("Main","SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10","SV11","SV12","SV13","SV14")

####Set up accessory objects for DE incorporating pre-post timepoints
design1<-as.data.frame(model.matrix(~0+Main+SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14,data=total))
colnames(design1)[1:2]<-c("Control","MCIAD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm <- makeContrasts(
	`MCIAD-Control` = MCIAD-Control,
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"

write.csv(T, file="DS_DE_ABCDS_10_16_19.csv")

colnames(T)[2]<-"Log2FC"

library(EnhancedVolcano)
pdf(file="WAVE1_ABCDS_Lipidyzer_Volcano_10_16_19.pdf")
EnhancedVolcano(T,lab=T$Feature,x="Log2FC",y="P.Value",xlim=c(-.5,.5),ylim=c(0,10),FCcutoff=.2)
dev.off()

library(WGCNA)
datExpr<-t(edata)
clean<-empiricalBayesLM(data=datExpr,removedCovariates=total[,-c(1)],retainedCovariates=total[,1])$adjustedData

### Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

###Call the network topology analysis function
sft <- pickSoftThreshold(clean, powerVector = powers, verbose = 5)

###Soft power is 3
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

###Define soft thresholding power to be used
softPower <- 3
adjacency <- adjacency(clean, power = softPower)

###Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency);
dissTOM <- 1-TOM

### Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

###Define module size
minModuleSize <- 25

### Identify modules using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

### Convert numeric labels into colors
set.seed(122)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

#### Calculate module eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

###Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



### Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs);

### Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average");

### Plot the result
pdf(file="Plot1.pdf")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()


### Define module dissimilarity threshold 
MEDissThres <- 0.6
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
moduleColors <- mergedColors

### Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50));
moduleLabels <- match(moduleColors, colorOrder)-1;
MEs <- mergedMEs


index<-as.factor(merge$colors)
levels(index)
colnames(clean[,index=="red"])

cor<-cor(clean,MEs,method="spearman")
hubLipids<-chooseTopHubInEachModule(clean, index)
















library(caret)
library(pROC)
library(WGCNA)
expr<-empiricalBayesLM(t(edata),removedCovariates=svobj,retainedCovariates=pheno$DiseaseState)$adjustedData
#### Set up Caret control object and outcome vector for classification

cctrl1 <- trainControl(method = "repeatedcv", repeats=10, number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

####Train random forest; 10-fold CV repeated 10 times
set.seed(122)
RandomForest<-train(modeling_abunds,DiseaseState, method="rf", trControl = cctrl1, metric = "ROC")
























Y<-read.csv("ptdemog.csv",na.strings="")
colnames(Y)[1]<-"SubjectID"

final_join<-dplyr::full_join(join,Y)
K<-read.csv("26046 plasma key_FINAL.csv")
colnames(K)[2:3]<-c("SubjectID","GeorgetownID")
fin_join<-dplyr::full_join(K,final_join)%>%arrange(GeorgetownID)
write.csv(fin_join,file="ABCDS_Pheno_7_11.csv")







write.csv(fin_join,file="ABCDS_Pheno_7_11.csv")















T<-read.csv("Pheno_Sheet_PreLim_7_11.csv")

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

combined<-inner_join(T,bind,by="GeorgetownID")
na<-apply(combined,2,is.na)
na_index<-apply(na,2,sum)/dim(combined)[1]
combined<-combined[,na_index<.333]

edata<-log2(t(apply(combined[,-c(1:10)],2,as.numeric)))
pheno<-combined[,c(1:10)]
pheno$HasDS<-as.factor(pheno$HasDS)
levels(pheno$HasDS)<-list("Dementia"=c("1"),"Control"=c("0"))

####Impute missing data with KNN imputation
edata<-impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
mod<-model.matrix(~as.factor(HasDS),data=pheno)
mod0<-model.matrix(~1,data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(as.data.frame(pheno$HasDS),as.data.frame(svobj))
colnames(total)<-c("Main","SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8")

####Set up accessory objects for DE incorporating pre-post timepoints
design1<-model.matrix(~0+Main+SV1+SV2+SV3+SV4+SV5,data=total)
colnames(design1)[1:2]<-c("Control","MCIAD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm <- makeContrasts(
	`MCIAD-Control` = MCIAD-Control,
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"

write.csv(T, file="DS_DE_7_17.csv")
