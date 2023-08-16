library(dplyr)
library(impute)
library(WGCNA)
library(sva)
library(limma)
library(vioplot)

A<-read.csv("MERGED_UPDATED_ABCDS_10_15_19_CHECK.csv",na.strings="")
colnames(A)[2]<-"Key"

B<-read.csv("Lipid_NIST.csv",check.names=FALSE)
B<-B%>%mutate(CV = gsub("%?","",B$`Coefficient of variation`))%>%filter(CV<=40)

E<-read.csv("Lipid_Data.csv",check.names=FALSE)
Eprime<-E%>%t()
colnames(Eprime)<-Eprime[1,]
SampleID<-colnames(E)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Eprime))
abunds<-bound[2:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

F<-read.csv("Lipid_PooledQC.csv",check.names=FALSE)
F<-F%>%mutate(CV = gsub("%?","",F$`Coefficient of variation`))%>%filter(CV<=40)
joint<-as.character(E$compound_name) %in% intersect(as.character(B$compound_name),as.character(F$compound_name))

abunds<-abunds[,joint==TRUE]
SampleIDprime<-colnames(E)[-c(1)]

Lipidizer_Use<-cbind(as.data.frame(SampleIDprime),as.data.frame(abunds))
colnames(Lipidizer_Use)[1]<-"Key"
Lipidizer_Use<-Lipidizer_Use%>%filter(Key!="Pooled QC" & Key!="NIST Plasma")

lipid_bound<-dplyr::inner_join(A,Lipidizer_Use,by="Key")
index<-duplicated(lipid_bound$`Key`)
lipid_bound<-lipid_bound[index==FALSE,]%>%filter(Consensus_DX!="Unable to determine")
write.csv(lipid_bound,file="ABCDS_Lipid_CV40_1024.csv")

##################################
##################################
##################################
B<-read.csv("Metab_NIST.csv",check.names=FALSE)
B<-B%>%mutate(CV = gsub("%?","",B$`Coefficient of variation`))%>%filter(CV<=40)

E<-read.csv("Metab_Data.csv",check.names=FALSE)
Eprime<-E%>%t()
colnames(Eprime)<-Eprime[1,]
SampleID<-colnames(E)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Eprime))
abunds<-bound[4:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

F<-read.csv("Metab_PooledQC.csv",check.names=FALSE)
F<-F%>%mutate(CV = gsub("%?","",F$`Coefficient of variation`))%>%filter(CV<=40)
joint<-as.character(E$Unique_name) %in% intersect(as.character(B$Unique_name),as.character(F$Unique_name))

abunds<-abunds[,joint==TRUE]
SampleIDprime<-colnames(E)[-c(1:3)]

Metab_Use<-cbind(as.data.frame(SampleIDprime),as.data.frame(abunds))
colnames(Metab_Use)[1]<-"Key"
Metab_Use<-Metab_Use%>%filter(Key!="Pooled QC" & Key!="NIST Plasma")

metab_bound<-dplyr::inner_join(A,Metab_Use,by="Key")
index<-duplicated(metab_bound$`Key`)
metab_bound<-metab_bound[index==FALSE,]%>%filter(Consensus_DX!="Unable to determine")
write.csv(metab_bound,file="ABCDS_Metab_CV40_1024.csv")

##################################
##################################
##################################

R<-read.csv("ABCDS_Metab_CV40_1024.csv",check.names=FALSE)
R$`Consensus_DX`<-as.factor(R$`Consensus_DX`)
levels(R$`Consensus_DX`)<-list("Control"=c("0-No dementia, no MCI"),"MCI.AD"=c("1-MCI","2-DEMENTIA"))
colnames(R)[9]<-"Sex"
edata<-t(log2(R[,-c(1:38)]))

pheno<-R[,c(5,7,9,15,12)]

####Impute missing data with KNN imputation
edata<-impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
mod<-model.matrix(~as.factor(Consensus_DX)
				  +as.factor(Sex)
				  +as.numeric(age),data=pheno)
mod0<-model.matrix(~as.factor(Sex)
				  +as.numeric(age),data=pheno)
n.sv<-num.sv(edata,mod,method="leek",seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(as.data.frame(pheno$Consensus_DX),as.data.frame(svobj))
colnames(total)<-c("Main")

design1<-as.data.frame(model.matrix(~0+Main,data=total))
colnames(design1)[1:2]<-c("Control","MCI.AD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm1 <-makeContrasts(
	`MCI.AD-Control`= MCI.AD-Control,
	levels=design1)



cm1 <-makeContrasts(`MCI-Control`= MCI-Control,
	levels=design1)
cm2<-makeContrasts(`AD-MCI`= AD-MCI,
	levels=design1)
cm3<-makeContrasts(`AD-Control` = AD-Control,
	levels=design1)

cmd<-makeContrasts(
	`MCI.AD-Control`= MCI.AD-Control,
	levels=design1)


fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)

fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)

T<-topTableF(fit1_F,adjust="BH",number=100000)
U<-topTableF(fit2_F,adjust="BH",number=100000)
Z<-topTableF(fit3_F,adjust="BH",number=100000)

plotting_abunds<-as.data.frame(t(edata))
vioplot(plotting_abunds$LACTATE_neg_1~pheno$Consensus_DX)

index<-colnames(plotting_abunds)%in% unique(c(T%>%filter(P.Value<.05)%>%rownames(),
				                              U%>%filter(P.Value<.05)%>%rownames(),	
                                              Z%>%filter(P.Value<.05)%>%rownames()))

pdf(file="plot_metab.pdf")
plotting_abunds_metab<-plotting_abunds[,index==TRUE]
for (i in 1:ncol(plotting_abunds_metab)){
	vioplot(plotting_abunds_metab[,i]~pheno$Consensus_DX, main=colnames(plotting_abunds_metab)[i])
}
dev.off()
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################

K<-read.csv("ABCDS_Lipid_CV40_1024.csv",check.names=FALSE)
K$`Consensus_DX`<-as.factor(K$`Consensus_DX`)
levels(K$`Consensus_DX`)<-list("Control"=c("0-No dementia, no MCI"),"MCI.AD"=c("1-MCI","2-DEMENTIA"))
colnames(K)[9]<-"Sex"
edata<-t(log2(K[,-c(1:38)]))
pheno<-K[,c(5,7,9,15,12)]

####Impute missing data with KNN imputation
edata<-impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
mod<-model.matrix(~as.factor(Consensus_DX)
				  +as.factor(Sex)
				  +as.numeric(age),data=pheno)
mod0<-model.matrix(~as.factor(Sex)
				  +as.numeric(age),data=pheno)
n.sv<-num.sv(edata,mod,method="leek",seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(as.data.frame(pheno$Consensus_DX),as.data.frame(svobj))
colnames(total)<-c("Main","SV1","SV2","SV3")
design1<-as.data.frame(model.matrix(~0+Main+SV1+SV2+SV3,data=total))
colnames(design1)[1:2]<-c("Control","MCI.AD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm1 <-makeContrasts(
	`MCI-Control`= MCI-Control,
	levels=design1)
cm2<-makeContrasts(`AD-MCI`= AD-MCI,
	levels=design1)
cm3<-makeContrasts(`AD-Control` = AD-Control,
	levels=design1)

cmd<-makeContrasts(
	`MCI.AD-Control`= MCI.AD-Control,
	levels=design1)


fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)

fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)

C<-topTableF(fit1_F,adjust="BH",number=100000)
P<-topTableF(fit2_F,adjust="BH",number=100000)
Q<-topTableF(fit3_F,adjust="BH",number=100000)

plotting_abunds<-as.data.frame(t(edata))
vioplot(plotting_abunds$`TAG56:6-FA18:2`~pheno$Consensus_DX)


plotting_abunds<-as.data.frame(t(edata))
vioplot(plotting_abunds$LACTATE_neg_1~pheno$Consensus_DX)

index<-colnames(plotting_abunds)%in% unique(c(C%>%filter(P.Value<.05)%>%rownames(),
				                              P%>%filter(P.Value<.05)%>%rownames(),	
                                              Q%>%filter(P.Value<.05)%>%rownames()))

pdf(file="plot_lipid.pdf")
plotting_abunds_lipid<-plotting_abunds[,index==TRUE]
for (i in 1:ncol(plotting_abunds_lipid)){
	vioplot(plotting_abunds_lipid[,i]~pheno$Consensus_DX, main=colnames(plotting_abunds_lipid)[i])
}
dev.off()

##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################
##############>>>>>>>>>>>>>>>>>>>>>>>>>###########################

R<-read.csv("ABCDS_Metab_CV40_1024.csv",check.names=FALSE)
R$`Consensus_DX`<-as.factor(R$`Consensus_DX`)
levels(R$`Consensus_DX`)<-list("Control"=c("0-No dementia, no MCI"),"MCI"=c("1-MCI"),"AD"=c("2-DEMENTIA"))
colnames(R)[9]<-"Sex"
edata_metab<-log2(R[,-c(1:38)])
MAD<-apply(edata_metab,2,mad)
median<-apply(edata_metab,2,median)
datExpr<-edata_metab###########(edata_metab-median)/sqrt(MAD)
pheno<-R[,c(5,7,9,15,12)]





K<-read.csv("ABCDS_Lipid_CV40_1024.csv",check.names=FALSE)
K$`Consensus_DX`<-as.factor(K$`Consensus_DX`)
levels(K$`Consensus_DX`)<-list("Control"=c("0-No dementia, no MCI"),"MCI"=c("1-MCI"),"AD"=c("2-DEMENTIA"))
colnames(K)[9]<-"Sex"
edata_lipid<-t(scale(log2(K[,-c(1:38)])))
pheno_lipid<-K[,c(5,7,9,15,12)]

####Should evaluate TRUE 
all.equal(pheno,pheno_lipid)

total_abunds<-rbind(edata_metab,edata_lipid)
####Impute missing data with KNN imputation
edata<-impute.knn(total_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
mod<-model.matrix(~as.factor(Consensus_DX)
				  +as.factor(Sex)
				  +as.numeric(age),data=pheno)
mod0<-model.matrix(~as.factor(Sex)
				  +as.numeric(age),data=pheno)
n.sv<-num.sv(edata,mod,method="leek",seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv

datExpr<-WGCNA::empiricalBayesLM(t(edata),removedCovariates=svobj,retainedCovariates=pheno$Consensus_DX)$adjustedData

####Run WGCNA functions to cluster lipids and plot resulting modules  
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 10
adjacency <- adjacency(datExpr, power = softPower,type="signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize <- 15
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                 deepSplit = 2, pamRespectsDendro = FALSE,
                 minClusterSize = minModuleSize);

####Name modules and replot; write out eigenlipid matrix 
table(dynamicMods)
set.seed(122)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05,
                     main = "Gene dendrogram and module colors")
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
All_eigenfeatures <- MEList$eigengenes
colorsUntargeted<-MEList$validColors


vioplot(All_eigenfeatures$`MEblack`~pheno$Consensus_DX)


for (i in 1:dim(All_eigenfeatures)[2]){
  bcorsis(all_eigenfeatures[,i]~as.numeric(pheno$pheno$Consensus_DX))
}



pdf(file="plot_wgcna.pdf")
for (i in 1:ncol(All_eigenfeatures)){
	vioplot(All_eigenfeatures[,i]~pheno$Consensus_DX, main=colnames(All_eigenfeatures)[i])
}
dev.off()
