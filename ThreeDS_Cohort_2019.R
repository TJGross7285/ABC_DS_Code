#########################################(11/15/2019) Combined DE of R56, Legacy, and ABCDS 
library(dplyr)
library(WGCNA)
library(impute)


#####Read-in data to produce merging key from A. Cheema GWU XCMS output files; write to intermed. file 
#####for manual edits  
Y<-read.csv("2.neg_3_batch_xcms_IS.csv",header=FALSE)
Z<-read.csv("2.pos_3_batch_xcms_is.csv",header=FALSE)
G<-cbind(as.data.frame(t(Y[1,])), as.data.frame(t(Z[1,])))
write.csv(G,file="COMPLETE_PHENO.csv")


####DO NOT OVERWRITE "COMPLETE_PHENO.csv"!!!!!!!!!!!!!!!!!!!
G<-read.csv("COMPLETE_PHENO.csv")
index<-grep("Plasma",as.character(G$ID1))
G<-as.data.frame(G[-index,])
G<-G[G$Study!="ABCDS_Sib",]
write.csv(G,file="MATCHING_TEMPLATE.csv")



###################
###################



####Generate Matching Key between Clinical and Metabolomics Data 
matching_template<-read.csv("MATCHING_TEMPLATE.csv",check.names=FALSE)[,-c(1)]
legacy<-matching_template[matching_template$Study=="LegacyUCI",]
abcds<-matching_template[matching_template$Study=="ABCDS",]
r56<-matching_template[matching_template$Study=="LegacyR56",]
colnames(legacy)[1:2]<-c("UCI ID","Georgetown ID")
colnames(r56)[1:2]<-c("ParID","IncrementedParID")
colnames(abcds)[1:2]<-c("GeorgetownID","ModGeorgetownID")

abcds_pheno<-read.csv("MERGED_UPDATED_ABCDS_10_15_19_CHECK.csv",check.names=FALSE)
legacy_pheno<-read.csv("pheno_table.csv",check.names=FALSE)
r56_pheno<-read.csv("R56_DS_Clinical_Variables_Metabolomics_2019.08.27.csv",check.names=FALSE)
legacy_pheno$`Georgetown ID`<-as.factor(legacy_pheno$`Georgetown ID`)
r56_pheno$`ParID`<-as.factor(r56_pheno$`ParID`)

join_r56<-r56%>%inner_join(r56_pheno,by="ParID")%>%filter(AD_status!="")
join_r56<-join_r56[duplicated(join_r56$ParID)==FALSE,]
join_abcds<-abcds%>%inner_join(abcds_pheno,by="GeorgetownID")
join_legacy<-legacy%>%inner_join(legacy_pheno,by="Georgetown ID")

colnames(join_r56)[c(1:2,16,12,14)]<-c("ID1","ID2","DiseaseState","Sex","Age")
colnames(join_abcds)[c(1:2,7,11,9)]<-c("ID1","ID2","DiseaseState","Sex","Age")
colnames(join_legacy)[c(1:2,9,6,5)]<-c("ID1","ID2","DiseaseState","Sex","Age")

H<-rbind(join_r56[,c(1:2,16,12,14)],join_abcds[,c(1:2,7,11,9)],join_legacy[,c(1:2,9,6,5)])
levels(H$Sex)<-list("Male"=c("Male","M"),"Female"=c("Female","F"))
write.csv(H,file="Key.csv")



###################
###################



####Ingest A. Cheema XCMS metabolomics data and write to useable analysis files; CALL AS "NA" MEASUREMENTS WITH ZERO ABUNDANCE  
neg_abunds<-read.csv("2.neg_3_batch_xcms_IS.csv",check.names=FALSE,na.strings="0")%>%t()
colnames(neg_abunds)<-paste(read.csv("2.neg_3_batch_xcms_IS.csv",header=FALSE)[-c(1),1],"NEG",sep="/")
write.csv(neg_abunds, file="threeGroup_NEG_abunds.csv")

pos_abunds<-read.csv("2.pos_3_batch_xcms_is.csv",check.names=FALSE,na.strings="0")%>%t()
colnames(pos_abunds)<-paste(read.csv("2.pos_3_batch_xcms_is.csv",header=FALSE)[-c(1),1],"POS",sep="/")
write.csv(pos_abunds, file="threeGroup_POS_abunds.csv")



###################
###################



####Merge expression and clinical data; Write to CSV; DROP FEATURES WITH >33.3% NA VALUES OVERALL
####SHOULD WE FILTER FEATURES BY STUDY ======> MERGE ? UNCLEAR 
X<-read.csv("Key.csv",check.names=FALSE)[,-c(1)]

T<-read.csv("threeGroup_NEG_abunds.csv",check.names=FALSE)
colnames(T)[1]<-"ID1"
join_a<-dplyr::inner_join(X,T,by="ID1")
neg_index<-apply(apply(join_a,2,is.na),2,sum)/dim(join_a)[1]
join_a<-join_a[,neg_index<.333]

U<-read.csv("threeGroup_POS_abunds.csv",check.names=FALSE)
colnames(U)[1]<-"ID2"
join_b<-dplyr::inner_join(X,U,by="ID2")
pos_index<-apply(apply(join_b,2,is.na),2,sum)/dim(join_b)[1]
join_b<-join_b[,pos_index<.333]

total_join<-cbind(join_a,join_b[,-c(1:4)])
total_join<-total_join[total_join$DiseaseState!="Unable to determine",]

write.csv(total_join,file="ThreeGroup_4_23.csv")


####POSSIBLE DATA TRANSFORMATION STEPS   
####Impute missing data with KNN imputation
edata<-log2(apply(t(total_join[,-c(1:4)]),2,as.numeric))
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data
datExpr<-t(as.data.frame(edata))

####adjust for study effects; residualize metabolite matrix using empirical Bayes linear model  
datExpr<-empiricalBayesLM(datExpr,removedCovariates=total_join$Study)$adjustedData






pdf(file="plots.pdf")
for (i in 1:length(J)){
	vioplot(datExpr[,J[i]]~pheno$DiseaseState,main=colnames(total_join)[J[i]])
}
dev.off()

pdf(file="plots_MEs.pdf")
for (i in 1:dim(MEs)[2]){
	vioplot(MEs[,i]~pheno$DiseaseState,main=colnames(MEs)[i])
}
dev.off()






Mode<-c(rep("ESI-",dim(join_a)[2]-4),rep("ESI+",dim(join_b)[2]-4))
MZ_RT<-c(colnames(T)[-c(1)],colnames(U)[-c(1)])

######*********NEED TO ADJUST BY features removed in thresholding 
Y<-read.csv("meta_three.csv")[,-c(1)]
RT<-gsub("^.+_","",Y$MZ_RT)
MZ<-gsub("_.+$","",Y$MZ_RT)
Y<-cbind(as.data.frame(seq(1,dim(Y)[1])),as.data.frame(MZ),as.data.frame(RT),as.data.frame(Mode))
colnames(Y)[1]<-"Feature"
Y$Feature<-as.factor(Y$Feature)


####Impute missing data with KNN imputation; n.sv=6
edata<-log2(apply(t(total_join[,-c(1:4)]),2,as.numeric))
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data









pheno<-total_join[,1:4]
levels(pheno$DiseaseState)<-list("DS-AD"=c("1-MCI","2-DEMENTIA","AD","MCI"),
								"DS-NAD"=c("0-No dementia, no MCI","No AD","Non-AD"))

mod<-model.matrix(~as.factor(DiseaseState),data=pheno)
mod0<-model.matrix(~1,data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=6)$sv

design_table<-cbind(as.data.frame(pheno$DiseaseState),as.data.frame(svobj))
colnames(design_table)<-c("pheno","SV1","SV2","SV3","SV4","SV5","SV6")
design1<-as.data.frame(model.matrix(~0+pheno+SV1+SV2+SV3+SV4+SV5+SV6,data=design_table))
colnames(design1)[1:2]<-c("DSAD","DSNAD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm <- makeContrasts(
 	`DSAD-DSNAD` = DSAD-DSNAD,
 	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"





T<-dplyr::inner_join(Y,T)
colnames(T)[c(2,3,5)]<-c("m/z","retention time","statistic")
write.csv(T, file="DS_DE_threeGroup_11_15.csv")


write.table(T%>%filter(Mode=="ESI+")%>%select("m/z","retention time","P.Value","statistic"),
			sep="\t",file="DS_ThreeGroup_POS_DSAD.DSNAD.txt",row.names=FALSE)
write.table(T%>%filter(Mode=="ESI-")%>%select("m/z","retention time","P.Value","statistic"),
			sep="\t",file="DS_ThreeGroup_NEG_DSAD.DSNAD.txt",row.names=FALSE)



cd "/Users/TGross/Desktop/ThreeGroup" 
source activate py27
mummichog -f DS_ThreeGroup_POS_DSAD.DSNAD.txt -o DS_ThreeGroup_POS_DSAD.DSNAD -m positive -u 7 
mummichog -f DS_ThreeGroup_NEG_DSAD.DSNAD.txt -o DS_ThreeGroup_NEG_DSAD.DSNAD -m negative -u 7 



