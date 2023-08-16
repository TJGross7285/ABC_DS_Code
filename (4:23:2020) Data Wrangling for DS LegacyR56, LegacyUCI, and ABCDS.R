############################(4/23/2020) Data Wrangling for DS LegacyR56, LegacyUCI, and ABCDS  
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
####DO NOT OVERWRITE "MATCHING_TEMPLATE.csv"!!!!!!!!!!!!!!!!!!!



###################
###################TO REPLICATE DOWNSTREAM WORK, START FROM HERE ON 
library(dplyr)
library(WGCNA)
library(impute)


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

Study<-as.factor(c(rep("R56",dim(join_r56)[1]),rep("ABCDS",dim(join_abcds)[1]),rep("Legacy",dim(join_legacy)[1])))

####NO APOE DATA FOR LEGACYUCI; MCI NOT CODED FOR LEGACYUCI; TALK TO/GET FROM ERIC DORAN? 
H<-rbind(join_r56[,c(1:2,16,12,14)],join_abcds[,c(1:2,7,11,9)],join_legacy[,c(1:2,9,6,5)])
H<-cbind(Study,H)
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
####SHOULD WE FILTER FEATURES BY STUDY ======> MERGE DATA ? UNCLEAR (but currently implemented)
X<-read.csv("Key.csv",check.names=FALSE)[,-c(1)]

T<-read.csv("threeGroup_NEG_abunds.csv",check.names=FALSE)
colnames(T)[1]<-"ID1"
join_a<-dplyr::inner_join(X,T,by="ID1")

neg_index_R56<-apply(apply(join_a[join_a$Study=="R56",],2,is.na),2,sum)/dim(join_a[join_a$Study=="R56",])[1]
neg_index_ABCDS<-apply(apply(join_a[join_a$Study=="ABCDS",],2,is.na),2,sum)/dim(join_a[join_a$Study=="ABCDS",])[1]
neg_index_Legacy<-apply(apply(join_a[join_a$Study=="Legacy",],2,is.na),2,sum)/dim(join_a[join_a$Study=="Legacy",])[1]

join_a<-join_a[,neg_index_R56<.333 & neg_index_ABCDS<.333 & neg_index_Legacy<.333]

U<-read.csv("threeGroup_POS_abunds.csv",check.names=FALSE)
colnames(U)[1]<-"ID2"
join_b<-dplyr::inner_join(X,U,by="ID2")

pos_index_R56<-apply(apply(join_b[join_b$Study=="R56",],2,is.na),2,sum)/dim(join_b[join_b$Study=="R56",])[1]
pos_index_ABCDS<-apply(apply(join_b[join_b$Study=="ABCDS",],2,is.na),2,sum)/dim(join_b[join_b$Study=="ABCDS",])[1]
pos_index_Legacy<-apply(apply(join_b[join_b$Study=="Legacy",],2,is.na),2,sum)/dim(join_b[join_b$Study=="Legacy",])[1]

join_b<-join_b[,pos_index_R56<.333 & pos_index_ABCDS<.333 & pos_index_Legacy<.333]


total_join<-cbind(join_a,join_b[,-c(1:6)])
total_join<-total_join[total_join$DiseaseState!="Unable to determine",]

write.csv(total_join,file="ThreeGroup_4_29.csv")


####POSSIBLE DATA TRANSFORMATION STEPS   
####Impute missing data with KNN imputation
E<-read.csv("ThreeGroup_5_6.csv",check.names=FALSE)
edata<-log2((t(E[,-c(1:6)])))
edata<-impute.knn(edata,k = 10,rng.seed=362436069)$data
datExpr<-t(as.data.frame(edata))

####adjust for study effects; residualize metabolite matrix using empirical Bayes linear model  
datExpr<-empiricalBayesLM(datExpr,removedCovariates=E$Study)$adjustedData