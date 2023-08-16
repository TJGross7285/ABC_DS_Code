E<-read.csv("abcds_wave_1_metabolomics_NEG_FINAL.csv",header=FALSE)
header_neg<-paste(paste(E[-c(1),1],E[-c(1),2],sep="_"),rep("NEG",dim(E)[1]-1),sep="/")
neg_abunds<-as.data.frame(t(E[-c(1),-c(1:2)]))
colnames(neg_abunds)<-header_neg
Cheema.Metabolomics.ID<-as.vector(unlist(E[1,-c(1:2)]))
total_NEG<-cbind(as.data.frame(Cheema.Metabolomics.ID),as.data.frame(neg_abunds))

H<-read.csv("abcds_wave_1_metabolomics_POS_FINAL.csv",header=FALSE)
header_pos<-paste(paste(H[-c(1),1],H[-c(1),2],sep="_"),rep("POS",dim(H)[1]-1),sep="/")
neg_abunds<-as.data.frame(t(H[-c(1),-c(1:2)]))
colnames(neg_abunds)<-header_pos
Cheema.Metabolomics.ID<-as.vector(unlist(H[1,-c(1:2)]))
total_POS<-cbind(as.data.frame(Cheema.Metabolomics.ID),as.data.frame(neg_abunds))

all.equal(total_NEG$Cheema.Metabolomics.ID,total_POS$Cheema.Metabolomics.ID)
total_combined<-cbind(total_NEG,total_POS[,-c(1)])

key<-read.csv("26046 plasma key_FINALII.csv")
joined<-dplyr::inner_join(key,total_combined,by="Cheema.Metabolomics.ID")[,-c(2:36)]
colnames(joined)[1]<-"ABCDS_ID"
write.csv(joined,file="Combined_ABCDS_WAVE1_Untargeted_Metabolomics.csv")