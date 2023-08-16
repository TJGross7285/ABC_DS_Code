

####Phenotype WAVE 1: "Unable to determine"/ NIAD-10055, NIAD-10039, NIAD-10069, NIAD-10075
####No APOE Genotype/ NIAD-10109  

library(dplyr)

A<-read.csv("ABCDS_Pheno_6_2.csv",check.names=FALSE)
B<-read.csv("ADDS_10.1.20.csv",check.names=FALSE)
C<-read.csv("NIAD_Wave.2.csv",check.names=FALSE)
D<-read.csv("ABCDS_Clinical_Proteomics_2019.10.21.csv",check.names=FALSE)%>%filter(`NCRAD Sample ID`!="")

ADDS.1<-B%>%dplyr::filter(Cycle=="1")%>%dplyr::select(SubjectID)%>%unique()
ADDS.2<-B%>%dplyr::filter(Cycle=="2")%>%dplyr::select(SubjectID)%>%unique()
intersect(ADDS.1,ADDS.2)

NIAD.1.SID<-D%>%dplyr::pull(`Study ID`)%>%unique()
ADDS.IDs<-B%>%dplyr::pull(`SubjectID`)%>%unique()
intersect(ADDS.IDs,NIAD.1.SID)



NIAD.2<-C%>%dplyr::pull(`Sequence Number`)%>%unique()
NIAD.1.MINE<-A%>%dplyr::pull(`GeorgetownID`)%>%unique()
NIAD.1.SID<-D%>%dplyr::pull(`NCRAD Sample ID`)%>%unique()
intersect(NIAD.2,NIAD.1.MINE)
intersect(NIAD.2,NIAD.1.SID)

setdiff(intersect(NIAD.2,NIAD.1.SID),intersect(NIAD.2,NIAD.1.MINE))
setdiff(intersect(NIAD.2,NIAD.1.MINE),intersect(NIAD.2,NIAD.1.SID))



setdiff(NIAD.2,NIAD.1.SID)

ADDS.1.SID<-D%>%dplyr::pull(`Study ID`)%>%unique()
ADDS.1.ME<-A%>%dplyr::pull(`SampleID`)%>%unique()
setdiff(setdiff(ADDS.1.SID,ADDS.IDs),setdiff(ADDS.1.ME,ADDS.IDs))











library(dplyr)

B<-read.csv("ADDS_10.1.20.csv",check.names=FALSE)
C<-read.csv("NIAD_Wave.2.csv",check.names=FALSE)
D<-read.csv("ABCDS_Clinical_Proteomics_2019.10.21.csv",check.names=FALSE)%>%filter(`NCRAD Sample ID`!="")


ADDS.1.SID<-D%>%dplyr::pull(`Study ID`)%>%unique()
NIAD.1.SID<-D%>%dplyr::pull(`NCRAD Sample ID`)%>%unique()
ADDS.2<-B%>%dplyr::pull(`SubjectID`)%>%unique()
NIAD.2<-C%>%dplyr::pull(`Sequence Number`)%>%unique()

write.csv(intersect(ADDS.1.SID,ADDS.2),file="intersect_ADDS_2_3.csv")
write.csv(intersect(NIAD.1.SID,NIAD.2),file="intersect_NIAD_2_3.csv")
write.csv(setdiff(ADDS.2,ADDS.1.SID),file="new_ADDS_2_3.csv")
write.csv(setdiff(NIAD.2,NIAD.1.SID),file="new_NIAD_2_3.csv")


