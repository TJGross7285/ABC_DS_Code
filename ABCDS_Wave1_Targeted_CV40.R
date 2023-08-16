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
cmd<-makeContrasts(
	`MCI.AD-Control`= MCI.AD-Control,
	levels=design1)

fit1_F <- contrasts.fit(fit1, cmd)
fit1_F <- eBayes(fit1_F,trend=TRUE)

T<-topTableF(fit1_F,adjust="BH",number=100000)


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

##################################
##################################
##################################
##################################

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
cmd<-makeContrasts(
	`MCI.AD-Control`= MCI.AD-Control,
	levels=design1)


fit1_F <- contrasts.fit(fit1, cmd)
fit1_F <- eBayes(fit1_F,trend=TRUE)

C<-topTableF(fit1_F,adjust="BH",number=100000)

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









##################################
##################################
##################################
##################################

R<-read.csv("ABCDS_Metab_CV40_1024.csv",check.names=FALSE)
R$`Consensus_DX`<-as.factor(R$`Consensus_DX`)
levels(R$`Consensus_DX`)<-list("Control.MCI"=c("0-No dementia, no MCI","1-MCI"),"AD"=c("2-DEMENTIA"))
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
colnames(design1)[1:2]<-c("Control.MCI","AD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cmd<-makeContrasts(
	`AD-Control.MCI`= AD-Control.MCI,
	levels=design1)


fit1_F <- contrasts.fit(fit1, cmd)
fit1_F <- eBayes(fit1_F,trend=TRUE)

X<-topTableF(fit1_F,adjust="BH",number=100000)

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

##################################
##################################
##################################
##################################

K<-read.csv("ABCDS_Lipid_CV40_1024.csv",check.names=FALSE)
K$`Consensus_DX`<-as.factor(K$`Consensus_DX`)
levels(K$`Consensus_DX`)<-list("Control.MCI"=c("0-No dementia, no MCI","1-MCI"),"AD"=c("2-DEMENTIA"))
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
colnames(total)<-c("Main","SV1")
design1<-as.data.frame(model.matrix(~0+Main+SV1,data=total))
colnames(design1)[1:2]<-c("Control.MCI","AD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)

cmd<-makeContrasts(
	`AD-Control.MCI`= AD-Control.MCI,
	levels=design1)


fit1_F <- contrasts.fit(fit1, cmd)
fit1_F <- eBayes(fit1_F,trend=TRUE)

Z<-topTableF(fit1_F,adjust="BH",number=100000)

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


##################################
##################################
##################################
##################################


R<-read.csv("ABCDS_Metab_CV40_1024.csv",check.names=FALSE)
R$`Consensus_DX`<-as.factor(R$`Consensus_DX`)
levels(R$`Consensus_DX`)<-list("Control"=c("0-No dementia, no MCI"),"MCI"=c("1-MCI"),"AD"=c("2-DEMENTIA"))
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
colnames(design1)[1:3]<-c("Control","AD","MCI")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cmd<-makeContrasts(
	`AD-Control`= AD-Control,
	levels=design1)


fit1_F <- contrasts.fit(fit1, cmd)
fit1_F <- eBayes(fit1_F,trend=TRUE)

G<-topTableF(fit1_F,adjust="BH",number=100000)

plotting_abunds<-as.data.frame(t(edata))

