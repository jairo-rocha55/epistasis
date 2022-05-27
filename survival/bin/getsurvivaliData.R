setwd("~/epistasis3/survival/")
X=read.csv("SurvivalTCGA-CDR-SupplementalTableS1.csv",sep="\t")
B=X[X$type=="BRCA",]
C=X[X$type=="COAD",]
L=X[X$type=="LUAD",]
Corig=read.csv("dataCOAD/coad_survival.tsv",sep="\t")
Z=merge(C,Corig,by.x="bcr_patient_barcode",by.y="submitter_id")
Z[,c("OS.time","OS","time","censored")]
B[,c("OS.time","OS")]

# 2c30dc20-18a8-44f9-ab10-558c1e5634dc	TCGA-COAD	true	6	1	TCGA-AY-6196	
# 26a12266-c6d8-42f9-bd4e-5093d827ac9a	TCGA-COAD	true	14	1	TCGA-AM-5820	

n=dim(B)[1]
filename=rep("filename",n)
censored=sapply(B$OS,function(x) if (x==1) "false" else "true")
score=rep(0.0,n)
write.table(data.frame(filename,B$type,censored,B$OS.time,score,B$bcr_patient_barcode),
            "dataBREAST/breast_survival.tsv",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
