setwd("~/epistasis3/survival/dataCOAD")

args=commandArgs(trailingOnly=TRUE)
args[1]="../../3cell/dataCOAD/GMat"
args[2]="SurvSolSorted"
args[3]="filelistTCGA422"
args[4]="../SurvivalTCGA-CDR-SupplementalTableS1.csv"

args[5]="SubjectSubList"

SNPxSampleMat=args[1]
SNPpairs=args[2]

FileList=args[3]
DataSub=args[4]
OutFile=args[5]

Mat=read.table(SNPxSampleMat,as.is=TRUE)
Pairs=read.table(SNPpairs,as.is=TRUE)[1:3,] # [1:100,]

subjectlist = read.table(FileList,as.is=TRUE)
subData=read.csv(DataSub,sep="\t")

  
d=dim(Mat)
  
#library(reticulate)
  
groupssubjects <- function (j){
  name1=Pairs[j,1]
  name2=Pairs[j,2]  
  pos1=which(Mat[,1]==name1)[1]
  pos2=which(Mat[,1]==name2)[1]
 

  G1=which((Mat[pos1,]==1 & Mat[pos2,]==1) | (Mat[pos1,]==1 & Mat[pos2,]==2) | (Mat[pos1,]==2 & Mat[pos2,]==1))-1  # ss, sb, bs
# G1=which((Mat[pos1,]>=1 & Mat[pos2,]>=1) )-1  # both mutated in tumor 

  return(subjectlist[G1,1])
}
j=2
Sol=merge(groupssubjects(j),subData,by.x=1,by.y="bcr_patient_barcode")
#Sol=sapply(1:dim(Pairs)[1], groupssubjects)
# FinalPairs=data.frame(Pairs[,1],Pairs[,2],Pairs[,dim(Pairs)[2]],Sol)
# Ind=order(FinalPairs$Sol)
# FinalPairs<- FinalPairs[Ind,]
# write.table(FinalPairs,OutFile,quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(Sol,paste("SubjectSupport",Pairs[j,1],Pairs[j,2],".csv",sep=""),sep="\t",quote = FALSE,row.names=FALSE)


#install.packages("tables")
# library(xtable)
# 
# XXX<- xtable(FinalPairs,file="table.tex", digits=-2)
# print(XXX)
