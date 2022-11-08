args=commandArgs(trailingOnly=TRUE)
SNPxSampleMat=args[1]
SNPpairs=args[2]
Surv=args[3]
FileList=args[4]
OutFile=args[5]

Mat=read.table(SNPxSampleMat,as.is=TRUE)
Pairs=read.table(SNPpairs,as.is=TRUE) # [1:100,]
surv=read.table(Surv,as.is=TRUE,header=TRUE)
subjectlist = read.table(FileList,as.is=TRUE)

  
d=dim(Mat)
  
library(reticulate)
  
groupssurvival <- function (j){
  name1=Pairs[j,1]
  name2=Pairs[j,2]  
  pos1=which(Mat[,1]==name1)[1]
  pos2=which(Mat[,1]==name2)[1]
 

 G2=which((Mat[pos1,]>=1 & Mat[pos2,]==0) | (Mat[pos1,]==0 & Mat[pos2,]>=1) ) # One mutated in tumor an other other not.  
#| (Mat[pos1,]==2 & Mat[pos2,]==1))-1  # ss, sb, bs
 G1=which((Mat[pos1,]>=1 & Mat[pos2,]>=1) )-1  # both mutated in tumor 
# G2=which((Mat[pos1,]<=0 & Mat[pos2,]<=0) )-1  # any not mutated in tumor 
 
# G1=setdiff( 1:(d[2]-1) , G2)
 # G2=setdiff( 1:(d[2]-1) , G1)

G1=merge(surv,subjectlist[G1,],by.x=6,by.y=1)[,c(1,4:5)]
G2=merge(surv,subjectlist[G2,],by.x=6,by.y=1)[,c(1,4:5)]

if (dim(G1)[1]==0) 
  return(1.0)
if (dim(G2)[1]==0) 
  return(1.0)

   G1=cbind(G1,"YES")
   G2=cbind(G2,"NO")
   colnames(G2)<-colnames(G1)
   Sol=rbind(G1,G2)
for (i in 1:(dim(Sol)[1])) {
   if (Sol[i,2]=="true")
       Sol[i,2]=1
   else Sol[i,2]=0
}
write.table(Sol[,c(1,4,2,3)],"TmpgroupsSurv.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

 py_run_file("survival.py")
#print(py$PVALUE)
 return (py$PVALUE)
}

#RASA4B JAG2 189 5 74 7 2 11 44 24 66 0 1.84 5.67 69.77 7.13987602374289e-09 5.67

Sol=sapply(1:dim(Pairs)[1], groupssurvival)
FinalPairs=data.frame(Pairs[,1],Pairs[,2],Pairs[,dim(Pairs)[2]],Sol)
Ind=order(FinalPairs$Sol)
FinalPairs<- FinalPairs[Ind,]
write.table(FinalPairs,OutFile,quote=FALSE,col.names=FALSE,row.names=FALSE)
#install.packages("tables")
library(xtable)

XXX<- xtable(FinalPairs,file="table.tex", digits=-2)
print(XXX)
