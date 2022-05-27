library("R.utils")
PARAM <- commandArgs(trailingOnly=TRUE,asValues=TRUE)

MatrixGeneName=PARAM[[1]]
MatrixGeneNoSingle=PARAM[[2]]
qvalue=as.numeric(PARAM[[3]])
print(MatrixGeneName)
print(MatrixGeneNoSingle)
print(qvalue)

M=read.table(MatrixGeneName)

d=dim(M)
d1=d[1]
d2=d[2]
d1m=(d2-1)/2

illness=as.factor(c(rep(0,d1m),rep(1,d1m)))

calc<- function(kn) {
  A=t(M[kn,2:d2])
  n=sum(A==0)
  s=sum(A==1)
  
  return(s/(n+s))
}


PvalM=sapply(1:d1,calc);

#write.table(PvalM, "dataCOAD/proboftumormutation")
indO=order(PvalM,decreasing=TRUE)

#genpval=data.frame(1:d1,PvalM)

#genpval=genpval[indO,]


TopS=length(PvalM)*qvalue
#TopV=genpval[TopS,2]     # Under this value the gen by itself is related to tumor

MNosingle=M[-indO[1:TopS],]

write.table(MNosingle,MatrixGeneNoSingle,quote=FALSE,row.names=FALSE,col.names=FALSE)
