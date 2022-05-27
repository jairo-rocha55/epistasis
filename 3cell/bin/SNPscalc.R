
#M= matrix(c(1,-2,0,4,-13,1,2,-8,1),nrow=3)

#N=matrix(c(1,-4,0,1,-4,1,0,1,0),nrow=3)

X=read.table("SNPsSubjects",as.is=TRUE)
Y=read.table("freqSnps",header=TRUE,as.is=TRUE)

Z=strsplit(Y[,2],":")

W=lapply(Z,function(x) x[2])
A=as.data.frame(W)
A=as.numeric(unlist(W))

B=merge(X,A,by.x=6,by.y=1)

write.table(B,"SNPSusedSubjects",col.names=FALSE,row.names=FALSE,quote=FALSE)
