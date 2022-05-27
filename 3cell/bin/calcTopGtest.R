args=commandArgs(trailingOnly=TRUE)
fileinput=args[1]
fileout=args[2]
qvalue=as.numeric(args[3])
experiments=as.numeric(args[4])

X=read.table(fileinput)
#MAX=pmax(X[,12],X[,13],X[,14])
#PV1=1-pnorm(X[,12])
#PV2=1-pnorm(X[,13])
#PV3=1-pnorm(X[,14])
#MIN=pmin(PV1,PV2,PV3)

PV = 1-pchisq(X[,12],df=3) # df=1 si no hay germline ceros en last column and row 
# one degree of freedom is fixed by minimization, should it be 3?4 is less conservative?

Y=cbind(X,PV) 

posPV=dim(Y)[2]

Bonferroni <- function(V,Q) {
  maxexperiments=length(V)
  cte=Q/(experiments+maxexperiments)
  for (i in 1:maxexperiments)
    if (V[i]  > cte) 
      return(i-1)
  return(maxexperiments)
}



Benjamini <- function(V,Q) {
  maxexperiments=length(V)
  cte=(experiments+maxexperiments)/Q
  for (i in 1:maxexperiments)
    if ((V[i]*cte)  > i) 
      return(i-1)
  return(maxexperiments)
}

indO=order(Y[,12],decreasing=TRUE)
Y=Y[indO,]
#TopS=Benjamini(Y[,posPV],qvalue)
TopS=Bonferroni(Y[,posPV],qvalue)
TopV=Y[TopS,posPV] 

print(TopS)
print(TopV)


write.table(Y,fileout,quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(Y[1:TopS,],paste(fileout,"Top",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
