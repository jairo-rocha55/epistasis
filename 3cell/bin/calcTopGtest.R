args=commandArgs(trailingOnly=TRUE)
fileinput=args[1]
fileout=args[2]
qvalue=as.numeric(args[3])
experiments=as.numeric(args[4])
Pvaluesfile=args[5]

X=read.table(fileinput)
#MAX=pmax(X[,12],X[,13],X[,14])
#PV1=1-pnorm(X[,12])
#PV2=1-pnorm(X[,13])
#PV3=1-pnorm(X[,14])
#MIN=pmin(PV1,PV2,PV3)

Pvalues=read.table(Pvaluesfile)  
sizePvalues<-dim(Pvalues)[1]

pvalue <- function(chisq) {
i=1
while (i < sizePvalues & Pvalues[i,2] > chisq) {
  i<- i+1 
} 
  if (Pvalues[i,2]==chisq)
    return(Pvalues[i,1])
  else if (i==1)
    return(0)
  else
    return(Pvalues[i-1,1])

}

PV=sapply(X[,12],pvalue)

#PV = 1-pchisq(X[,12],df=3) # df=1 si no hay germline ceros en last column and row 
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
