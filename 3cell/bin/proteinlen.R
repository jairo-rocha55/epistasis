setwd("/home/jairo/epistasis3/3cell/dataCOAD/")
L=read.csv("protlen.tab",sep="\t",as.is=TRUE)
setlen <- function(a,len){
  return ((c(a,len)))
}
explodelen <- function(x) {
  y=strsplit(as.character(x[2])," ")
  return(lapply(y[[1]],setlen,x[3]))
}

S=as.data.frame(t(as.data.frame(unlist(apply(L,1,explodelen),recursive = FALSE))))

S[,2]=as.numeric(S[,2])


M=read.table("MatrixGene")
Len=merge(M[,1],S,by=1)
T=read.table("tmp1")
Len2=merge(T[,1],S,by=1)
mean(Len2[,2])
mean(Len[,2])
