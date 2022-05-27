


S=read.table("BREAST/gdc_sample_sheet.2019-12-19.tsv", as.is=TRUE,sep=",")

L=read.table("BREAST/sizes",as.is=TRUE)

chopname <- function (x,rem) {
    len=nchar(x)
    return(substr(x,1,len-rem))
}

chop.gz <- function (x) {
  len=nchar(x)
  y=substr(x,len-2,len)
  if (y==".gz")
    return(substr(x,1,len-3))
  else 
    return(x)
}

rem=nchar(".exonic_variant_function")

N=sapply(L[,4],chopname,rem)

M=sapply(S[,2],chop.gz)


L1=data.frame(as.vector(N),L)

S1=data.frame(as.vector(M),S)
M=merge(L1,S1,by.x=1,by.y=1)

X=aggregate(M$V1.x,by=list(M$V6),max)

colnames(X)=c("names","size")
colnames(M)[c(2,11)]=c("size","names")
#Sol=merge(X,M,by.x=c(2,1),by.y=c(2,11))[,2:3]
library("plyr")
Sol=join(X,M, type = "left", match = "first")


write.table(file="BREAST/vcfList1",Sol[,c(1,3)],quote = FALSE,col.names = FALSE,row.names = FALSE)



