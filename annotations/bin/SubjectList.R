
S=read.table("COADnew/gdc_sample_sheet.tsv", as.is=TRUE,sep="\t",header=TRUE)

chopname <- function (x,rem) {
    len=nchar(x)
    return(substr(x,1,len-rem),".gz",sep="")
}

chop.gz <- function (x) {
  len=nchar(x)
  y=substr(x,len-2,len)
  if (y==".gz")
    return(substr(x,1,len-3))
  else 
    return(x)
}


M=sapply(S[,2],chop.gz)

Sol=data.frame(as.vector(M),S)[,c(1,7)]

write.table(file="COADnew/vcfList",Sol,quote = FALSE,col.names = FALSE,row.names = FALSE)

