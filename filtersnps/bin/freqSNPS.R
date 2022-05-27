args <- commandArgs(trailingOnly=TRUE)
name<- args[1]
M=read.table(name)
 d=dim(M)
 d2=(d[2]-1)/2
 normal=apply(M[,2:d2+1],1,sum)
 tumor=apply(M[,(d2+2):d[2]],1,sum)
 freq=data.frame(M[,1],normal,tumor)
 write.table(freq,paste("freq",name,sep=""))
 
