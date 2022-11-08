setwd("~/public/epistasis/3cell/dataLUAD")
X=read.table("DistributionSortedShortShort")
size=dim(X)[1]
plot(X[,1],1-(size:1)/size,col="red",log=c("y"))
points(X[,1],1-pchisq(X[,1],df=2.257))

# ks.test(X[,1],pchisq,df=2.257)

setwd("~/public/epistasis/3cell/dataCOAD")
X=read.table("GpvaluesSortedShort")
size=dim(X)[1]
plot(X[,1],1-(1:size)/size,col="red",log=c("y"))
points(X[,1],1-pchisq(X[,1],df=2.1),log=c("y"))

i= 9647869023  * 0.05/51098733   #9.8e-10 i=9  position 9 de la distribution which has the value 40.2 
# Se toman los mayores de 40.2 
# Los primeros 449 pasan. 

setwd("~/public/epistasis/3cell/dataLUAD")
X=read.table("Distpvalues")
size=dim(X)[1]
X[1,]=as.numeric(X[1,])
for( i in 2:size)  {
  X[i,1]<-X[i,1]+X[i-1,1]
}

for (i in 1:size) {
  X[i,1]<- X[i,1]/X[size,1]
}

write.table(X,"PValues",quote=FALSE,col.names=FALSE,row.names=FALSE)



setwd("~/public/epistasis/3cell/dataCOAD")
X=read.table("Distpvalues")
size=dim(X)[1]
X[1,]=as.numeric(X[1,])
for( i in (size-1):1)  {
  X[i,1]<-X[i,1]+X[i+1,1]
}

for (i in size:1) {
  X[i,1]<- X[i,1]/X[1,1]
}

Y=data.frame(X[size:1,1],X[size:1,2])
write.table(Y,"PValues",quote=FALSE,col.names=FALSE,row.names=FALSE)
