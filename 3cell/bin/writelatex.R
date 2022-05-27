library(xtable)
Pairs=read.table("../dataCOAD/GPairsPvalues")

XXX<- xtable(Pairs[1:20,c(1,2,19)], digits=-2)

print(XXX)
