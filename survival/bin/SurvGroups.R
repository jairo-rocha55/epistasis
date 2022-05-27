args=commandArgs(trailingOnly=TRUE)
name1=args[1]
name2=args[2]

  
Mat=read.table("../data/Mat2",as.is=TRUE)
d=dim(Mat)
#Pairs=read.table("../data/Pairs2",as.is=TRUE)
surv=read.table("../data/luad_survival.tsv",as.is=TRUE,header=TRUE)
subjectlist = read.table("../data/files_list_405_TCGA.txt",as.is=TRUE)


  
  
groupssurvival <- function (name1,name2){
  
  pos1=which(Mat[,1]==name1)[1]
  pos2=which(Mat[,1]==name2)[1]
 
  G1=which((Mat[pos1,]==1 & Mat[pos2,]==1) | (Mat[pos1,]==1 & Mat[pos2,]==2) | (Mat[pos1,]==2 & Mat[pos2,]==1))-1
 
 G2=which(Mat[pos1,]==0 & Mat[pos2,]==0)
 # G2=setdiff( 1:(d[2]-1) , G1)


G1=merge(surv,subjectlist[G1,],by.x=6,by.y=1)[,c(1,4:5)]
G2=merge(surv,subjectlist[G2,],by.x=6,by.y=1)[,c(1,4:5)]

G1=cbind(G1,"YES")
G2=cbind(G2,"NO")
colnames(G2)<-colnames(G1)
Sol=rbind(G1,G2)

for (i in 1:(dim(Sol)[1])) {
   if (Sol[i,2]=="true")
       Sol[i,2]=0
   else Sol[i,2]=1
}

write.table(Sol[,c(1,4,2,3)],"../data/groupsSurv.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

}



groupssurvival(name1,name2)

