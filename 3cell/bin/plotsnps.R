 

#setwd("/home/jairo/epistasis2/filtersnps/dataCOAD/")  #COAD or LOAD
setwd("/home/jairo/epistasis2/filtersnps/dataLUAD/")  #COAD or LOAD

# NUMSUBJ=422
# gene1="TMEM178A"
# gene2="UGGT1"
# gene3="LYAN"
# chr1="Chr2"
# chr2="Chr2"
# chr3="Chr12"

NUMSUBJ=405
gene1="FRYL"
gene2="GEMIN2"
gene3="STK38"
chr1="Chr4"
chr2="Chr14"
chr3="Chr6"

splitsnp <- function(x) 
{unlist(strsplit(x,c(":","\n")))}


plotgene <- function(gene1,chr1,NUMSUBJ) {
  gene1file=paste("freq",gene1,"snps",sep="")
  SNPs=read.table(gene1file,as.is=TRUE)


A=sapply(SNPs[,1],splitsnp)

D=data.frame(A[1,],as.numeric(A[2,]),as.numeric(A[2,])+1,as.numeric(A[6,]),SNPs[,2],SNPs[,3])

names(D) <- c("chrom","start","end","score","freqnormal","freqtumor")
chromstart = min(D$start)
chromend = max(D$end)

maxy=max(D$freqnormal, D$score*NUMSUBJ)
maxy<-max(maxy,D$freqtumor)*1.1



offset=(chromend-chromstart)/200
plot(D$start+offset,D$freqtumor,  col = "red", type="h", xlim=c(chromstart,chromend),ylim=c(1,maxy),log="",axes=FALSE,xlab="",ylab="",lwd=2)
par(new=TRUE)
plot(D$start,D$freqnormal,  col = "green", type="h", xlim=c(chromstart,chromend),ylim=c(1,maxy),log="",axes=FALSE,xlab="",ylab="",lwd=2)
par(new=TRUE)
plot(D$start-offset,D$score * NUMSUBJ, main=paste("gene", gene1), col = "blue", type="h",xlim=c(chromstart,chromend),ylim=c(1,maxy),log="",xlab=paste(chr1,"pos"),ylab=paste("freq in", NUMSUBJ, "subjects"),lwd=2)

legend(x=chromend-(chromend-chromstart)*0.35,y=0.9*maxy,legend=c("MAF", "Normal tissue","Tumor tissue"),col=c("blue","green","red"), lty=c(1,1,1), cex=0.8)

}

pdf(paste("MutLandscape",gene1,gene2,gene3,".pdf",sep=""))
par(mfcol=c(3,1))
plotgene(gene1,chr1,NUMSUBJ)
plotgene(gene2,chr2,NUMSUBJ)
plotgene(gene3,chr3,NUMSUBJ)

dev.off()


