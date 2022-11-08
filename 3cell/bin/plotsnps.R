 

#setwd("/home/jairo/epistasis2/filtersnps/dataCOAD/")  #COAD or LOAD
setwd("/home/jairo/epistasis3/3cell/dataLUAD/")  #COAD or LOAD

#CCDC73  and HTR2B (left), and DDX4 and KCNJ16

# NUMSUBJ=405
# gene1="ACSM5"
# gene2="FLG2"
# gene3="MAGEB2"
# gene4="WASF3"
# chr1="Chr16"
# chr2="Chr1"
# chr3="ChrX"
# chr4="Chr13"

# NUMSUBJ=422
# gene1="CCDC73"
# gene2="HTR2B"
# gene3="DDX4"
# gene4="KCNJ16"
# chr1="Chr11"
# chr2="Chr2"
# chr3="Chr5"
# chr4="Chr17"

#(CHRM2, SLC6A15) and (PSD2, SUGP1)

NUMSUBJ=405
gene1="CHRM2"
gene2="SLC6A15"
gene3="PSD2"
gene4="SUGP1"
chr1="Chr7"
chr2="Chr12"
chr3="Chr5"
chr4="Chr17"

splitsnp <- function(x) 
{unlist(strsplit(x,c(":","\n")))}


plotgene <- function(gene1,chr1,NUMSUBJ,plotlegend=FALSE) {
  gene1file=paste("freq",gene1,"snps",sep="")
  SNPs=read.table(gene1file,as.is=TRUE)


A=sapply(SNPs[,1],splitsnp)

D=data.frame(A[1,],as.numeric(A[2,]),as.numeric(A[2,])+1,as.numeric(A[6,]),SNPs[,2],SNPs[,3])

names(D) <- c("chrom","start","end","score","freqnormal","freqtumor")
chromstart = min(D$start)
chromend = max(D$end)
print(D$chrom[1])
print(c(chromstart,chromend))
#maxy=max(D$freqnormal, D$score*NUMSUBJ)
#maxy<-max(maxy,D$freqtumor)*1.5

maxy=max(D$freqtumor)*1.05


offset=(chromend-chromstart)/200
offset=0
plot(D$start+offset,D$freqtumor,  col = "red", type="h", xlim=c(chromstart,chromend),ylim=c(1,maxy),#main=paste("gene", gene1), 
                ylab=paste("freq in", NUMSUBJ, "subj"),lwd=2, xaxt = "n",xlab="")#,xlab=paste("Chr",D$chrom[1],"pos"))
     #axes=FALSE,xlab="",ylab="",lwd=2)
#par(new=TRUE)
#plot(D$start,D$freqnormal,  col = "green", type="h", xlim=c(chromstart,chromend),ylim=c(1,maxy),log="",axes=FALSE,xlab="",ylab="",lwd=2)
#par(new=TRUE)
#plot(D$start-offset,D$score * NUMSUBJ, main=paste("gene", gene1), col = "blue", type="h",xlim=c(chromstart,chromend),ylim=c(1,maxy),log="",xlab=paste(chr1,"pos"),ylab=paste("freq in", NUMSUBJ, "subjects"),lwd=2)
axis(3)
title( xlab=paste("chr",D$chrom[1],"gene", gene1), line = 1.5)
if (plotlegend) #legend(x=chromend-(chromend-chromstart)*0.15,y=0.9*maxy,legend=c("MAF", "Normal tissue","Tumor tissue"),col=c("blue","green","red"), lty=c(1,1,1), cex=0.8)
  legend("topright",legend=c("MAF", "Normal tissue","Tumor tissue"),col=c("blue","green","red"), lty=c(1,1,1), cex=0.8)

}

# plotgene(gene4,chr1,NUMSUBJ)
# 
# 
# pdf(paste("MutLandscape",gene1,gene2,gene3,gene4,".pdf",sep=""))
# par(mfcol=c(4,1))
# plotgene(gene1,chr1,NUMSUBJ)
# plotgene(gene2,chr2,NUMSUBJ)
# plotgene(gene3,chr3,NUMSUBJ)
# plotgene(gene4,chr4,NUMSUBJ)
# 
# dev.off()

#setEPS()
#par(mar=c(0,0,0,0)+0.1)
pdf(paste("MutLandscape",gene1,".pdf",sep=""),height=2.5)
#par()
plotgene(gene1,chr1,NUMSUBJ)
dev.off()
pdf(paste("MutLandscape",gene2,".pdf",sep=""),height=2.5)
#par()
plotgene(gene2,chr1,NUMSUBJ)
dev.off()
pdf(paste("MutLandscape",gene3,".pdf",sep=""),height=2.5)
#par()
plotgene(gene3,chr1,NUMSUBJ)
dev.off()
pdf(paste("MutLandscape",gene4,".pdf",sep=""),height=2.5)
#par()
plotgene(gene4,chr1,NUMSUBJ)
dev.off()
