## Setting directory & aktifkan package ##
setwd("E:/seqadit/newDEG")

## Gabungngkan data gene abundance ##
read.table("R.txt",sep="\t",header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)->Resistan
read.table("S.txt",sep="\t",header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)->Sakit
merge(Resistan,Sakit, by="V1")->RS
write.table(RS,"ResSus.txt",sep="\t",row.names=FALSE, col.names=FALSE,quote = FALSE)

## jalankan edgeR ##
library(edgeR)
read.table("ResSus.txt", header = F, row.names = 1) -> deg
head(deg)
c("Resistan", "Sakit")-> names(deg)
factor(c("Resistan", "Sakit"))->group
DGEList(counts=deg, genes=rownames(deg), group=group)->y
calcNormFactors(y)->y1
y1$samples
0.2->bcv
exactTest(y1, dispersion=bcv^2)->et
et$table->DEGfix
summary(de<-decideTestsDGE(et))
topTags(et, 21819)$table->DE2 ## angka sesuaikan dengan jumlah Down + NotSig + Up
write.table(DE2, file="deg_rarb.txt", sep="\t", na="NA", quote=FALSE)

## merge anotasi dan deg ##
read.delim("deg_resistan_sakit.txt", header=T)->d ## file deg
read.delim("swissprot.txt", header=T)->ref ## file anotasi
merge(d, ref, by="genes")->gabung
head(gabung)
write.table(gabung, file="deg_anot_ResSak.txt", sep="\t", na="NA", quote=FALSE)
