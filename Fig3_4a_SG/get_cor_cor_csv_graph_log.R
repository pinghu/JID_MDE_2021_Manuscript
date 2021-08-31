#options(echo=TRUE) # if you want see commands in output file
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
print(args)


filename=args[1]
outname=args[2]
filename ="ShortGenus6Data.SG.Clean"
outname ="RA6SG"
rm(args)

A<-read.table(filename, sep="\t", header=TRUE)

d <- dim(A);
B <- data.matrix(A[1:d[1],2:d[2]]);
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
test_data=log10(B+ZZ)
xx<-t(test_data)

#B<-cor(xx)
colnames(xx) =A[,1]
rownames(xx)=colnames(A)[2:d[2]]

library("Hmisc")
rCOR=rcorr(as.matrix(xx), type="spearman")
write.table(rCOR$P,file=paste0(outname,'.corP.csv'),sep=",",row.names=A[,1], col.names = A[,1])
write.table(rCOR$r,file=paste0(outname,'.corR.csv'),sep=",",row.names=A[,1], col.names = A[,1])
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

ff<-flattenCorrMatrix(rCOR$r, rCOR$P)
write.table(ff,file=paste0(outname,'.corFlat.tsv'),sep="\t")

library(corrplot)
res <- cor(xx)
round(res, 2)
#corrplot(res, type = "upper", order = "hclust", 
#         tl.col = "black", tl.srt = 45)


tiff(filename = paste0("Fig4a.tiff"), width = 2660, height = 2660, res=300)
corrplot.mixed(rCOR$r, p.mat=rCOR$P, sig.level=0.05, diag='u', insig = "blank",  order="hclust", upper = "ellipse", lower.col = "black", tl.pos="lt")
dev.off()
