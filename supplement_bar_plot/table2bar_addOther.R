############################################
#need to try split panel
#https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html
#########################################
rm(list=ls())
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
#filename="RA2.ave.top.data"
#outname="RA2MDE"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
newsum=A[,2:d[2]]
bugs=A[,1]
row.names(newsum)=A[,1]
ttt=colnames(A)[1]
ddd=dim(newsum)
library(RColorBrewer)
col1=brewer.pal(9, "Set1")
col2=brewer.pal(12, "Set3")
col3=brewer.pal(9, "Pastel1")
col4=brewer.pal(8, "Accent")
col5=brewer.pal(8, "Dark2")
cols=c(col1, col2, col3, col4, col5)[1:(ddd[1]+1)]

Other=1-colSums(newsum)
newsum2<-rbind(newsum, Other)

tiff(filename=paste0(outname,".bar.tiff"), width=3600, height=1800,res=300 ) 

#png(filename=paste0(outname,".bar.png"), width=1800, height=1800,res=300 ) 
par(mar=c(11.8,4.1,4.1,2.1))

barplot(as.matrix(newsum2), col=cols, width=2, ylab="relative abundance", las=2, cex.names=1)

dev.off()

tiff(filename=paste0(outname,".lengend.tiff"), width=1800, height=1800,res=300)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("topleft",bugs, cex=0.8, fill=cols, title=ttt)
dev.off()

