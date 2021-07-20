##########################

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#filename="L2RA"
library(ggplot2)
library("ggpubr")
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
truefc<-function(VVV){
	XXX=VVV
	if(VVV==0){
	    XXX=NA
   	}else if(VVV<1){
	    XXX=-1/VVV
    	}
	return(XXX)
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=A[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=log10(B+ZZ)
#C=B

Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Site=rep("NA", Clen)
AgeGrp=rep("NA", Clen)
Age=rep("NA", Clen)
SID=rep("NA", Clen)
Order=rep("NA", Clen)
SiteAgeGrp=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
       Site[mm]=splitname[[mm]][1]
       AgeGrp[mm]=splitname[[mm]][3]
       Age[mm]=splitname[[mm]][4]
       SID[mm]=splitname[[mm]][5]
       Order[mm]=splitname[[mm]][6]
       SiteAgeGrp[mm]=paste0(Site[mm], AgeGrp[mm])
}


listMean<-c( "Mean_All", "shortTaxa", "CorP_Face_SG", "Cor_Face_SG")
m = matrix(, nrow = d[1], ncol = length(listMean))
rownames(m)=A[,1]
colnames(m)=listMean

Age=as.numeric(Age)
SG=as.numeric(B[1,])
mydata=data.frame(AgeGrp, Age,SG, Site,  SiteAgeGrp)


tiff(filename="SG.age.boxplot.tiff", width=1200, height=1600,res=300)  
p1<-ggplot(mydata, aes(x=AgeGrp, y=SG, color=AgeGrp), alpha=0.02) +geom_boxplot(color="black")+geom_jitter(position=position_jitter(0.2))+theme_bw()+geom_smooth(data=mydata,aes(x = AgeGrp, y = SG, group=1), color="orange", alpha=0.4, level=0.95, size = 0.8, method = "loess", span = 1) +ggtitle("SG Area vs. Age Group")+ labs(x = "", y="Sebaceous gland (SG) area measurements (square um)" )+ theme(legend.position = "none") 
print(p1)
dev.off()
	
tiff(filename="SG.age.cor.tiff", width=1200, height=1600, res=300)
p2<-ggplot(mydata, aes(x=Age, y=log10(SG))) +theme_bw()+ggtitle("log10SG ~ Age Correlation")+geom_point(color=AgeGrp) +geom_smooth(method='loess', color="orange")
print(p2)
dev.off()


for (i in 1:d[1]){
        genename=A[i,1]
	gene0=as.numeric(B[i,])
        gene<-as.numeric(C[i,])
        splitG<-strsplit(as.character(genename), "[;]")
        LLL=length(splitG[[1]])
        genus=splitG[[1]][LLL]
	
        mydata=data.frame(gene,gene0, AgeGrp, Age,SG, Site, SID, Order, SiteAgeGrp)


       tiff(filename=paste0(genus,".log10SG.cor.tiff"), width=1200, height=1600, res=300)
       p2<-ggplot(mydata, aes(x=SG, y=gene)) +theme_bw()+ggtitle(paste0(genus, ":RA ~ SG Correlation"))+geom_point(color=AgeGrp) +geom_smooth(method='loess')
       print(p2)
       dev.off()
       	          
       count=+1
       m[i, count]=mean(gene0)
       count=count+1
       m[i, count]=genus
      
       count=count+1
       m[i,count]<-cor.test(gene, SG, method="spearman")$p.value
       count=count+1
       m[i,count]<-cor(gene, log10(SG), method="spearman")
      
}
write.table(m, file=paste0(filename, ".stat"),sep="\t")
