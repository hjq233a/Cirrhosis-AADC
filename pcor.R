library(psych)
data<-read.csv("pcordata.txt",sep="\t")

library(ppcor)


x<-as.data.frame(data[data$HE!="0",])
gene_name1<-c()
gene_name2<-c()
cor_r<-c()
pvalue<-c()
qvalue<-c()
ciu<-c()
cil<-c()

for (i in 4:4)
{ for (j in 5:ncol(x))
{ 
  g1=colnames(x)[i]
  g2=colnames(x)[j]
  y1<-pcor.test(x[,i],x[,j],x[,3],method ="spearman")
  gene_name1=c(gene_name1,g1)
  gene_name2=c(gene_name2,g2)
  cor_r=c(cor_r,y1$estimate)
  pvalue=c(pvalue,y1$p.value)
  ciu=c(ciu,y1$conf.int[1])
  cil=c(cil,y1$conf.int[2])
}
}

qvalue=p.adjust(pvalue,'BH')
data_cor<-data.frame(gene_name1,gene_name2,cor_r,pvalue,qvalue)
write.csv(data_cor,"data.pcor.spearman.alb.csv")
