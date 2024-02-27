data <- read.csv("cordata.txt", sep = "\t",check.names=F, header = TRUE)
x<-as.data.frame(data)
#以下步骤是大循环的关键//--约两天6460*6460#

gene_name1<-c()
gene_name2<-c()
cor_r<-c()
pvalue<-c()
qvalue<-c()
ciu<-c()
cil<-c()


for (i in 3:ncol(x))
{ 
  g1=colnames(x)[i]
  g2=colnames(x)[2]
  y1=cor.test(as.numeric(x[,i]),as.numeric(x[,2]),method ="spearman")
  gene_name1=c(gene_name1,g1)
  gene_name2=c(gene_name2,g2)
  cor_r=c(cor_r,y1$estimate)
  pvalue=c(pvalue,y1$p.value)
  ciu=c(ciu,y1$conf.int[1])
  cil=c(cil,y1$conf.int[2])
}


qvalue=p.adjust(pvalue,'BH')
data_cor<-data.frame(gene_name1,gene_name2,cor_r,pvalue,qvalue)
write.csv(data_cor,file="spearman.cor.csv",sep='\t',row.names=FALSE)

