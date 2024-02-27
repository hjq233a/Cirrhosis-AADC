library(tidyverse)
library(dplyr)
#step1零值处理

origin<-read.csv("r.txt",header = T,row.names = 1)
col=ncol(origin);col
origin1<-origin[,3:ncol(origin)]+0.000000001
origin2<-cbind(origin[1],origin[2],origin1)
##########################################################
#step2计算效应值原始value计算，disease=n1, con=n2
filename <- data.frame(filename=c("NAFLD-Europe","PRJEB47902","PRJNA373901"))
result<-data.frame()
for(i in 1:nrow(filename))
{
  a1<-origin2[which(origin2$cohort==filename$filename[i] & origin2$group=="NAFLD"),]
  a2<-a1[,-c(1:2)]
  a3<- gather(a2, factor_key=TRUE)
  a4<-a3%>% group_by(key)%>%
    summarise(x1= mean(value), s1= sd(value), n1=length(value))
  b1<-origin2[which(origin2$cohort==filename$filename[i] & origin2$group=="Con"),]
  b2<-b1[,-c(1:2)]
  b3<- gather(b2, factor_key=TRUE)
  b4<-b3%>% group_by(key)%>%
    summarise(x2= mean(value), s2= sd(value), n2=length(value))
  c<-merge(a4,b4,by="key")
  result<-rbind(result,c)
}
f<-c("NAFLD-Europe","PRJEB47902","PRJNA373901")
x<-data.frame(cohort=f)
y=nrow(result)/nrow(x);y
x1<-result[1:y,]
x2<-result[(y+1):(2*y),-1]
x3<-result[(2*y+1):(3*y),-1]
x4<-result[(3*y+1):(4*y),-1]
x5<-result[(4*y+1):(5*y),-1]
result2<-cbind(x1,x2,x3,x4,x5)
write.csv(result2,"NAFLD_merge.csv",row.names=FALSE)

library(tidyr)
library(meta)
inputdata=read.csv("input3.txt",sep='\t',header=TRUE)
a<-as.data.frame(inputdata)
x<-a[,1:2]
result <- data.frame()
for (i in 3:ncol(a))
  #for (i in 3:10)
{
  y<-a[,i]
  c<-cbind(x,y)
  c<-as.data.frame(c)
  d<-spread(data=c,key=key,value=y)
  metawsd = metacont(n1,x1,s1,n2,x2,s2,data=d,sm="SMD",
                     comb.fixed = FALSE,comb.random = TRUE,
                     studlab = cohort)
  r <- data.frame(colnames(a)[i],metawsd$I2,
                  metawsd$pval.Q,metawsd$pval.random,
                  metawsd$TE.random,
                  metawsd$lower.random,
                  metawsd$upper.random,
                  metawsd$TE[1],metawsd$TE[2],
                  metawsd$TE[3] #队列数量#
                  
  )
  result <- rbind(result,r)
}
metawsd$studlab
names(result)<-c("Species","I^2","Q-Test","P-value","Coefficient",
                 "Confidence Interval Lower Limit","Confidence Interval Upper Limit",
                 "NAFLD-Europe","PRJEB47902","PRJNA373901"
)
write.csv(result, file = "nafld_coenficience_all.csv",row.names=FALSE)