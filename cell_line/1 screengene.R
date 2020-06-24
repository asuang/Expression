setwd("D:\\研究生\\密码子\\multi_cellline")
df = read.csv("rna_celline.csv", header = T)
df = df[,c(1,2,3,4)]
df1<-reshape(df,timevar="Sample",idvar="Gene",direction="wide")
s1<-function(DF1){
MAX<-data.frame()
for (n in c(1:56)) {
  #取不同细胞系基因表达量
  N<-DF1[,c((2*n),(2*n+1))]
  #对不同基因表达量排序，降序
  N<-N[order(-N[,2]),]
  #筛选出前1000个表达量高的基因
  N<-N[c(1:1000),]
  #从前1000个表达量高的基因中随机筛选出200个基因,返回的SN是向量
  SN<-sample(N[,1],400,replace = F)
  SND<-as.data.frame(SN)
  colnames(SND)<-c("Gene.name")
  MAX<-rbind(MAX,SND)
  }
return(MAX)
}


n=0
Max<-s1(df1)
#去除Max中重复的基因
Max<-unique(Max)
while (n<200) {
 
  M<-sample(Max[,1],100,replace = F)
  M<-as.data.frame(M)
  M<-t(M)
  write.table(M,"total.txt",quote=F,row.names = F,col.names = F,append = T)
  n=n+1
}

s2<-function(DF2){
  MIN<-data.frame()
  for (n in c(1:56)) {
    #取不同细胞系基因表达量
    N<-DF2[,c((2*n),(2*n+1))]
    #对不同基因表达量排序，升序
    N<-N[order(N[,2]),]
    #筛选出表达量大于1且小于等于5的基因
    N<-N[N[,2]>1&N[,2]<=10,]
    #从剩余的基因中随机筛选出200个基因,返回的SN是向量
    SN<-sample(N[,1],400,replace = F)
    SND<-as.data.frame(SN)
    colnames(SND)<-c("Gene.name")
    MIN<-rbind(MIN,SND)
  }
  return(MIN)
}

Min<-s2(df1)
#去除Max中重复的基因
Min<-unique(Min)
n=0
while (n<200) {
  
  #筛选出200个基因
  I<-sample(Min[,1],100,replace = F)
  I<-as.data.frame(I)
  I<-t(I)
  write.table(I,"total.txt",quote=F,row.names = F,col.names = F,append = T)
  n=n+1
}
