setwd("D:\\�о���\\������\\multi_cellline")
df = read.csv("rna_celline.csv", header = T)
df = df[,c(1,2,3,4)]
df1<-reshape(df,timevar="Sample",idvar="Gene",direction="wide")
s1<-function(DF1){
MAX<-data.frame()
for (n in c(1:56)) {
  #ȡ��ͬϸ��ϵ���������
  N<-DF1[,c((2*n),(2*n+1))]
  #�Բ�ͬ������������򣬽���
  N<-N[order(-N[,2]),]
  #ɸѡ��ǰ1000���������ߵĻ���
  N<-N[c(1:1000),]
  #��ǰ1000���������ߵĻ��������ɸѡ��200������,���ص�SN������
  SN<-sample(N[,1],400,replace = F)
  SND<-as.data.frame(SN)
  colnames(SND)<-c("Gene.name")
  MAX<-rbind(MAX,SND)
  }
return(MAX)
}


n=0
Max<-s1(df1)
#ȥ��Max���ظ��Ļ���
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
    #ȡ��ͬϸ��ϵ���������
    N<-DF2[,c((2*n),(2*n+1))]
    #�Բ�ͬ�����������������
    N<-N[order(N[,2]),]
    #ɸѡ������������1��С�ڵ���5�Ļ���
    N<-N[N[,2]>1&N[,2]<=10,]
    #��ʣ��Ļ��������ɸѡ��200������,���ص�SN������
    SN<-sample(N[,1],400,replace = F)
    SND<-as.data.frame(SN)
    colnames(SND)<-c("Gene.name")
    MIN<-rbind(MIN,SND)
  }
  return(MIN)
}

Min<-s2(df1)
#ȥ��Max���ظ��Ļ���
Min<-unique(Min)
n=0
while (n<200) {
  
  #ɸѡ��200������
  I<-sample(Min[,1],100,replace = F)
  I<-as.data.frame(I)
  I<-t(I)
  write.table(I,"total.txt",quote=F,row.names = F,col.names = F,append = T)
  n=n+1
}