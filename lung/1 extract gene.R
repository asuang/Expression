setwd("D:\\�о���\\������\\rna_tissue_consensus.tsv")
#ʹ��read.table������ȡtsv�ļ�������
all<-read.table("rna_tissue_consensus.tsv",header = T,sep ="\t" )
all1<-reshape(all,timevar="Tissue",idvar="Gene",direction="wide")
lung<-all1[,c("Gene.name.lung","NX.lung")]
lung<-lung[order(-lung[,2]),]
#��ȡ�ٷֱȷ�λ��
quantile(lung$NX.lung,seq(0,1,0.05))
High<-as.data.frame(lung[c(1:1000),1])
Low<-lung[lung[,2]>=4.1&lung[,2]<=6.7,]
Low<-as.data.frame(sample(Low[,1],1000,replace = F))
for (i in seq(0,990,by=10)) {
  Tem<-t(as.data.frame(High[(i+1):(i+10),]))
  write.table(Tem,"lung_test.txt",append = T,quote = F,row.names = F,col.names = F)
}
for (i in seq(0,990,by=10)) {
  Tem<-t(as.data.frame(Low[(i+1):(i+10),]))
  write.table(Tem,"lung_test.txt",append = T,quote = F,row.names = F,col.names = F)
}