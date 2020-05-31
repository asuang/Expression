setwd("D:\\研究生\\tRNAseq\\GSM16248")
total=read.csv("merge21AA.csv",header=T)
n=1
total<-transform(total,normalized=rep(0,54))
N<-data.frame()
New<-data.frame()
while(n<54){
  for (i in seq(5)) {
    a=substr(total[n,1],1,3)
    b=substr(total[(n+i),1],1,3)
    if(a==b){
      next
    }else{
      N=total[n:(n+i-1),]
      N<-N[order(-N[,2]),]
      n=n+i
      break
    }
  
  }
  for (m in seq(nrow(N))) {
    N[m,3]=N[m,2]/N[1,2]
  }

  New<-rbind(New,N)
  N<-data.frame()
}  
N=total[52:54,]
N<-N[order(-N[,2]),]
for (m in seq(nrow(N))) {
  N[m,3]=N[m,2]/N[1,2]
}
New<-rbind(New,N)
write.csv(New,"normalize21.csv")

#汇总
total18=read.csv("merge18.csv",header=T)
total19=read.csv("merge19.csv",header=T)
total20=read.csv("merge20.csv",header=T)
total21=read.csv("merge21.csv",header=T)
S<-cbind(total18,total19,total20,total21)
colnames(S)<-c("control1_anti_codon","control1_num","control2_anti_codon","control2_num","treated1_anti_codon","treated1_num","treated2_anti_codon","treated2_num")
write.csv(S,"total_merge.csv")

total18=read.csv("normalize18.csv",header=T)
total18[,1]=NULL
total19=read.csv("normalize19.csv",header=T)
total19[,1]=NULL
total20=read.csv("normalize20.csv",header=T)
total20[,1]=NULL
total21=read.csv("normalize21.csv",header=T)
total21[,1]=NULL
S<-cbind(total18,total19,total20,total21)
colnames(S)<-c("control1_anti_codon","control1_num","control1_normalize","control2_anti_codon","control2_num","control2_normalize","treated1_anti_codon","treated1_num","treated1_normalize","treated2_anti_codon","treated2_num","treated2_normalize")
write.csv(S,"total_nomalize.csv")

