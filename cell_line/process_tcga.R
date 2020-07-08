setwd("/home/song/Data/download/lungRNAseq_TCGA/download_dir")
list<-read.table("txtlist.txt")
d<-read.table(list[1,],row.names = 1)
for (i in seq(2,nrow(list))) {
  a<-read.table(list[i,],row.names = 1)
  d<-cbind(d,a)
}
colnames(d)<-c(1:30)
write.csv(d,"merge.csv")
d<-read.csv("merge.csv",header = T)
tem<-row.names(d)
for ( i in seq(1,length(tem))){
  tem[i]<-sub("\\..*","",tem[i])
}
row.names(d)<-tem
write.csv(d,"merge.csv")
library(org.Hs.eg.db)
library(AnnotationHub)
hs<-org.Hs.eg.db
keytypes(org.Hs.eg.db)#see the keytype you can convert in org.Hs.eg.db
a<-AnnotationDbi::select(hs, keys = tem,columns= c("SYMBOL","ENSEMBL"),keytype = "ENSEMBL")
library(clusterProfiler)
b<-AnnotationDbi::mapIds(hs, keys = tem,column= c("SYMBOL","ENSEMBL"),keytype = "ENSEMBL")
symbol <- bitr(tem, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
BiocManager::install("ensembldb")
library(ensembldb)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
symbol<-ensembldb::select(EnsDb.Hsapiens.v86,keys=tem,keytype = "GENEID",columns= c("SYMBOL","GENEID"))
data<-merge(symbol,d,by.x="GENEID",by.y="X")
a<-cbind(a,d[,1:30])
a<-unique(a)
write.csv(data,"lung.csv")
##STEP1
data<-read.csv("lung.csv",header = T,row.names = 1)
s1<-function(DF1){
  MAX<-data.frame()
  for (n in c(3:32)) {
    #取不同细胞系基因表达量
    N<-DF1[,c(2,n)]
    #对不同基因表达量排序，降序
    N<-N[order(-N[,2]),]
    #筛选出前1000个表达量高的基因
    N<-N[c(1:600),1]
    ND<-as.data.frame(N)
    colnames(ND)<-c("Gene.name")
    MAX<-rbind(MAX,ND)
  }
  return(MAX)
}


n=0
Max<-s1(data)
#去除Max中重复的基因
Max<-unique(Max)
while (n<100) {
  
  M<-sample(Max[,1],100,replace = F)
  M<-as.data.frame(M)
  M<-t(M)
  write.table(M,"total.txt",quote=F,row.names = F,col.names = F,append = T)
  n=n+1
}


s2<-function(DF2){
  MIN<-data.frame()
  for (n in c(3:32)) {
    #取不同细胞系基因表达量
    N<-DF2[,c(2,n)]
    #对不同基因表达量排序，升序
    N<-N[order(N[,2]),]
    #筛选出表达量大于1且小于等于5的基因
    N<-N[N[,2]>5&N[,2]<=10,]
    #从剩余的基因中随机筛选出200个基因,返回的SN是向量
    SN<-sample(N[,1],400,replace = F)
    SND<-as.data.frame(SN)
    colnames(SND)<-c("Gene.name")
    MIN<-rbind(MIN,SND)
  }
  return(MIN)
}

Min<-s2(data)
#去除Max中重复的基因
Min<-unique(Min)
n=0
while (n<100) {
  
  #筛选出200个基因
  I<-sample(Min[,1],100,replace = F)
  I<-as.data.frame(I)
  I<-t(I)
  write.table(I,"total.txt",quote=F,row.names = F,col.names = F,append = T)
  n=n+1
}

##STEP3
codon<-read.csv("count.csv",header=T)
#remove ATG TGA TAA TAG TGG
codon<-transform(codon,ATG=NULL,TGA=NULL,TAA=NULL,TAG=NULL,TGG=NULL)

fre<-data.frame()
sum=0

for (i in seq(200)) {
  for (m in seq(59)) {
    sum=sum+codon[i,m]
  }
  for (n in seq(59)){
    a<-codon[i,n]/sum
    fre[i,n]<-a
  }
  sum=0
  
}
colnames(fre)<-c('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT','ATA', 'ATC', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT','TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT')
write.csv(fre,"fre.csv",row.names = F)


##STEP4
fre<-read.csv("fre.csv",header =T )
fre<-transform(fre,expression=rep(c(1,0),each=100))
#使用glm()进行回归分析
FF=glm(formula = expression~.,data=fre,binomial(link='logit'),control=list(maxit=100))
print(summary(FF))

#画列线图
#使用rms包lrm()进行回归分析
library("rms")
#设定nomogram的参数,对数据进行打包
Y<-rep(c(1,0),each=40)
ddist<-datadist(fre)
options(datadist='ddist')
#logistics 回归
f<-lrm(formula = expression ~ ATA+ATT,fre)
#绘图
nom<-nomogram(f,lp.at=seq(-3.5,2,by=0.5),fun=plogis,fun.at = c(0.05, seq(0.0,1.0, by=0.2),0.95),funlabel = "Expression")
plot(nom,col.grid = gray(c(0.8,0.95)))
jpeg(filename = "lungAAGAAA_test.jpeg")

#对模型进行预测，训练组70%，预测组30%（随机抽样）
set.seed(1234)
ind<-sample(x=2,size=nrow(fre),replace = T,prob = c(0.7,0.3))
#分别生成训练组和预测组
train<-fre[ind==1,]
test<-fre[ind==2,]
FF=glm(formula = expression~.,data=train,binomial(link='logit'),control=list(maxit=100))
summary(FF)
real<-test$expression
pre<-predict(FF_back,type="response",newdata = test)
#用ROCR包
library("ROCR")
#把预测结果和真实结果形成一个prediction对象
pred<-prediction(pre,real)
auc<-performance(pred,"auc")@y.values
perf<-performance(pred,"tpr","fpr")
X11(width=10,height=10)
plot(perf,col='red',xlim=c(0,1),ylim=c(0,1),main="ATA&ATT ROC",lwd=2)
abline(a=0,b=1,lty=2,col='black')

#用pROC包
library("pROC")
mroc<-roc(real,pre,smooth=T)
X11(width=10,height=10)
plot(mroc,print.auc=T,auc.polygon=F,legacy.axes=T,grid=c(0.1,0.2),
     grid.col=c("green","red"),max.auc.polygon=T,
     print.thres=F,main="Lung",font.label=2)

#将lung的模型汇合在一张图上
plot(mroc,print.auc=T,auc.polygon=T,legacy.axes=T,grid=c(0.1,0.2),
     grid.col=c("green","red"),max.auc.polygon=T,auc.polygon.col="linen",
     print.thres=F,main="Lung",font.label=2,col='yellowgreen')

par(font.lab=3,cex.axis=0.7,cex.lab=1,cex.main=2,font.main=1)
FF=glm(formula = expression~ATA+ATC+ATT,data=train,binomial(link='logit'),control=list(maxit=100))
summary(FF)
pre<-predict(FF,type="response",newdata = test)
pred<-prediction(pre,real)
performance(pred,"auc")@y.values
mroc<-roc(real,pre,smooth=T)
plot(mroc,print.auc=F,auc.polygon=F,legacy.axes=T,grid=c(0.1,0.2),main="Lung",
     font.label=2,col='yellowgreen')
#add1
FF=glm(formula = expression~CAC+CAT,data=train,binomial(link='logit'),control=list(maxit=100))
summary(FF)
pre<-predict(FF,type="response",newdata = test)
pred<-prediction(pre,real)
performance(pred,"auc")@y.values
mroc<-roc(real,pre,smooth=T)
plot(mroc,print.auc=F,auc.polygon=F,legacy.axes=T,grid=c(0.1,0.2),add=TRUE,
     font.label=2,col='hotpink')
#add2
FF=glm(formula = expression~GGA+GGC+GGG+GGT,data=train,binomial(link='logit'),control=list(maxit=100))
summary(FF)
pre<-predict(FF,type="response",newdata = test)
pred<-prediction(pre,real)
performance(pred,"auc")@y.values
mroc<-roc(real,pre,smooth=T)
plot(mroc,print.auc=F,auc.polygon=F,legacy.axes=T,grid=c(0.1,0.2),add=TRUE,
     font.label=2,col='darkslategray3')
legend("bottomright",legend = c("I  0.771","H 0.833","G 0.760"),
       col=c("yellowgreen","darkslategray3","hotpink"),lwd = 2,bty="n")


#通过方差膨胀系数判断共线性问题
library("car")
vif(FF_back)
kappa(fre[,1:59])

#主成分回归解决共线性
prin<-princomp(~ACC + AGT + CGA + GAA + GCT,data=fre,cor=T,scores=T)
summary(prin,loadings = T)
fre<-cbind(fre,prin$scores[,1:4])
fre<-transform(fre,expression=rep(c(1,0),each=200))
FP<-glm(expression~Comp.1+Comp.2+Comp.3,data=fre,binomial(link="logit"))
summary(FP)     

#逐步回归分析
library("MASS")
FF=glm(formula = expression~ACC + AGT + CGA + GAA ,data=fre,binomial(link='logit'),control=list(maxit=100))
FF_back<-stepAIC(FF,direction = "backward",steps=3000)
FF=glm(formula = expression~1,data=fre,binomial(link='logit'),control=list(maxit=100))
FF_for<-stepAIC(FF,direction="forward",scope=list(upper=~AAA+AAC+AAG+AAT+ACA+ACC+ACG+ACT+
                                                    AGA+AGC+AGG+AGT+ATA+ATC+ATT+CAA+CAC+CAG+
                                                    CAT+CCA+CCC+CCG+CCT+CGA+CGC+CGG+CGT+CTA+
                                                    CTC+CTG+CTT+GAA+GAC+GAG+GAT+GCA+GCC+GCG+
                                                    GCT+GGA+GGC+GGG+GGT+GTA+GTC+GTG+GTT+TAC+
                                                    TAT + TCA + TCC +TCG+TCT+TGC+TGT+TTA+TTC+TTG+
                                                    TTT,lower=~1),steps = 2000)
summary(FF_for)
#或者通过step（）函数
FF_for<-step(FF,scope=list(upper=~AAA+AAC+AAG+AAT+ACA+ACC+ACG+ACT+AGA+AGC+AGG+AGT+ATA+ATC+ATT+CAA+CAC+CAG+CAT+CCA+CCC+CCG+CCT+CGA+CGC+CGG+CGT+CTA+CTC+CTG+CTT+GAA+GAC+GAG+GAT+GCA+GCC+GCG+GCT+GGA+GGC+GGG+GGT+GTA+GTC+GTG+GTT+TAC+TAT + TCA + TCC +TCG+TCT+TGC+TGT+TTA+TTC+TTG+TTT,lower=~1),direction="forward")
FF_back<-step(FF,direction = "both",test="F")

#绘制logistic S型曲线
train<-transform(train,z=4964+14966*AAG-255312*AGT-150133*CAC)
train<-transform(train,z=-8347+208115*ACC-533196*AGT+306466*CGA+103238*GAA+304393*GCT)
fre<-transform(fre,z=3892-112660*CAC-42249*GAA-141723*GTA)
fre<-transform(fre,y=1/(1+exp(-z)))
plot(y~z,data=fre)

#皮尔逊系数,查看各个变量之间的关系
x<-fre[,c('ACC','AGT','CGA','GAA','GCT')]
cor(x,x)