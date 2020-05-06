fre<-read.csv("lungfre_test.csv",header =T )
fre<-transform(fre,expression=rep(c(1,0),each=100))
#使用glm()进行回归分析
FF=glm(formula = expression~ATA+ATT,data=fre,binomial(link='logit'),control=list(maxit=100))
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
FF=glm(formula = expression~ATA+ATT,data=train,binomial(link='logit'),control=list(maxit=100))
summary(FF)
real<-test$expression
pre<-predict(FF,type="response",newdata = test)
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

