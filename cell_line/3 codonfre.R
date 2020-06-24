codon<-read.csv("count.csv",header=T)
#È¥³ýATG£¬TGA£¬TAA£¬TAG£¬TGG
codon<-transform(codon,ATG=NULL,TGA=NULL,TAA=NULL,TAG=NULL,TGG=NULL)

fre<-data.frame()
sum=0

for (i in seq(400)) {
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
