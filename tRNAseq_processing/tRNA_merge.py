import re
import csv
anti={}
c1=open("D:\\研究生\\tRNAseq\\GSM16248\\GSM1624821_TP4_TGACCA_L005_R1_001.txt",'r+')
for line in c1:
    line=line.strip('\n')
    s=re.match(r".*-.*(\w\w\w)\s(\d*)",line)
    if s:
        codon=s.group(1)
        num=s.group(2)
        num=eval(num)
    else:
        continue
    anti.setdefault(codon,0)
    if anti[codon]:

        anti[codon]=anti[codon]+num
    else:

        anti[codon]=num
with open("D:\\研究生\\tRNAseq\\GSM16248\\merge21.csv",'w+',newline='') as out1:
    csvwriter=csv.writer(out1)
    csvwriter.writerow(['anti_codon','num'])
    Santi=sorted(anti.items(),key=lambda x:x[1],reverse=True)
    for a in range(len(Santi)):
        csvwriter.writerow([Santi[a][0],Santi[a][1]])
c1.close()
out1.close()