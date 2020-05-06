import csv
import re
seq={}
codon={}
dict=[]
list=[]
list2=[]
file=open("D:\\研究生\\密码子\\rna_tissue_consensus.tsv\\lung_test.txt",'r+')
gene=open("D:/研究生/密码子/human_cell_lines/CDS_DNA.fa",'r+')
headers=['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
for line in gene:
    line=line.strip('\n')
    gr=re.match(r">(.+)_NM_\d*",line)
    if gr:
        gn=gr.group(1)
    else:
        seq[gn]=line
with open("D:/研究生/密码子/rna_tissue_consensus.tsv/lungcount_test.csv", 'w+',newline='') as count:
    writer=csv.writer(count)
    writer.writerow(headers)
for line in file:
    line=line.strip('\n')
    list=line.split(' ')
    for i in list:
        try:
            if seq[i]:
                for n in range(100,len(seq[i])-1,3):
                    cod=seq[i][n:n+3]
                    codon.setdefault(cod,0)
                    codon[cod]+=1
        except:
            continue
    with open("D:/研究生/密码子/rna_tissue_consensus.tsv/lungcount_test.csv",'a',newline='') as count:
        writer=csv.writer(count)
        dict= sorted(codon.items(), key=lambda x: x[0], reverse=False)
        for m in range(len(dict)):
           str(dict[m][1])
           list2.append(dict[m][1])
        writer.writerow(list2)
    list2.clear()
    codon.clear()
gene.close()
file.close()
count.close()
