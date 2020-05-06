#!/usr/bin/perl



open (GENE,"D:/研究生/密码子/human_cell_lines/CDS_DNA.fa") || die "cannot open:$!";
open (HIGH,"D:/研究生/密码子/human_cell_lines/high.csv") || die "cannot open:$!\n";
open (COUNT_HIGH,">D:/研究生/密码子/human_cell_lines/count_high.txt") || die "cannot open:$!\n";

my %seq;
my $gene;
my %codon;
while (<GENE>){
chomp;
if(/>(.+)_NM_\d*/){
$gene=$1;
next;
}
$seq{$gene}.= $_;
}
close GENE;

while (<HIGH>){
chomp;
countCodon($_,\%seq);
}


sub countCodon{
my ($high_gene,$cDNA)=@_;
my $cod;
if (exists $$cDNA{$high_gene}){
for ($i=100;$i<=length($$cDNA{$high_gene})-103;$i+=3){
$cod=substr($$cDNA{$high_gene},$i,3);
$codon{$cod}++;
}
}
}
print %codon;
while (($key,$value)=each %codon){
print COUNT_HIGH "$key $value\n";
}
close HIGH;
close COUNT_HIGH;