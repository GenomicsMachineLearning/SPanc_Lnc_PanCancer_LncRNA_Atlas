cd /nfs_master/prakrithi/lncRNA_atlas/phastcons/
conda activate ucsc
awk '{print $1"\t"$2"\t"$3"\t"$4"\t1\t"$6}' ../Q4386/polyA_FreshFrozen/Visium15_HeadNeck/D1/uTARs_uniq.bed > HND.bed
bigWigAverageOverBed /nfs_master/pallavi/databases/phastCons30way_hg38/hg38.phastCons30way.bw your_regions.bed your_regions_phastcons.tab

#extract fasta for coding potential 
bedtools getfasta -fi ../../test/ST/refdata-gex-GRCh38-2020-A/fasta/genome.fa -bed HNB.bed -s > HNB.fa

# CPAT - Coding potential
nohup ~/miniconda3/envs/cpat_new/bin/cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData -g ../HNB.fa -o HNB_cpat 2> HNB_cpat.log &

awk '$11<0.364 {print }' HNB_cpat.ORF_prob.best.tsv > HNB_noncoding.tsv
#Human coding probability (CP) cutoff: 0.364 (CP >=0.364 indicates coding sequence, CP < 0.364 indicates noncoding sequence) 

conda activate cpat
#RNA fasta for MFE
seqtk seq -r HNB.fa | sed 's/T/U/g'  > HNB_RNA.fa

RNAfold HNB_RNA.fa > HNB_RNAfold.txt
awk '/^>/ {id=$1} / \(-[0-9]+\.[0-9]+\)$/ {mfe=$NF; print id, mfe}' HNB_RNAfold.txt | sed 's/>//g' | sed 's/(//g' | sed 's/)//g' > HNB_MFE.txt


awk -F"_" '{print $1":"$2"-"$3$4}' HNB_phastcons.tab > hnb1
awk '{print $6}' HNB_phastcons.tab > hnb2
paste hnb1 hnb2 > HNB_phastcons_score.txt
rm hnb*


join HNB_phastcons_score.txt HNB_MFE.txt > HNB_phastcons_MFE_data.txt



in R:
d<-read.csv("HNC_phastcons_MFE_data.txt",sep=" ",header=FALSE)
ggplot(d,aes(x=V2,y=V3))+geom_point()

