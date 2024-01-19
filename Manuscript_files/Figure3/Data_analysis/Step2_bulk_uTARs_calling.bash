#!/bin/bash
#nohup bash scripts/SingleCellHMM_MW_SS2.bash bam/SRR7309700Aligned.sortedByCoord.out.bam gencode.vM27.primary_assembly.annotation.gtf 50 scripts/

INPUT_BAM=$1 #pbmc4k_possorted_genome_bam.bam
refFlat=gencode.v43.chr_patch_hapl_scaff.annotation.refflat.refFlat #gene annotations in refFlat format
CORE=$3
MERGEBP=$4
THRESH=$5
PL=$6

CORE="${CORE:-5}"
#MINCOV="${MINCOV:-5}"
MERGEBP="${MERGEBP:-500}"
THRESH="${THRESH:-10000000}"
CURDIR=`pwd`
PL="${PL:-${CURDIR}/scripts}"

reads=`samtools view -q 255 $INPUT_BAM | wc -l`
echo "Number of aligned reads is $reads"
minCovReads=`expr $reads / ${THRESH}`
MINCOV=$minCovReads

PREFIX=`echo ${INPUT_BAM} | rev | cut -d / -f 1 |cut -d . -f 2- |rev`
tmp="HMM_features"
TMPDIR=${PREFIX}_${tmp}
mkdir ${TMPDIR}

exec > >(tee SingleCellHMM_Run_${TMPDIR}.log)
exec 2>&1
echo "Path to SingleCellHMM.R   $PL" 
echo "INPUT_BAM                 $INPUT_BAM"
echo "temp folder               $TMPDIR"
echo "number Of thread          $CORE"
echo "minimum coverage          $MINCOV"
echo "thresholded at 1 in $THRESH reads"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Spliting and sorting reads..."
#bedtools bamtobed -i ${INPUT_BAM} -split |LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk '$1!="." && !/GL|JH|MU/ {print $0}' | gzip > ${TMPDIR}/${PREFIX}_split.sorted.bed.gz 
bedtools bamtobed -i ${INPUT_BAM} -split -bedpe |cut -f1,2,6,7,8,9 |LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk '$1!="\." && !/GL|JH|MU/ {print $0}' | gzip > ${TMPDIR}/${PREFIX}_split.sorted.bed.gz

cp scripts/SingleCellHMM.R ${TMPDIR}/
cd ${TMPDIR}
zcat ${PREFIX}_split.sorted.bed.gz  |awk '{print $0 >> "chr"$1".bed"}'
cd /scratch/user/s4716765/TAR-scRNA-seq/from_fastqs/${TMPDIR}
for f in chr*.bed; do  R --vanilla --slave --args ./ ${f}  < SingleCellHMM.R  > ${f}.log 2>&1 & pids+=($!); done
