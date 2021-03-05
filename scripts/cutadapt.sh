#!/usr/bin/env bash

sed -e 1d ../data/CD8/SraRunTable.txt | while read line
do

SRR=$(echo ${line} | cut -d , -f1)

indir="../data/CD8/sequencing/${SRR}"
outdir="../output/CD8/trimmed/${SRR}"
[[ -d ${outdir} ]] || mkdir ${outdir}

echo "SRR: ${SRR}"

call="cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-j 8 --trim-n -m 25 -q 20,20 \
-o ${outdir}/${SRR}_1.fastq.gz \
-p ${outdir}/${SRR}_2.fastq.gz \
${indir}/${SRR}_1.fastq.gz \
${indir}/${SRR}_2.fastq.gz"

echo ${call}
eval ${call} 2>&1 | tee ../output/CD8/trimmed/log/${SRR}.log

echo "${SRR} done"

done
