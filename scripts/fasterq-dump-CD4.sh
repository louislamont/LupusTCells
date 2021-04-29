#!/usr/bin/env bash

sed -e 1d ../data/CD4/SraRunTable.txt | while read line
do

SRR=$(echo ${line} | cut -d , -f1)

outdir="../data/CD4/sequencing"
[[ -d ${outdir} ]] || mkdir ${outdir}

echo "SRR: ${SRR}"

echo "prefetch -O ${outdir} ${SRR}"

prefetch -O ${outdir} ${SRR}

echo "fasterq-dump ${outdir}/${SRR}/${SRR}.sra --split-3 --progress --outdir ${outdir}/${SRR}"

fasterq-dump ${outdir}/${SRR}/${SRR}.sra --split-3 --progress --outdir ${outdir}/${SRR}

mv ${outdir}/${SRR}/${SRR}.sra_1.fastq ${outdir}/${SRR}/${SRR}_1.fastq
mv ${outdir}/${SRR}/${SRR}.sra_2.fastq ${outdir}/${SRR}/${SRR}_2.fastq

echo "Compressing"

pigz ${outdir}/${SRR}/${SRR}_1.fastq ${outdir}/${SRR}/${SRR}_2.fastq

echo "${SRR} done"

done
