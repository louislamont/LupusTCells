#!/usr/bin/env bash

lib="/Users/matt/projects/hisat-indexes/gencode/gencode.v37"

sed -e 1,2d ../data/CD8/SraRunTable.txt | while read line
do

SRR=$(echo ${line} | cut -d , -f1)

indir="../output/CD8/trimmed/${SRR}"
outdir="../output/CD8/bams"
[[ -d ${outdir} ]] || mkdir ${outdir}

echo "SRR: ${SRR}"

#### Run HISAT2
### Indicate strandedness (RF)
### Use 20 cores
### Output any unaligned reads
### Stringtie needs --dta
### Pipe directly to samtools view for quality filtering (-q 15) and
### sorting
call="hisat2 -p 6 --rna-strandness RF \
--dta -x ${lib} \
--summary-file ../output/CD8/bams/logs/${SRR}.log \
-1 ${indir}/${SRR}_1.fastq.gz \
-2 ${indir}/${SRR}_2.fastq.gz | \
samtools view -q 15 -Su | samtools sort -m 1G -@ 4 - \
-o ${outdir}/${SRR}.sorted.bam"

echo $call
eval $call

echo "${SRR} done"

done
