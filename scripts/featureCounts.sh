#!/usr/bin/env bash

### Set in and out directories
indir="../output/CD8/bams"
outdir="../output/CD8/counts"
### Get bam files in a variable
bamfiles=$(ls ${indir}/*.bam)
#echo ${bamfiles}

###### Run featureCounts
call="featureCounts -p -a ../data/gencode.v37.annotation.gtf.gz \
-s 2 -o ${outdir}/counts.txt -T 6 ${bamfiles}"

echo $call
eval $call
