# Gene expression in Lupus CD4 and CD8 T cells

This is an analysis of some publically available data, including scripts for data download, clean up, alignment, analysis, and reporting.

Sequencing files were downloaded from SRA, accessions GSE97264 (CD8 T cells) and GSE97263 (CD4 T cells). Thank you to the M. Botto and C. Shen for making their data available.

### Methods

*Data acquisition*

Five SRA files from each disease group (Control, Less Active, Active) were downloaded using `prefetch` and converted into fastq files using `fasterq-dump` from the SRA toolkit .

*Sequence data processing*

FastQC (v0.11.9) was used to assess quality metrics of sequencing data. Reads were trimmed of Illumina universal adapters using cutadapt (v3.2) and the following flags: -m 25 --trim-n -q 20,20. Trimmed reads were aligned to the GENCODE v37 human genome (GRCh38.p13) using hisat2 (v2.2.1) and the following flags: --rna-strandedness RF --dta. Multimapping reads were removed with samtools view (-q 15) and reads were sorted with samtools sort (v1.10). Read counts were determined using featureCounts (v2.0.1) with the flags: -p -s 2.

*Data analysis*

Counts were read into R (v4.0.4) and normalized using the TMM method with the edgeR package (v3.32.1). edgeR was also used for differential expression between the three disease groups. Significantly differentially expressed genes were determined as having FDR < 0.05 and abs(log2(FC)) >= 2. Normalized counts were exported using the cpm() function. Co-expression analysis was performed on normalized counts using WGCNA (v1.70.3) with minModuleSize=30, networkType="signed hybrid", corType="bicor", mergeCutHeight = 0.25, and maxPOutliers=0.1. Modules enriched in DE genes were determined using a hypergeometric test. Functional enrichment was performed using WebGestaltR (v0.4.4) for over-representation analyses, clusterProfiler (v3.18.1) and ReactomePA (v1.34.0) for GSEA, and [Enrichr](https://maayanlab.cloud/Enrichr/). Gene networks were visualized in Cytoscape (v3.8.1).

### Results

*Differential Expression*

*Co-expression*

Genes were grouped into 32 co-expression modules (plus the "grey" uncategorized module). Of these, 4 modules were enriched in genes DE either in Active or Less Active disease compared to control, or DE between Active and Less Active disease groups.

![Fig 2: Association with disease group and modules with enriched genes](plots/wgcna/modules-w-DE-genes.png)


