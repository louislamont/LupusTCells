# Gene expression in Lupus CD8 T cells (and possibly later, CD4 T cells)

This repository contains an analysis of some publically available data, including scripts for data download, clean up, alignment, analysis, and reporting (mostly shell and R in the form of R markdown files). 

Sequencing files were downloaded from SRA, accessions GSE97264 (CD8 T cells) and GSE97263 (CD4 T cells). Thank you to the M. Botto and C. Shen for making their data available.

# Summary

Accessions contain RNA-Seq data from 3 lupus disease groups: Active (AL), Less Active (LAL), and Control (Ctrl). Differential expression was performed between each of the groups (i.e., AL v Ctrl, LAL v Ctrl, AL v LAL), as well as downstream analyses such as co-expression, functional enrichment, and network analysis.

### Methods

*Data acquisition*

Five SRA files from each disease group (Control, Less Active, Active) were downloaded using `prefetch` and converted into fastq files using `fasterq-dump` from the SRA toolkit .

*Sequence data processing*

FastQC (v0.11.9) was used to assess quality metrics of sequencing data. Reads were trimmed of Illumina universal adapters using cutadapt (v3.2) and the following flags: -m 25 --trim-n -q 20,20. Trimmed reads were aligned to the GENCODE v37 human genome (GRCh38.p13) using hisat2 (v2.2.1) and the following flags: --rna-strandedness RF --dta. Multimapping reads were removed with samtools view (-q 15) and reads were sorted with samtools sort (v1.10). Read counts were determined using featureCounts (v2.0.1) with the flags: -p -s 2.

*Data analysis*

Counts were read into R (v4.0.4) and normalized using the TMM method with the edgeR package (v3.32.1). edgeR was also used for differential expression between the three disease groups. Significantly differentially expressed genes were determined as having FDR < 0.05 and abs(log2(FC)) >= 2. Normalized counts were exported using the cpm() function. Co-expression analysis was performed on normalized counts using WGCNA (v1.70.3) with minModuleSize=30, networkType="signed hybrid", corType="bicor", mergeCutHeight = 0.25, and maxPOutliers=0.1. Modules enriched in DE genes were determined using a hypergeometric test. Functional enrichment was performed using WebGestaltR (v0.4.4) for over-representation analyses, clusterProfiler (v3.18.1) and ReactomePA (v1.34.0) for GSEA, and [Enrichr](https://maayanlab.cloud/Enrichr/). Gene networks were visualized in Cytoscape (v3.8.1).

### Results

Currently, results are only for CD8 T cells.

**Differential Expression**

Differential expression at the gene level was examined between LAL and Ctrl, AL and Ctrl, and AL and LAL. The results are summarized in table 1. Genes were determined to be DE if they had an FDR less than 0.05 and a fold change greater than 2 (either up or down).

**Table 1: DE genes in all categories (FDR<0.05, abs(FC)>=2)**

|               | LAL v Ctrl | AL v Ctrl | AL v LAL |
|---------------|------------|-----------|----------|
| Upregulated   | 10         | 504       | 109      |
| Downregulated | 33         | 459       | 8        |

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/compareDE/all-groups.png" alt="Fig 1a: Overlap of DE genes in LAL v Ctrl, AL v Ctrl, and AL v LAL" width="600"/>

**Fig 1a: Overlap of DE genes in LAL v Ctrl, AL v Ctrl, and AL v LAL**

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/compareDE/FC.png" alt="Fig 1b: Comparison of gene-level fold changes between AL v Ctrl and LAL v Ctrl" width="600"/>

**Fig 1b: Comparison of gene-level fold changes between AL v Ctrl and LAL v Ctrl**

43 genes are DE between Less Active Lupus (LAL) and Control (Ctrl) groups (Table 1). Of these, 37 are also DE between Active Lupus (AL) and Ctrl (Fig 1a). All genes DE in LAL v Ctrl are expressed in the same direction as in AL v Ctrl (Fig 1b).

On the other hand, there are 963 genes DE between AL and Ctrl (Table 1, Fig 1a). Of note, there are a handful of genes DE in AL v Ctrl that change in the opposite direction as LAL v Ctrl, including LTF, ITGAD, and IGHG1 (upregulated in AL), and UTS2 and ZFP57 (downregulated in AL) (Fig 1b).

Finally, I found 117 genes DE between AL and LAL, with a strong majority upregulated in AL (Table 1). Of these, nearly all are also DE with respect to the AL v Ctrl analysis (Fig 1a).

**Co-expression**

Genes were grouped into 32 co-expression modules (plus the "grey" uncategorized module). Of these, 4 modules were enriched in genes DE either in Active or Less Active disease compared to control, or DE between Active and Less Active disease groups.

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/wgcna/modules-w-DE-genes.png" alt="Fig 2: Association with disease group and modules with enriched genes" width="600"/>

**Fig 2: Association of module expression with disease state. Number of DE genes in each module and enrichment FDR indicated**

Here, we can identify co-expression modules related to lupus (enriched in DE genes in both lupus groups compared to control), modules related to active lupus (enriched in DE genes in AL, but not LAL compared to control), and modules potentially related to lupus activity (modules enriched in DE genes between AL and LAL).

For Less Active lupus vs. control, we see enrichment in the skyblue (5/52, FDR=1e-05) and turquoise (24/3652, FDR=1e-04) modules. Both of these modules are downregulated compared to controls. For Active lupus vs. control, we also see enrichment in the skyblue (17/52, FDR=1e-07) and turquoise (322/3652, FDR=5e-11) modules, as well as the blue (385/2649, FDR<2.2e-16) and darkgrey (11/67, FDR=0.02) modules. Of these, blue is upregulated compared to control, while the other three modules are downregulated. Strikingly, when examining genes DE between AL and LAL groups, most (104/117 DE genes, 89%) fall into the blue module (2649 total blue genes, FDR<2.2e-16). All 104 genes are upregulated in AL compared to LAL. Expression and co-expression module of all DE genes is shown in Fig 3.

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/heatmaps/allDE-genes.png" alt="Fig 3: Heatmap of expression of DE genes, with co-expression module indicated " width="600"/>

**Fig 3: Expression and module membership of differentially expressed genes in lupus CD8 T cells**

(note: Small modules are interesting, look into skyblue and darkgrey)

**Functional analysis**

Here, I performed GSEA to examine whole transcriptome differences in expression between groups.

*Active lupus vs. controls*

The top 30 enriched Reactome pathways in AL vs. Ctrl are shown in Fig 4a and include many related cell cycle and cytoskeleton pathways (cell cycle checkpoints, mitotic metaphase and anaphase, RHO GTPase effectors), interferon pathways, neutrophil degranulation, and hemostatis. Although not shown, these pathways are all upregulated with respect to control samples.

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/gsea/Reactome/Active-v-Ctrl-pathways.png" alt="Fig 4a: GSEA enrichment of Reactome pathways between AL and Ctrl patients" width="600"/>

**Fig 4a: GSEA enrichment of Reactome pathways between AL and Ctrl patients**

*Active lupus vs. less active lupus*

Likewise, in AL vs. LAL, we see a similar set of pathways related to cell cycle (Fig 4b) which are also upregulated in AL with respect to LAL patients. However, we also observe a set of pathways related to translation (Eukaryotic translation elongation, peptide chain elongation), with genes in these pathways downregulated in AL with respect to LAL patients.

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/gsea/Reactome/Active-v-LessActive-pathways.png" alt="Fig 4b: GSEA enrichment of Reactome pathways between AL and LAL patients" width="600"/>

**Fig 4b: GSEA enrichment of Reactome pathways between AL and LAL patients**

**Network analysis**

It was suprising that nearly 90% of genes DE between AL and LAL fell into the blue module, so I wanted to examine these further. These 104 genes were all upregulated in 102 out of 104 genes were annotated in the STRING database; strikingly, nearly all of these were known to interact in a very well connected network (Fig 5).

<img src="https://github.com/louislamont/LupusTCells/blob/main/plots/CD8/cytoscape/AL-v-LAL-DE-Blue-w-legend.png" alt="Fig 5: Known protein-protein interactions of DE genes in the blue co-expression network" width="600"/>

**Fig 5: Known protein-protein interactions of DE genes in the blue co-expression network**

These genes are also enriched for many GO terms relating to the cell cycle, including cell cycle phase transition (FDR<2.2e-16), positive regulation of cell cycle process (FDR<2.2e-16), cell cycle checkpoint (FDR<2.2e-16), and regulation of cell cycle G2/M phase transition (FDR=4.4e-14), as well as many similar Reactome pathways, such as Amplification of signal from the kinetochores (FDR=2.1e-12) and Activation of E2F1 target genes at G1/S (FDR=1.7e-9).

Additionally, Enrichr found enriched binding near the TSS of these genes for many TFs, including E2F4 (FDR=6.9e-68), FOXM1 (FDR=1.1e-33), NFYA (FDR=4.1e-21), NFYB (FDR=1.9e-18), and SIN3A (FDR=4.9e-16). IRF3 also shows increased binding near these genes (FDR=5.2e-16). Enrichr also predicted enriched interactions between genes in this group and hsa-miR-193b-3p (FDR=9.2e-35).

### Thoughts and discussion

* LAL shows kind of an "inbetween" phenotype with regards to differentially expressed genes. 2 AL patients have a set of genes in the blue co-expression module which look more like LAL/Ctrl (i.e. low expression vs. the other 3 patients). Could this be related to disease phenotype? Maybe the other 3 patients have higher disease activity, nephritis, etc. Can see if there is any metadata on this.

LAL and AL: 

* Cytoskeleton and cell cycle related terms are upregulated in active lupus vs. less active lupus - increased activity of CD 8 T cells contibutes to increased lupus activity.

* [miR-193b "within the peripheral mononuclear cells \[was\] capable of differentiating between lupus patients with nephritis and those without nephritis"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6021342/)

