## Notch3 regulates pericyte phenotypic plasticity in colorectal cancer

This repository contains the code used to create all the plots associated with our manuscript: **Notch3 regulates pericyte phenotypic plasticity in colorectal cancer** 

In-house raw data from this work are deposited in GEO with the following accession numbers: 

**scRNA-seq**
  - Control colon: **GSE296873**
  - AOM/DSS-induced tumor: **GSE296872**

## Methods Summary

### Preprocessing

-   Raw scRNA FASTQ files were processed using **10X Genomics CellRanger v7.1.0** (`cellranger count`, default parameters)
-   Mouse reference transcriptome: `mm10` 

### Downstream Analysis

-   scRNA-seq QC, Harmony integration, clustering, plotting: **Seurat v4**

<sub>For more details on our methodology, see the corresponding manuscript.:...