# Integrative Multiomics for FOUNDIN-PD


## Overview

This repository contains the R code for an integrative multiomics bimodal data analysis pipeline, focusing on single-cell RNA sequencing (scRNA-seq) and single-cell ATAC sequencing (scATAC-seq) data integration. The analysis is designed for neuronal cell types with special focus on dopaminergic neurons under different mutation conditions (SNCA, LRRK2, GBA1) compared to healthy controls.

## Features

- Complete scRNA-seq and scATAC-seq data processing pipeline
- Quality control and doublet removal for both modalities
- Integration and batch correction using Harmony
- Peak calling with MACS2
- Differential accessibility analysis and motif enrichment
- Peak-to-gene linkage visualization
- Multi-omics factor analysis using MOFA2

## Dependencies

### R Libraries

```r
library(biovizBase)              # Version: 1.50.0
library(irlba)                   # Version: 2.3.5.1
library(Matrix)                  # Version: 1.6-5
library(BSgenome.Hsapiens.UCSC.hg38) # Version: 1.4.5
library(BSgenome)                # Version: 1.70.2
library(rtracklayer)             # Version: 1.62.0
library(BiocIO)                  # Version: 1.12.0
library(Biostrings)              # Version: 2.70.3
library(XVector)                 # Version: 0.42.0
library(EnsDb.Hsapiens.v86)      # Version: 2.99.0
library(ensembldb)               # Version: 2.26.0
library(AnnotationFilter)        # Version: 1.26.0
library(GenomicFeatures)         # Version: 1.54.4
library(AnnotationDbi)           # Version: 1.64.1
library(BiocGenerics)            # Version: 0.48.1
library(fastmap)                 # Version: 1.2.0
library(Signac)                  # Version: 1.13.0
library(dplyr)                   # Version: 1.1.4
library(lifecycle)               # Version: 1.0.4
library(Seurat)                  # Version: 5.0.3
library(SeuratObject)            # Version: 5.0.1
library(sp)                      # Version: 2.1-3
library(MOFA2)                   # Required for multi-omics integration
```

### Python Dependencies

- MACS2 (for peak calling)
- reticulate (for R-Python interface)

## Workflow Overview

1. **scRNA-seq Processing**
   - Quality control and filtering
   - Doublet removal
   - Integration of multiple samples

2. **scATAC-seq Processing**
   - Quality control and filtering
   - Doublet removal
   - Integration of multiple samples

3. **RNA Data Integration and Clustering**
   - SCTransform normalization
   - Harmony batch correction
   - Clustering and visualization
   - Cell type annotation

4. **ATAC Data Processing and Integration**
   - Peak annotation and dimensionality reduction
   - RNA to ATAC label transfer
   - MACS2 peak calling
   - Differential accessibility analysis
   - Motif enrichment analysis
   - Peak-to-gene linkage

5. **Multiomics Integration with MOFA2**
   - Data preparation and cell matching
   - MOFA2 model creation and training
   - Factor visualization and interpretation


