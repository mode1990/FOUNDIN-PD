##### Integrative Multiomics Bimodal Data Analysis #####
# by Mo Dehestani, March 2025

#-----------------------------#
# 0. Load Libraries
#-----------------------------#

# Load necessary libraries (with version information)
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

#-----------------------------#
# 1. scRNA-seq Processing
#-----------------------------#

# 1.1 Load and QC Single scRNA-seq Sample
rna_counts <- Read10X("/path/to/the/sample#/filtered_feature_bc_matrix/")$`Gene Expression`
pbmc_rna <- CreateSeuratObject(counts = rna_counts)
pbmc_rna[["percent.mt"]] <- PercentageFeatureSet(pbmc_rna, pattern = "^MT-")

# Pre-QC plot
VlnPlot(pbmc_rna, features = c("nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

# Filter cells based on QC metrics
pbmc_rna <- subset(pbmc_rna, subset = nCount_RNA < 25000 & nCount_RNA > 1000 & percent.mt < 20)

# Post-QC plot
VlnPlot(pbmc_rna, features = c("nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

# Save QC-ed object
saveRDS(pbmc_rna, file = "/path/to/the/output/s#.rds")

# 1.2 Merge Multiple scRNA-seq Samples and Remove Doublets
rna_samples <- list.files(path = "/path/to/the/output/s#", pattern = "*.rds", full.names = TRUE)
rna_objects <- lapply(rna_samples, readRDS)

# Add sample metadata
sample_names <- paste0("Sample", seq_along(rna_objects))
for (i in seq_along(rna_objects)) {
  rna_objects[[i]]@meta.data$Sample <- sample_names[i]
}

# Doublet removal
filtered_rna_objects <- lapply(rna_objects, function(seurat_obj) {
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$scDblFinder_class <- sce$scDblFinder.class
  subset(seurat_obj, subset = scDblFinder_class == "singlet")
})

# Merge filtered objects and save
merged_rna <- do.call(merge, filtered_rna_objects)
saveRDS(merged_rna, file = "path/to/the/merged_GEX_QCied_DblRemoved.rds")

#-----------------------------#
# 2. scATAC-seq Processing
#-----------------------------#

# 2.1 Load and QC Single scATAC-seq Sample
atac_counts <- Read10X("/path/to/the/sample#/filtered_feature_bc_matrix/")$Peaks

# Create Seurat object for ATAC data
pbmc_atac <- CreateSeuratObject(counts = atac_counts)

# Add ATAC-seq data with standard chromosomes
grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
atac_counts <- atac_counts[as.vector(grange_use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
genome(annotations) <- "hg38"
frag_file <- "/path/to/the/sample/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"), fragments = frag_file, min.cells = 10, annotation = annotations)
pbmc_atac[["ATAC"]] <- chrom_assay

#TSS and Nucleosome QC
pbmc_atac <- NucleosomeSignal(pbmc_atac)
pbmc_atac <- TSSEnrichment(pbmc_atac)

#TSS and Nucleosome plot
DensityScatter(pbmc_atac, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Pre-QC plot
VlnPlot(
  object = pbmc_atac,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
# Filter cells based on QC metrics
pbmc_atac <- subset(pbmc_atac, subset = nCount_ATAC > 2000 & nCount_ATAC < 30000 & 
    nucleosome_signal < 2 & TSS.enrichment > 1)

pbmc
# Post-QC plot
VlnPlot(
  object = pbmc_atac,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
# Save QC-ed ATAC object
saveRDS(pbmc_atac, file = "/path/to/the/output/s#.rds")

# 2.2 Merge Multiple scATAC-seq Samples and Remove Doublets
atac_samples <- list.files(path = "/path/to/the/output/s#", pattern = "*.rds", full.names = TRUE)
atac_objects <- lapply(atac_samples, readRDS)

# Doublet removal
filtered_atac_objects <- lapply(atac_objects, function(seurat_obj) {
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$scDblFinder_class <- sce$scDblFinder.class
  subset(seurat_obj, subset = scDblFinder_class == "singlet")
})

# Merge filtered objects and save
merged_atac <- do.call(merge, filtered_atac_objects)
saveRDS(merged_atac, file = "path/to/the/merged_ATAC_QCied_DblRemoved.rds")

#-----------------------------#
# 3. RNA Data Integration and Clustering
#-----------------------------#

# Load merged RNA object
pbmc_rna <- readRDS("path/to/the/merged_GEX_QCied_DblRemoved.rds")

# Split object by sample
pbmc_rna <- SplitObject(pbmc_rna, split.by = "Sample")

# SCT normalization with mitochondrial regression
pbmc_rna <- lapply(pbmc_rna, function(sample) {
  SCTransform(sample, vars.to.regress = "percent.mt", verbose = TRUE)
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = pbmc_rna, nfeatures = 2000)

# Merge the split objects
pd21 <- pbmc_rna$Sample21
pbmc_rna$Sample21 <- NULL
pbmc_rna <- merge(pd21, y = pbmc_rna, project = "pd", merge.data = TRUE)

# Run PCA and determine the number of PCs to use
pbmc_rna <- RunPCA(object = pbmc_rna, npcs = 30, assay = "SCT", features = features)
DimHeatmap(pbmc_rna, dims = 1:30, cells = 1000, balanced = TRUE)
ElbowPlot(pbmc_rna)

# Sample-wise Harmony integration
pbmc_rna <- IntegrateLayers(
  object = pbmc_rna, method = HarmonyIntegration,
  orig.reduction = "pca", layers = "Sample", new.reduction = "harmony",
  verbose = FALSE
)

# KNN graph and clustering
pbmc_rna <- RunUMAP(object = pbmc_rna, assay = "SCT", reduction = "harmony", dims = 1:30)
pbmc_rna <- FindNeighbors(object = pbmc_rna, assay = "SCT", reduction = "harmony", graph.name = "test", dims = 1:30)
pbmc_rna <- FindClusters(object = pbmc_rna, graph.name = "test", resolution = 0.3)

# Inspect and annotate clusters
DimPlot(pbmc_rna, group.by = "Sample")
DimPlot(pbmc_rna, label = TRUE)

# Marker gene visualization for cell type annotation
DefaultAssay(pbmc_rna) <- "SCT"
feature_list <- list(
  EarlyNeuronProgenitor = c("RFX4", "HES1", "SLIT2"),
  LateNeuronProgenitor = c("DLK1", "LGALS1", "VCAN"),
  DopaminergicNeurons = c("TH", "ATP1A3", "ZCCHC12", "MAP2", "SYT1"),
  ImmatureDopaminergicNeurons = c("TPH1", "SLC18A1", "SLC18A2", "SNAP25"),
  ProliferatingFloorPlateProgenitors = c("HMGB2", "TOP2A", "MKI67")
)
for (name in names(feature_list)) {
  FeaturePlot(pbmc_rna, features = feature_list[[name]], label = TRUE)
}

saveRDS(pbmc_rna, file = "/path/to/the/merged_GEX_integrated_annotated_final.rds")

#-----------------------------#
# 4. ATAC Data Processing and Integration
#-----------------------------#

# Load merged ATAC object
pbmc_atac <- readRDS("path/to/the/merged_ATAC_QCied_DblRemoved.rds")

# 4.1 ATAC Processing and Dimensionality Reduction
DefaultAssay(pbmc_atac) <- "peaks"
pbmc_atac <- AnnotateChromatinPeaks(object = pbmc_atac, annotation = annotations)
pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = "q0")
pbmc_atac <- RunTFIDF(pbmc_atac)
pbmc_atac <- RunSVD(pbmc_atac)

# UMAP and clustering for ATAC data
pbmc_atac <- RunUMAP(pbmc_atac, reduction = "lsi", dims = 2:30)
pbmc_atac <- FindNeighbors(pbmc_atac, reduction = "lsi", dims = 2:30)
pbmc_atac <- FindClusters(pbmc_atac, verbose = FALSE, resolution = 0.5)

# 4.2 RNA to ATAC Label Transfer
transfer_anchors <- FindTransferAnchors(reference = pbmc_rna, query = pbmc_atac, reduction = "cca")
pbmc_atac <- TransferData(anchorset = transfer_anchors, refdata = pbmc_rna$celltype, weight.reduction = "lsi")

# 4.3 MACS2 Peak Calling
# Create a Python environment for MACS2
reticulate::conda_create("macs2_env")
# Activate the environment and install MACS2
reticulate::use_condaenv("macs2_env")
reticulate::py_install("MACS2", pip = TRUE)
# Check Python configuration
reticulate::py_config()
Sys.which("python")
system("which macs2")
# Set up environment with explicit paths
use_condaenv("macs2_env", required = TRUE)
# Set Python path explicitly before calling peaks
python_path <- reticulate::conda_python("macs2_env")
Sys.setenv(RETICULATE_PYTHON = python_path)
Sys.setenv(PATH = paste(dirname(python_path), Sys.getenv("PATH"), sep=":"))
# Run peak calling
peaks <- CallPeaks(
  object = pbmc_atac,
  group.by = "Mutation"
)
saveRDS(peaks, "/path/to/macs2_peaks.rds")

# 4.4 Differentially Accessible Peaks Analysis
da_peaks <- FindMarkers(
  object = pbmc_atac,
  ident.1 = 'SNCA', #LRRK2 or GBA1
  ident.2 = "HC",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)
# Select top peaks given pvalue and pct
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.05 & da_peaks$pct.1 > 0.2, ])
# Find motifs 
enriched.motifs <- FindMotifs(object = pbmc_atac, features = top.da.peak)
# Plot the motif 
MotifPlot(object = pbmc_atac, motifs = head(rownames(enriched.motifs)))
# Save motifs
saveRDS(enriched.motifs, "/path/to/enriched_motifs.rds")

# 4.5 Link Peaks to Genes
# First compute the GC content for each peak
gene.activities <- GeneActivity(pbmc_atac)
DefaultAssay(pbmc_atac) <- 'ATAC'
pbmc_atac <- RegionStats(pbmc_atac, genome = BSgenome.Hsapiens.UCSC.hg38)
pbmc_atac <- LinkPeaks(
  object = pbmc_atac,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = "TH", #as example
  peak.slot = "counts",
  expression.slot = "data",
)
# Visualize peak tracks 
CoveragePlot(
  object = pbmc_atac,
  group.by = "Mutation",
  region = "TH",
  peaks = macs2peaks
)

#-----------------------------#
# 5. Multiomics Integration with MOFA2
#-----------------------------#

library(Seurat)
library(Signac)
library(MOFA2)

# 5.1 Load and Prepare Data
rna <- readRDS("/path/to/rna.rds")
atac <- readRDS("/path/to/atac.rds")
# Check data consistency
summary(rna@meta.data$nCount_RNA)
summary(atac@meta.data$nCount_ATAC)

# 5.2 Subset to Dopaminergic Neurons
Idents(pbmc_rna) <- "BroadCellType"
rna <- subset(pbmc_rna, idents = c("DAN"))

Idents(pbmc_atac) <- "BroadCellType"
atac <- subset(pbmc_atac, idents = c("DAN"))

# 5.3 Extract and Normalize Data
rna_counts <- GetAssayData(rna, slot = "data")  # RNA expression (normalized)
atac_counts <- GetAssayData(atac, slot = "data")  # ATAC peaks (normalized)

# Save matrices to CSV
write.csv(as.data.frame(rna_counts), "/path/to/rna_expression.csv", row.names = TRUE)
write.csv(as.data.frame(atac_counts), "/path/to/atac_expression.csv", row.names = TRUE)

# 5.4 Process Cell Barcodes
colnames(rna_counts) <- gsub("-1_.*", "", colnames(rna_counts))  # Remove "-1_XX_XX"
colnames(atac_counts) <- gsub("-1_.*", "", colnames(atac_counts))  # Remove "-1_XX"

# Check barcodes
head(colnames(rna_counts))
head(colnames(atac_counts))

# 5.5 Match Cells Across Modalities
common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
rna_counts <- rna_counts[, common_cells]
atac_counts <- atac_counts[, common_cells]

# 5.6 Create and Run MOFA2 Model
data_list <- list(
  RNA = as.matrix(rna_counts),
  ATAC = as.matrix(atac_counts)
)

mofa_obj <- create_mofa(data_list)

# Set training options
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 15

mofa_obj <- prepare_mofa(mofa_obj,
                         model_options = model_opts
)

# Run MOFA in python
mofa_obj <- run_mofa(mofa_obj, use_basilisk = TRUE)

saveRDS(mofa_obj, "/path/to/mofa_model.rds")

# 5.7 Add Metadata and Visualize Results
rownames(rna@meta.data) <- make.unique(gsub("-1_.*", "", rownames(rna@meta.data)))
aligned_metadata <- rna@meta.data[unlist(samples_names(mofa_obj)), , drop = FALSE]
samples_metadata(mofa_obj) <- aligned_metadata %>%
  tibble::rownames_to_column("sample") %>%
  as.data.table

rownames(atac@meta.data) <- make.unique(gsub("-1_.*", "", rownames(atac@meta.data)))
aligned_atac_metadata <- atac@meta.data[unlist(samples_names(mofa_obj)), , drop = FALSE]
samples_metadata(mofa_obj) <- cbind(
  samples_metadata(mofa_obj),
  aligned_atac_metadata
)

# Plot data overview
plot_data_overview(mofa_obj)

# Plot factors 
plot_factor(mofa_obj, 
            factor = 1:3,
            color_by = "Mutation",
            shape_by = "BroadCellType"
)

# Plot variance explained and factor loadings
plot_variance_explained(mofa_obj, x="view", y="factor")
