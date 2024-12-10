################################################################################
##### ANNEX 1: Initial data processing (Burdziak et all. dataset)          #####
################################################################################

### GENERAL DATA INFORMATION ###
# Burdziak et al. (Science, 2023)
# Progression Cohort: single-cell RNA seq data
# Accession: GSE207938

# Design: 	Characterization of single-cell transcriptional profiles of 
# lineage-traced pancreatic epithelial cells (mKate2+) isolated directly 
# by FACS-sorting from normal (N1), regenerating (N2), pre-neoplastic 
# (K1, K2), bening neoplastic (PanIN, K3, K4), malignant (PDAC, K5) and 
# metastasis (K6) murine pancreatic tissues

wd <- r"(C:\Users\javim\OneDrive - UNIVERSIDAD DE SEVILLA\TFM_UPO)"
setwd(wd)

### Libraries
library(tidyverse)
library(glue)
library(Seurat)    # scRNAseq analysis
library(schard)    # h5ad file transformation
library(infercnv)  # CNV analysis
library(patchwork) # plot tools

# ==============================================================================
### Data import and H5AD transformation to seurat
# Raw counts and processed data available in file: GSE207938_ProgressionCohort.h5ad.gz
# Data manually downloaded and unzipped from:
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207938&format=file&file=GSE207938%5FProgressionCohort%2Eh5ad%2Egz
h5ad.path <- file.path(wd,"/data/GSE207938_ProgressionCohort.h5ad")
# Transforming anndata raw data to seurat (w/ schard)
# When importing data from the h5ad file, both counts and data layers will have 
# the same data. I import two copies and save the normalised data to add it later
burdziak_seu <- h5ad2seurat(h5ad.path, use.raw = TRUE)   
burdziak_seu

# Both counts and data assays contain the same data (raw counts)
identical(burdziak_seu@assays$RNA$counts,burdziak_seu@assays$RNA$data)
head(burdziak_seu@assays$RNA$counts) # raw counts

# Importing a new copy with only norm data
normdata <- h5ad2seurat(h5ad.path, use.raw = FALSE, load.obsm = FALSE)
identical(normdata@assays$RNA$counts,normdata@assays$RNA$data)
head(normdata@assays$RNA$counts) # normalized counts
norm_cts <- normdata@assays$RNA$counts # storing norm counts
rm(normdata)

# Add tissue state in metadata.
unique(burdziak_seu@meta.data$condition) # check conditions
burdziak_seu@meta.data <- burdziak_seu@meta.data %>% mutate(tissue = case_when( # add tissue column, based on reference paper
  condition == "N1" ~ "Normal",
  condition == "N2" ~ "Regeneration",
  condition %in% c("K1", "K1.5", "K2") ~ "Pre-neoplasia",
  condition %in% c("K3", "K4") ~ "PanIN",
  condition == "K5"  ~ "PDAC",
  condition == "K6" ~ "Metastasis"))
burdziak_seu$tissue <- factor(burdziak_seu$tissue, # ordering tissue states
                              levels =c("Normal","Regeneration", "Pre-neoplasia",
                                        "PanIN","PDAC", "Metastasis"))

# ==============================================================================
### Checking data pre-processing
# Uploaded data is already filtered and pre-processed
burdziak_seu[["percent.mt"]] <- PercentageFeatureSet(object = burdziak_seu, pattern = "^MT-")
summary(burdziak_seu@meta.data$percent_mt) # lenght = 0; no MT genes

# QC metrics
glue("Min nCount = {min(burdziak_seu@meta.data$nCount_RNA)}\n")
glue("Min nFeat = {min(burdziak_seu@meta.data$nFeature_RNA)}\n")
VlnPlot(burdziak_seu, features = c("nCount_RNA", "nFeature_RNA"),
        ncol = 2, pt.size = 0)

# Metadata cols
glue("Metadata:\n")
colnames(burdziak_seu@meta.data)

# Reduced dims
glue("tSNE plot:\n") #tSNE data imported in XEMBED_ dimensional reduction
DimPlot(burdziak_seu, reduction = "XEMBED_", group.by = "cluster", label=TRUE) +
  theme_minimal() + theme(panel.grid = element_blank()) + ggtitle(NULL) + labs(x = "tSNE_1", y = "tSNE_2") 


# ==============================================================================

### Gene symbol formatting

# The dataset originates from mouse tissue However, due to previous treatment, the 
# uploaded format of the genes in the count matrix is not standard (for mice, only
# the first letter is capitalised, e.g: Kras19)

original_genes <- rownames(burdziak_seu) # copy of original gene names 
head(original_genes,10) 
tail(original_genes,10)

# I confirm that the dataset contains converted to uppercase rather than human 
# ortologs. Firts, by checking mouse-specific genes such as Trp53, mouse 
# ortholog of human TP53:
glue("TP53 is not present in the database:", "\n \"TP53\" %in% original_genes = ",
     "TP53" %in% original_genes, "\n\n")
glue("The mouse ortholog, Trp53, is also not present:", "\n \"Trp53\" %in% original_genes = ", 
     "Trp53" %in% original_genes, "\n\n")
glue("However, TRP53 is:", "\n \"TRP53\" %in% original_genes = ", 
     "TRP53" %in% original_genes, "\n\n")

# These genes are easy to process with stringr, using the str_to_sentence() 
# function that maintain only the first letter capitalized. However, the dataset
# also includes genes from the Riken mouse genome encyclopedia project 
# (https://www.sciencedirect.com/science/article/pii/S1631069103002166?via%3Dihub). 
# These symbols are annotated in different ways. While the majority of these genes 
# end with the RIK suffix, some have an additional number (RIKx). Also, there are
# GRIK genes that interfere with the grep functions.
original_genes[grep(x= original_genes, pattern = "RIK")] # checking all gene names that contain RIK

# Genes that may generate conflict when correcting the ID are:
original_genes[grep(x= original_genes, pattern = "RIK\\d")]

# I manually correct the gene symbols by checking the mentioned cases:
corrected_genes <- character(length(original_genes)) # empty vector to save symbols
pattern_grik = "^(GRIK)(\\d)"   # start with GRIK followed by any digit
pattern_riken = "RIK(\\d+)?$"   # ends with RIK followed (or not) by any digit

# Iterate each gene
for (i in seq_along(original_genes)) { # seq_along to use positions
  gene <- original_genes[i]            # assign i symbol to gene var
  
  # Check conditions. First, Grik genes
  if (grepl(pattern = pattern_grik, x = gene)) {
    corrected_genes[i] <- str_to_sentence(gene) # store corrected symbol in position
    
    # Then, correct Riken genes
  } else  if (grepl(pattern = pattern_riken, x = gene)) {
    corrected_genes[i] <- gsub(pattern = "RIK", 
                               replacement= "Rik",   # Correcting only Rik suffix, 
                               x= gene)              # mantaining the rest of uppercases
    # The resting genes are in SYMBOL
  } else {
    corrected_genes[i] <- str_to_sentence(gene)
  }
}

# Checking results
head(corrected_genes,10)
tail(corrected_genes,10)
corrected_genes[grepl(pattern = "Rik", x = corrected_genes, ignore.case = TRUE)]

# Loop to see if all genes are in the correct position (case insensitive)
matches <- logical(length(original_genes))
for (i in seq_along(original_genes)){
  original <- original_genes[i]
  corrected <- corrected_genes[i]
  matches[i] <- grepl(corrected, original, ignore.case = TRUE)
}
any(!matches) # I check for any FALSE match with !. FALSE = all values are match (TRUE)

# ==============================================================================
### Saving results
# I can't rename the features in the seurat object, so I create a new corrected one
# Extracting counts matrix and correcting gene symbols in raw and norm matrices
cts <- LayerData(object= burdziak_seu, layer = "counts") 
identical(rownames(cts), rownames(norm_cts)) # same gene symbols
rownames(cts) <- corrected_genes 
rownames(norm_cts) <- corrected_genes

# Following authors recommendations, transforming unlogged normalized data
# to log2 with 0.1 pseudocount
log2_norm <- log2(norm_cts + 0.1)
log2_norm <- as(log2_norm, "CsparseMatrix") # Seurat integration 

# Creating the corrected Seurat object
burdziak_seu_corrected <- CreateSeuratObject( 
  counts = cts,      # counts
  data = NULL,       # empty norm counts
  meta.data = burdziak_seu@meta.data) # metadata

burdziak_seu_corrected <- SetAssayData(burdziak_seu_corrected, 
             layer = "data",
             new.data = log2_norm,
             assay = "RNA")
burdziak_seu_corrected@reductions[["tsne"]] <- burdziak_seu@reductions$XEMBED_ # dim reductions
burdziak_seu_corrected@reductions[["pca"]] <- burdziak_seu@reductions$XPCA_
burdziak_seu_corrected

# Save RDS for future use
rds <- file.path(wd, "/data/burdziak_seu_corrected.rds")
saveRDS(burdziak_seu_corrected, file = rds)



