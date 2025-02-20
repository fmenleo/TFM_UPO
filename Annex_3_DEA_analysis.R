################################################################################
##### ANNEX 3: DEA analysis (Burdziak et all. dataset)                     #####
################################################################################

# Working directory
wd <- r"(C:\Users\javim\OneDrive - UNIVERSIDAD DE SEVILLA\TFM_UPO)"
setwd(wd)

# Libraries 
library(Seurat) 
library(tidyverse)
library(scCustomize) #scRNAseq data visualization
library(clusterProfiler) # functional enrichments
library(glue)
library(patchwork) # plots
library(org.Mm.eg.db) # Genome annotation for mouse

# Load processed data (see Annex_1)
burdziak_seu <- LoadSeuratRds("./data/burdziak_seu_corrected.rds")
burdziak_seu

### Standard Seurat workflow ===================================================
# The dataset already contains the log2 normalized  counts in data layer
# Continuing the standard workflow after the normalization step
burdziak_seu <- FindVariableFeatures(burdziak_seu, # Variable genes
                                     selection.method = "vst", 
                                     nfeatures = 8000,
                                     verbose = FALSE) %>%
  ScaleData(.,verbose = FALSE) # centering and scaling (z-scores)

### Customized function for DEA ================================================
# Customized function based on FindMarkers that orders the results based on log2FC
# and p adjusted value. The results are written in a txt file.

custom_FindMarkers <- function(object, ident.1 ,ident.2 = NULL, min.pct = 0.01, 
                               logfc.threshold = 0.1,out.dir, file.name = "markers.txt") {
  
  # using Seurat's FindMarkers() and scCustomize's Add_Pct_Diff()
  cat("Performing DEA with FindMarkers.\n")
  markers <- FindMarkers(object = object, 
                         ident.1 = ident.1,         
                         ident.2 = ident.2,    
                         min.pct = min.pct,      
                         logfc.threshold = logfc.threshold,
                         test.use = "wilcox") %>%
    Add_Pct_Diff() %>% # add %diff between each group 
    arrange(desc(avg_log2FC), desc(p_val_adj)) # ordering by log2FC and pval
  
  # Write results
  cat("Writing results\n")
  file_out <- file.path(out.dir, file.name)
  write.table(markers, file_out , sep = "\t", quote = FALSE, col.names = TRUE)
  cat("Results available in file", file.name, "\n")
  return(markers)
}

### Differential Expression Analysis (DEA) =====================================
# outputs dir
dea_dir <- file.path(wd, "dea_outputs")
if (!dir.exists(dea_dir)) {dir.create(dea_dir)}

# Subset of clusters of interest: acinar, ADM and PDAC1
Idents(burdziak_seu) <- factor(burdziak_seu$cluster)
burdziak_subset <- subset(burdziak_seu, 
                          subset = cluster %in% c("2","7","3","5")) 
burdziak_subset

## Markers of each cluster of interest vs the rest
# Acinar vs all (ADM & PDAC)
acinar_all_markers <- custom_FindMarkers(object = burdziak_subset,
                                         ident.1 = 2,         
                                         out.dir = dea_dir, 
                                         file.name = "acinar_all_markers.txt")

# ADM vs Acinar & PDAC
adm_all_markers <- custom_FindMarkers(object = burdziak_subset,
                                      ident.1 = c(3,7),         
                                      out.dir = dea_dir, 
                                      file.name = "ADM_all_markers.txt")

early_adm_all_markers <- custom_FindMarkers(object = burdziak_subset,
                                            ident.1 = 7,         
                                            out.dir = dea_dir, 
                                            file.name = "earlyADM_all_markers.txt")

late_adm_all_markers <- custom_FindMarkers(object = burdziak_subset,
                                           ident.1 = 3,
                                           out.dir = dea_dir,
                                           file.name = "lateADM_all_markers.txt")
# PDAC vs Acinar & ADM
pdac_all_markers <- custom_FindMarkers(object = burdziak_subset,
                                       ident.1 = 5,         
                                       out.dir = dea_dir, 
                                       file.name = "earlyPDAC_all_markers.txt")


### GO enrichments =============================================================
# list of DEA markers results
markers_list <- list("Acinar" = acinar_all_markers,
                     "ADM" = adm_all_markers,
                     "early_ADM" = early_adm_all_markers,
                     "late_ADM" = late_adm_all_markers,
                     "early_PDAC" = pdac_all_markers)

plot_list <- list()

# GO enrichment loop for each marker result
for (i in seq_along(markers_list)) {
  markers <- markers_list[[i]]          # extract markers dataframe
  list_name <- names(markers_list[i])   # extract name
  
  # thresholds
  log2FC_thr <- 0.5
  pval_thr <- 0.05
  
  # Get differentially expressed genes
  degs <- markers %>%  dplyr::select(avg_log2FC, p_val_adj)
  degs_list <- list(
    upregulated = rownames(degs)[degs$avg_log2FC > log2FC_thr &
                                   degs$p_val_adj < pval_thr],
    
    downregulated = rownames(degs)[degs$avg_log2FC < -log2FC_thr &
                                     degs$p_val_adj < pval_thr])
  
  #GO enrichment in all ontologies (BP, MF and CC)
  goALL <- compareCluster(
    gene = degs_list,
    OrgDb = "org.Mm.eg.db",
    keyType = "SYMBOL",
    ont = "ALL")
  
  # plots
  plotALL <- dotplot(goALL, showCategory = 10) +
    facet_grid(Cluster ~ ONTOLOGY, scales = "free") + # group by up/downreg and onts
    theme_minimal() +  
    theme(strip.text.x = element_text(size = 18, face = "bold"), # highlight onts
          strip.text.y = element_blank(),           # no labs in Y facet groups
          axis.title.x = element_blank(),           # no X lab
          axis.text.y = element_text(size = 14)) +  # Y labs size
    ggtitle(paste(list_name, "GO enrichment"))
  
  filename = paste0(list_name, "_enrich.png")
  ggsave(filename = filename, plot = plotALL, width = 12, height = 10, path = dea_dir)
  plot_list[[i]] <- plotALL
}
#ADM cells show alterations in nuclear envelope and stressors response related categories:
adm_plot <- plot_list[[2]]
adm_plot

