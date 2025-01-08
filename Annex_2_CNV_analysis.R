################################################################################
##### ANNEX 2: CNV analysis via inferCNV (Burdziak et all. dataset)        #####
################################################################################

wd <- r"(C:\Users\javim\OneDrive - UNIVERSIDAD DE SEVILLA\TFM_UPO)"
setwd(wd)

### Libraries
library(tidyverse)
library(glue)
library(Seurat)  
library(infercnv)
library(patchwork)

# Load processed data (see Annex 1)
burdziak_seu <- LoadSeuratRds("./data/burdziak_seu_corrected.rds")

###  Clusters of interest ======================================================
# Setting clusters as identifiers
Idents(burdziak_seu) <-  factor(burdziak_seu$cluster)

# Checking the distribution of clusters among the different tissue states
# plots dir
plot_dir <- file.path(wd, "cls_plots")
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}

# Barplot: tissue states
colors <- c("chartreuse3", "cadetblue3", "gold", "chocolate2", "brown1", "deeppink4")
tissue_cls_barplot <- ggplot(burdziak_seu@meta.data, aes(x = factor(cluster), fill = tissue)) +
  geom_bar(position = "fill") + 
  labs( x = "Clusters",y = "Relative percentage", fill = "Tissue states") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors) + # personalized color
  theme_minimal(base_size = 16) # theme & text size
tissue_cls_barplot
ggsave(plot = tissue_cls_barplot, 
       filename = file.path(plot_dir, "barplot_cls_states.png"), 
       height = 4, width = 6)

# Based on cell annotations from reference paper:
# Acinar cells (cls 2) come from healthy pancreas tissue
# ADM cells from cls 7 can be considered as an early ADM,and come mainly from 
# regenerative tissue (N2 conditions)
# ADM cells from cls 3 can be considered as a later ADM state (CAE+KRAS),
# and come from pre-neoplasic and PanIN tissue
# PDAC cells from cls 5 can be considered as an earlier PDAC, coming only from 
# non-metastasic tumorigenic tissue

# Checking marker genes expression in clusters of interest
# Violin plot: markers
burdziak_subset <- subset(burdziak_seu, # Subset of clusters of interest
                          subset = cluster %in% c("2","7","3","5")) 

features <- c("Zg16", "Cpa1",       # PDAC Markers
              "Krt19", "Sox9",      # Ductal Markers
              "Kras")               #PDAC markers

violin_plot <- VlnPlot(burdziak_subset, features, stack = TRUE, flip = TRUE) +
  labs(x = "Clusters") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        strip.text.y = element_text(face = "bold", angle = 0))
violin_plot
ggsave(plot = violin_plot,
       filename = file.path(plot_dir, "violinplot_cls.png"), 
       height = 6, width = 6)

# Table summary
tissue_cls_summary <- burdziak_seu@meta.data %>%
  as_tibble() %>%
  count(cluster, tissue, condition) %>% 
  filter(cluster %in% c("2","7","3","5")) %>% 
  split(.$cluster)
tissue_cls_summary

# ==============================================================================

### CNV analysis
# Create dir for CNV analysis
cnv_path <- "./cnv_analysis/"
if (!dir.exists( cnv_path)) {dir.create( cnv_path)}
setwd(cnv_path)

## Inputs
# Raw counts
cts_subset <- LayerData(burdziak_subset, layer = "counts") # counts matrix

# Reference gene ordering file from Trinity Cancer Transcriptome Analysis Toolkit (CTAT)
url <- "https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions"
genome_pos_path <- "./GRCm39.vM32_gene_id_positions.txt"
download.file(url, destfile = genome_pos_path)

# Annotations file: gene name and annotation (cluster) cols
annot_path <- "./pdac_annotations.txt"

burdziak_subset@meta.data %>% mutate(
  cell_name = rownames(.),   # col with cell names
  annotation = case_when(    # col with annot based on cluster
    cluster == "2" ~ "Acinar cells",
    cluster == "7" ~ "Early ADM",
    cluster == "3" ~ "Late ADM",
    cluster == "5" ~ "PDAC")) %>%
  dplyr::select(cell_name, annotation) %>% 
  write_tsv(annot_path, col_names = FALSE) # txt file

# Creating InferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = cts_subset,      # raw counts
  annotations_file = annot_path,       # path to conditions annotations 
  gene_order_file = genome_pos_path,   # path to genome positions file
  delim = "\t",
  ref_group_names = "Acinar cells")   # reference group (Acinar cells)

# InferCNV run
infercnv::run(infercnv_obj, 
              cutoff = 0.1, # recommended for 10X Genomics
              output_format = "png", 
              out_dir = "./infercnv_run")