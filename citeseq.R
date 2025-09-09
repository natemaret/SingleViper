setwd("/Users/natemaretzki/Desktop/Obradovic Lab/standard-workflow-main")

source("classes/citesequtils.r")
source("classes/converter.r")
source("classes/loader.r")
source("classes/preprocessor.r")
source("classes/integrator.r")
source("classes/citeseqclusterer.r")
source("classes/plotter.r")
source("classes/metacell_generator.r")
source("classes/runner.r")
library(matrixStats)
library(Seurat)
library(SingleR)
library(R6)
library(celldex)
library(cluster)
library(pheatmap)
library(Hmisc)
library(dplyr)
library(umap)
library(ggplot2)
library(scales)


plot_output_path <- "/Users/natemaretzki/Desktop/Obradovic Lab/CiteseqOutput"

rna_counts <- as.sparse(read.csv("/Users/natemaretzki/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi 5.csv.gz", 
                       sep = ",", header = TRUE, row.names = 1))
rna_counts <- as.matrix(rna_counts)

adt_counts <- as.sparse(read.csv("/Users/natemaretzki/Downloads/GSE100866_CBMC_8K_13AB_10X-ADT_umi 4.csv.gz", 
                       sep = ",", header = TRUE, row.names = 1))
adt_counts <- as.matrix(adt_counts)

cite_obj <- CreateSeuratObject(counts = rna_counts)
cite_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)

cite_obj <- NormalizeData(cite_obj, assay = "ADT", normalization.method = "CLR")
cite_obj <- ScaleData(cite_obj, assay = "ADT")
VariableFeatures(cite_obj[["ADT"]]) <- rownames(cite_obj[["ADT"]])




DefaultAssay(cite_obj) <- "ADT"


clusterer <- Clusterer$new(cite_obj, verbose = FALSE,dims = 1:10)
clusterer$run_clustering()

silhouette_results <- clusterer$calc_silhouette_scores()
best_resolution <- silhouette_results$best_resolution

integrated_seurat <- clusterer$set_best_clusters(best_resolution)

top_genes <- clusterer$find_top_genes()

plotter <- Plotter$new(integrated_seurat, plot_output_path)

plotter$plot_silhouette_scores(silhouette_results$mean_scores,
                               silhouette_results$sd_scores)

plotter$plot_umap_clusters()






adt_features <- rownames(integrated_seurat[["ADT"]])

DoHeatmap(
  object = integrated_seurat,
  features = adt_features,
  assay = "ADT"
) + scale_fill_gradientn(
  colors = c("navy", "white", "firebrick"),  # Blue → White → Red
  name = "Expression"
)+
theme(legend.position = "right")




table(Idents(integrated_seurat))




cluster_map <- list(
  "Cluster1" = "CD4+ T-cells",
  "Cluster2" = "Myeloid",
  "Cluster3" = "NK Cells",
  "Cluster4" = "HSC",
  "Cluster5" = "CD8+ T-cells",
  "Cluster6" = "B Cells",
  "Cluster7" = "Myeloid",
  "Cluster8" = "T Cell/Myeloid Doublet?"
)


# Get current cluster identities (as integers or factors)
clusters <- as.character(Idents(integrated_seurat))

# Convert to ClusterX format
cluster_ids <- paste0("Cluster", as.integer(clusters) + 1)

# Map cluster IDs to cell type labels
named_labels <- unlist(cluster_map[cluster_ids])

# Assign names to match cell names in Seurat
names(named_labels) <- colnames(integrated_seurat)

# Add to Seurat metadata
integrated_seurat$mapped_idents <- named_labels

# Make a table
print(table(integrated_seurat$mapped_idents))



# Set the identity class to your mapped labels
Idents(integrated_seurat) <- "mapped_idents"

# Plot the UMAP
DimPlot(integrated_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Colored by mapped_idents")



print(table(integrated_seurat$mapped_idents))


marker_map <- list(
  "CD3" = "T cell",
  "CD4" = "CD4+ T cell",
  "CD8" = "CD8+ T cell",
  "CD14" = "Myeloid",
  "CD19" = "B cell",
  "CD56" = "NK cell",
  "CD34" = "HSC",
  "CD11c" = "Myeloid",
  "CD16" = "Myeloid",
  "CD45RA" = "T Cell",
  "CCR7" = "T Cell",
  "CCR5" = "T Cell"
)
