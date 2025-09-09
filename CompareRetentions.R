#Run citeseq.R before
setwd("/Users/natemaretzki/Desktop/Obradovic Lab/standard-workflow-main")
library(viper)
library(SingleR)
library(celldex)
library(umap)
library(ggplot2)
library(Seurat)
library(matrixStats)
library(dplyr)
library(cluster)
library(pheatmap)
library(Hmisc)
library(dplyr)
library(scales)
library(zellkonverter)
library(DropletUtils)

source("classes/utils.r")
source("classes/loader.r")
source("classes/SVpreprocessor.r")
source("classes/metacell_generator.r")
source("classes/runnercustom.r")
source("classes/clusterer.r")
source("classes/plotter.r")

SingleR_custom <- function(test, trained_labels, ref_mat, metric = "pearson") {
  common_genes <- sort(intersect(rownames(test), rownames(ref_mat)))
  test <- test[common_genes, , drop=FALSE]
  ref_mat <- ref_mat[common_genes, , drop=FALSE]
  stopifnot(identical(rownames(test), rownames(ref_mat)))
  
  # Correlation similarity
  similarity <- cor(as.matrix(ref_mat), as.matrix(test), method = "pearson")
  
  # Assign best label
  best_indices <- apply(similarity, 2, which.max)
  best_labels <- trained_labels[best_indices]
  
  DataFrame(
    scores = I(t(similarity)),
    labels = best_labels
  )
}

SingleR_custom2 <- function(test, trained_labels, ref_mat, metric = "pearson") {
  common_genes <- sort(intersect(rownames(test), rownames(ref_mat)))
  test <- test[common_genes, , drop = FALSE]
  ref_mat <- ref_mat[common_genes, , drop = FALSE]
  stopifnot(identical(rownames(test), rownames(ref_mat)))
  
  # Drop zero-variance rows
  test_sd <- matrixStats::rowSds(as.matrix(test))
  ref_sd  <- matrixStats::rowSds(as.matrix(ref_mat))
  keep <- (test_sd > 1e-6) & (ref_sd > 1e-6)
  test <- test[keep, , drop = FALSE]
  ref_mat <- ref_mat[keep, , drop = FALSE]
  
  # Pearson correlation
  similarity <- cor(as.matrix(ref_mat), as.matrix(test), method = metric)
  
  # Ensure 2D matrix
  if (is.null(dim(similarity))) {
    similarity <- matrix(similarity, nrow = 1)
    colnames(similarity) <- colnames(test)
    rownames(similarity) <- colnames(ref_mat)
  }
  
  best_indices <- apply(similarity, 2, function(x) which.max(x)[1])
  best_labels <- as.character(trained_labels)[best_indices]
  
  DataFrame(
    scores = I(t(similarity)),
    labels = best_labels
  )
}

recode_label <- function(x, mapping) {
  for (canon in names(mapping)) {
    if (x %in% mapping[[canon]]) return(canon)
  }
  return("T Cell/Myeloid Doublet?")  # Default fallback
}

calculate_purity <- function(predicted, ground_truth) {
  tab <- table(predicted, ground_truth)
  sum(apply(tab, 1, max)) / length(predicted)
}

calculate_overall_accuracy <- function(predicted, ground_truth) {
  mean(predicted == ground_truth)
}

per_label_accuracy <- function(predicted, ground_truth) {
  df <- data.frame(predicted = predicted, ground_truth = ground_truth)
  df %>%
    dplyr::group_by(ground_truth) %>%
    dplyr::summarise(
      accuracy = mean(predicted == ground_truth),
      .groups = "drop"
    )
}

set_fixed_depth <- function(expr_mat, target_depth) {
  cell_depths <- Matrix::colSums(expr_mat)
  keep_cells <- which(cell_depths >= target_depth)
  expr_mat <- expr_mat[, keep_cells, drop = FALSE]
  
  out <- matrix(0, nrow = nrow(expr_mat), ncol = ncol(expr_mat))
  rownames(out) <- rownames(expr_mat)
  colnames(out) <- colnames(expr_mat)
  
  for (i in seq_len(ncol(expr_mat))) {
    counts <- expr_mat[, i]
    total <- sum(counts)
    
    if (total <= target_depth) {
      out[, i] <- counts  # leave unchanged
    } else {
      sampled_genes <- sample(
        rep(seq_along(counts), counts),
        size = target_depth
      )
      out[, i] <- tabulate(sampled_genes, nbins = nrow(expr_mat))
    }
  }
  
  if (inherits(expr_mat, "dgCMatrix")) {
    out <- Matrix(out, sparse = TRUE)
  }
  
  return(out)
}

label_map <- list(
  "B Cells" = c("B-cells"),
  "CD4+ T-cells" = c("CD4+ T-cells"),
  "CD8+ T-cells" = c("CD8+ T-cells"),
  "HSC" = c("HSC"),
  "Myeloid" = c("Monocytes", "Neutrophils"),
  "NK Cells" = c("NK cells"),
  "T Cell/Myeloid Doublet?" = c("Eosinophils", "Erythrocytes", "Endothelial cells", "Epithelial cells", "Skeletal muscle")
)

BP <- BlueprintEncodeData(rm.NA = "rows")

allowed_labels <- unlist(label_map)
keep_ref_idx <- which(BP$label.main %in% allowed_labels)
labels_subset <- as.character(BP$label.main[keep_ref_idx])

true_labels <- integrated_seurat$mapped_idents

#--- PATHS AND SUCH ---

query_path <- "/Users/natemaretzki/Desktop/Obradovic Lab/SingleViper_Data/query"
ref_path <- "/Users/natemaretzki/Desktop/Obradovic Lab/SingleViper_Data/ref"
viper_output <- "/Users/natemaretzki/Desktop/Obradovic Lab/viper_query"

regulator_dir_path <- paste0("/Users/natemaretzki/Desktop/Obradovic Lab/regulators/human_hugo")
regulator_files <- list(
  cotfs = file.path(regulator_dir_path, "cotfs-hugo.txt"),
  surface = file.path(regulator_dir_path, "surface-hugo.txt"),
  sig = file.path(regulator_dir_path, "sig-hugo.txt"),
  tfs = file.path(regulator_dir_path, "tfs-hugo.txt")
)
aracne_exec <- "/Users/natemaretzki/Desktop/ARACNe3/build/src/app/ARACNe3_app"

query_aracne_outdir <- "/Users/natemaretzki/Desktop/Obradovic Lab/output_query"

depths <- c(8000, 6000, 4000, 2000, 1000, 500)

purity_SV_list <- list()
purity_SR_list <- list()
acc_SV_list <- list()
acc_SR_list <- list()
overall_acc_SV <- list()
overall_acc_SR <- list()

for(depth in depths) {

#--- LOAD ---

files <- list.files(query_aracne_outdir, full.names = TRUE)
unlink(files, recursive = TRUE)

ref = readRDS("/Users/natemaretzki/Desktop/Obradovic Lab/SingleViper_Data/ref/ref.rds")
query = readRDS("/Users/natemaretzki/Desktop/Obradovic Lab/SingleViper_Data/query/query.rds")

ref_expr <- SummarizedExperiment::assay(ref, "logcounts")
query_expr <- as.matrix(query)

rownames(query_expr) <- gsub("^HUMAN_", "", rownames(query_expr))
rownames(query_expr) <- gsub("^ERCC_", "", rownames(query_expr))

query_expr <- set_fixed_depth(query_expr, depth)

ref_cpm <- sweep(ref_expr, 2, colSums(ref_expr), FUN = "/") * 1e6
query_cpm <- sweep(query_expr, 2, colSums(query_expr), FUN = "/") * 1e6

query_seurat_for_aracne <- CreateSeuratObject(counts = query_cpm)
query_seurat_for_aracne <- SetAssayData(query_seurat_for_aracne, slot = "data", new.data = query_cpm)

overlap_genes <- sort(intersect(rownames(ref_cpm), rownames(query_cpm)))
ref_cpm <- ref_cpm[overlap_genes, ]
query_cpm <- query_cpm[overlap_genes, ]

ref_seurat <- CreateSeuratObject(counts = ref_cpm)
VariableFeatures(ref_seurat) <- NULL
ref_seurat <- FindVariableFeatures(ref_seurat, selection.method = "vst", nfeatures = 3000)
hvg_genes <- VariableFeatures(ref_seurat)
top_genes <- hvg_genes

ref_cpm <- ref_cpm[top_genes, ]
query_cpm <- query_cpm[top_genes, ]

ref_log <- log1p(ref_cpm)
query_log <- log1p(query_cpm)

ref_log <- as.matrix(ref_log)
query_log <- as.matrix(query_log)

ref_means <- rowMeans(ref_log)
ref_sds <- apply(ref_log, 1, sd)
ref_sds[ref_sds == 0] <- 1e-6

ref_z <- sweep(ref_log, 1, ref_means, "-")
ref_z <- sweep(ref_z, 1, ref_sds, "/")

query_z <- sweep(query_log, 1, ref_means, "-")
query_z <- sweep(query_z, 1, ref_sds, "/")

keep_genes <- (ref_sds > 0.2) & (abs(ref_means) < 5)
query_z <- query_z[keep_genes, ]

keep_genes <- matrixStats::rowSds(as.matrix(query_z)) > 1e-6
query_z <- query_z[keep_genes, , drop = FALSE]

#--- CLUSTER/METACELL ---

query_seurat_clust <- CreateSeuratObject(counts = query_expr)
query_seurat_clust <- NormalizeData(query_seurat_clust, verbose = FALSE)
query_seurat_clust <- FindVariableFeatures(query_seurat_clust, verbose = FALSE)
query_seurat_clust <- ScaleData(query_seurat_clust, verbose = FALSE)

clusterer_query <- Clusterer$new(query_seurat_clust, verbose = FALSE)
clusterer_query$run_clustering()

silhouette_results_query <- clusterer_query$calc_silhouette_scores()

best_resolution_query <- silhouette_results_query$best_resolution
query_seurat_clust <- clusterer_query$set_best_clusters(best_resolution_query)

top_genes <- clusterer_query$find_top_genes()

plotter <- Plotter$new(query_seurat_clust, query_aracne_outdir)
plotter$plot_silhouette_scores(silhouette_results_query$mean_scores,
                               silhouette_results_query$sd_scores)

plotter$plot_umap_clusters()

query_seurat_clust <- plotter$seurat_obj

query_metacell_gen <- MetacellGenerator$new(query_seurat_clust, query_aracne_outdir)
query_metacells <- query_metacell_gen$generate_metacell_matrices()

#--- ARACNE/VIPER ---

query_runner <- Runner$new(
  query_metacells,
  query_aracne_outdir,
  aracne_exec,
  query_aracne_outdir,
  regulator_files,
  query_seurat_for_aracne,
  viper_output,
  threads = 4,
  seed = 42
)

query_runner$run_aracne()
query_runner$create_rds()

regulon_files <- list.files(query_aracne_outdir, pattern = "_pruned\\.rds$", full.names = TRUE)
regulon_list <- lapply(regulon_files, readRDS)
all_regulons <- do.call(c, regulon_list)

all_regulons <- lapply(all_regulons, function(reg) {
  valid_targets <- intersect(names(reg$tfmode), rownames(ref_z))
  
  if (length(valid_targets) == 0) return(NULL)
  
  reg$tfmode <- reg$tfmode[valid_targets]
  
  # Check if likelihood exists and is the right length
  if (!is.null(reg$likelihood) && length(reg$likelihood) == length(reg$tfmode)) {
    reg$likelihood <- reg$likelihood[valid_targets]
  } else {
    # Fallback: infer from tfmode magnitude
    inferred_likelihood <- abs(reg$tfmode)
    inferred_likelihood <- inferred_likelihood / max(inferred_likelihood, na.rm = TRUE)
    reg$likelihood <- inferred_likelihood
  }
  
  return(reg)
})
all_regulons <- all_regulons[!sapply(all_regulons, is.null)]

ref_viper <- viper(ref_z, all_regulons)
query_viper <- viper(query_z, all_regulons)

#--- SingleViper ---
ref_mat_subset <- ref_viper[, keep_ref_idx]

common_genes <- intersect(rownames(query_viper), rownames(ref_mat_subset))
query_viper <- query_viper[common_genes, , drop = FALSE]
ref_mat_subset <- ref_mat_subset[common_genes, , drop = FALSE]

SingleViperPred <- SingleR_custom2(
  test = query_viper,
  trained_labels = labels_subset,
  ref_mat = ref_mat_subset,
  metric = "pearson"
)
SingleViperPred$labels_canon <- sapply(SingleViperPred$labels, recode_label, mapping = label_map)

#--- SingleR ---
ref_log_subset <- ref_log[, keep_ref_idx]

SingleRPred <- SingleR_custom2(
  test = query_log,
  trained_labels = labels_subset,
  ref_mat = ref_log_subset,
  metric = "pearson"
)
SingleRPred$labels_canon <- sapply(SingleRPred$labels, recode_label, mapping = label_map)

#--- CALCULATE PURITIES/ACCURACIES ---

SVLabels <- SingleViperPred$labels_canon
SRLabels <- SingleRPred$labels_canon
names(SVLabels) <- colnames(query_viper)
names(SRLabels) <- colnames(query_viper)

common_cells <- intersect(colnames(query_viper), colnames(integrated_seurat))
SVLabels_aligned <- SVLabels[common_cells]
SRLabels_aligned <- SRLabels[common_cells]
true_labels_aligned <- integrated_seurat$mapped_idents[common_cells]

purity_SV <- calculate_purity(SVLabels_aligned, true_labels_aligned)
acc_SV <- per_label_accuracy(SVLabels_aligned, true_labels_aligned)
overall_accuracy_SV <- calculate_overall_accuracy(SVLabels_aligned, true_labels_aligned)

purity_SR <- calculate_purity(SRLabels_aligned, true_labels_aligned)
acc_SR <- per_label_accuracy(SRLabels_aligned, true_labels_aligned)
overall_accuracy_SR <- calculate_overall_accuracy(SRLabels_aligned, true_labels_aligned)

purity_SV_list[[as.character(depth)]] <- purity_SV
purity_SR_list[[as.character(depth)]] <- purity_SV
acc_SV_list[[as.character(depth)]] <- acc_SV
acc_SR_list[[as.character(depth)]] <- acc_SR
overall_acc_SV[[as.character(depth)]] <- overall_accuracy_SV
overall_acc_SR[[as.character(depth)]] <- overall_accuracy_SR

cat(("finished with"), depth)
}

#--- PLOT ---

cat(paste("Finished!"))
k

valid_keys <- Reduce(intersect, list(
  names(purity_SV_list),
  names(purity_SR_list),
  names(acc_SV_list),
  names(acc_SR_list)
))

valid_retentions <- sort(as.numeric(valid_keys))
valid_keys <- as.character(valid_retentions)

purity_SV_vals <- unlist(purity_SV_list[valid_keys])
purity_SR_vals <- unlist(purity_SR_list[valid_keys])
acc_SV_vals <- sapply(valid_keys, function(k) acc_SV_list[[k]]["accuracy"])
acc_SR_vals <- sapply(valid_keys, function(k) acc_SR_list[[k]]["accuracy"])

plot(valid_retentions, purity_SV_vals, type = "b", col = "blue",
     xlab = "Retention", ylab = "Score", main = "Purity and Accuracy vs Retention",
     xlim = c(1, 0.2),
     ylim = c(0,1))
lines(unlist(valid_retentions), unlist(purity_SR_vals), type = "b", col = "red")
lines(unlist(valid_retentions), unlist(acc_SV_vals), type = "b", col = "green")
lines(unlist(valid_retentions), unlist(acc_SR_vals), type = "b", col = "purple")
legend("bottomleft",
       legend = c("Purity SV", "Purity SR", "Accuracy SV", "Accuracy SR"),
       col = c("blue", "red", "green", "purple"),
       lty = 1, pch = 1)

umap_coords <- umap(t(query_z))
umap_df <- as.data.frame(umap_coords$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$label <- SingleViperPred$labels_canon

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label)) +
  geom_point(size = 0.6) +
  theme_minimal() +
  labs(title = "VIPER-based UMAP with SingleR labels")


