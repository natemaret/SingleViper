Preprocessor <- R6::R6Class(
  "Preprocessor",
  public = list(
    seurat_obj = NULL,
    mt_threshold = 25,
    min_rna = 1000,
    max_rna = 15000,
    verbose = FALSE,
    blueprint_encode = NULL,
    
    initialize = function(seurat_obj,
                          mt_threshold = 25,
                          min_rna = 1000,
                          max_rna = 15000,
                          verbose = FALSE) {
      self$seurat_obj <- seurat_obj
      self$mt_threshold <- mt_threshold
      self$min_rna <- min_rna
      self$max_rna <- max_rna
      self$verbose <- verbose
      self$blueprint_encode <- celldex::BlueprintEncodeData()
    },
    
    preprocess_data = function() {
      message("Preprocessing Seurat object...")
      tryCatch({
        self$seurat_obj <- private$preprocess_seurat(self$seurat_obj)
        message("Preprocessing completed.")
        return(self$seurat_obj)
      }, error = function(e) {
        stop("Error during preprocessing: ", e$message)
      })
    }
  ),
  
  private = list(
    preprocess_seurat = function(seurat_object) {
      seurat_object <- private$calculate_percent_mt(seurat_object)
      seurat_object <- private$filter_cells(seurat_object)
      seurat_object <- private$normalize_data(seurat_object)
      seurat_object <- private$annotate_cells_with_singler(seurat_object)
      return(seurat_object)
    },
    
    calculate_percent_mt = function(seurat_object) {
      mitochondrial_genes <- grep("^MT-", rownames(seurat_object), value = TRUE)
      seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, features = mitochondrial_genes)
      return(seurat_object)
    },
    
    filter_cells = function(seurat_object) {
      seurat_object <- subset(
        seurat_object,
        subset = percent.mt < self$mt_threshold &
          nCount_RNA > self$min_rna &
          nCount_RNA < self$max_rna
      )
      return(seurat_object)
    },
    
    normalize_data = function(seurat_object) {
      seurat_object <- SCTransform(
        seurat_object,
        vars.to.regress = c("nCount_RNA", "percent.mt"),
        return.only.var.genes = FALSE,
        verbose = self$verbose,
        conserve.memory = TRUE
      )
      return(seurat_object)
    },
    
    annotate_cells_with_singler = function(seurat_object) {
      singler_results <- SingleR(
        test = seurat_object[["SCT"]]@counts,
        ref = self$blueprint_encode,
        labels = self$blueprint_encode$label.main
      )
      seurat_object$blueprint_labels <- singler_results$labels
      seurat_object$blueprint_pvals <- singler_results$scores
      return(seurat_object)
    }
  )
)