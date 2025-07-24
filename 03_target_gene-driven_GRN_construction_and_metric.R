library(Seurat)
library(scRank)
library(dplyr)
library(harmony)
library(tidyr)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(tidygraph)
library(enrichR)
library(ggplot2)
library(forcats)
library(knitr)
library(readr)
library(ggtext)
library(tidyverse)
library(Matrix)


gex_metadata_f <- "./rna_metadata.tsv"
gex_metadata <- read.table(gex_metadata_f, header = TRUE, sep = "\t")
gex_metadata$batch <- gsub("_", "-", gex_metadata$batch)

#scRNA-seq Seurat object

folder <- './in/'
flist <- list.files(path = folder)

seurat_list <- list()

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

for (file_path in flist) {
  name <- gsub("[[:blank:]]", "", paste(folder, file_path))
  counts <- Read10X_h5(name)

  obj <- CreateSeuratObject(counts = counts$"Gene Expression", assay = "RNA", min.cells = 10, min.features = 0, project = gsub("_raw_featureflist_bc_matrix.h5", "", file_path))
  obj@meta.data$sample <- gsub("-raw-feature-bc-matrix.h5", "", gsub("_", "-", file_path))
  obj@meta.data$cells <- rownames(obj@meta.data)
  obj@meta.data$cell_samp <- paste(obj@meta.data$cells, obj@meta.data$sample, sep = "_")
  obj$mitoRatio <- PercentageFeatureSet(object = obj, pattern = "^MT-")
  obj$mitoRatio <- obj@meta.data$mitoRatio / 100
  obj <- subset(x = obj, subset= (mitoRatio < 0.10))
  obj@meta.data$cell_samp <- paste(obj@meta.data$cells, obj@meta.data$sample, sep = "_")
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")  
  obj <- NormalizeData(obj)
  genesPresentInSeurat <- rownames(obj@assays$RNA)
  s.genes <- s.genes[s.genes %in% genesPresentInSeurat]
  g2m.genes <- g2m.genes[g2m.genes %in% genesPresentInSeurat]
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  seurat_list[[gsub("_raw_feature_bc_matrix.h5", "", file_path)]] <- obj
}

rna <- merge(seurat_list[[1]], y = seurat_list[-1])

new_cell_names <- rna@meta.data$cell_samp
names(new_cell_names) <- colnames(rna)
rna <- RenameCells(object = rna, new.names = new_cell_names)

dict <- setNames(gex_metadata$cell_type_global, gex_metadata$cell_samp)
rna@meta.data <- rna@meta.data %>% mutate(cell_type_global = dict[cell_samp])
dict <- setNames(gex_metadata$ASDAS_response_CRP, gex_metadata$cell_samp)
rna@meta.data <- rna@meta.data %>% mutate(ASDAS_response_CRP = dict[cell_samp])
dict <- setNames(gex_metadata$time_point, gex_metadata$cell_samp)
rna@meta.data <- rna@meta.data %>% mutate(time_point = dict[cell_samp])

rna <- subset(rna, subset = cell_samp %in% gex_metadata$cell_samp)
rna[["RNA"]] <- JoinLayers(rna[["RNA"]])
rna <- subset(rna, subset = time_point %in% c('1'))

rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 3000)
rna <- ScaleData(rna, vars.to.regress = "G2M.Score")
rna <- RunPCA(rna, pc.genes = seurat_combined@var.genes, npcs = 50, verbose = FALSE)
rna <- RunHarmony(rna, "sample")
rna <- RunUMAP(object = rna, dims = 1:50)

new_cell_names <- rna@meta.data$cell_samp
names(new_cell_names) <- colnames(rna)
rna <- RenameCells(object = rna, new.names = new_cell_names)

dict <- setNames(gex_metadata$cell_type_global, gex_metadata$cell_samp)
rna@meta.data <- rna@meta.data %>% mutate(cell_type_global = dict[cell_samp])
dict <- setNames(gex_metadata$ASDAS_response_CRP, gex_metadata$cell_samp)
rna@meta.data <- rna@meta.data %>% mutate(ASDAS_response_CRP = dict[cell_samp])
dict <- setNames(gex_metadata$time_point, gex_metadata$cell_samp)
rna@meta.data <- rna@meta.data %>% mutate(time_point = dict[cell_samp])

rna <- subset(rna, subset = cell_samp %in% gex_metadata$cell_samp)
rna[["RNA"]] <- JoinLayers(rna[["RNA"]])
rna <- subset(rna, subset = time_point %in% c('1'))

rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 3000)
rna <- ScaleData(rna, vars.to.regress = "G2M.Score")
rna <- RunPCA(rna, pc.genes = seurat_combined@var.genes, npcs = 50, verbose = FALSE)
rna <- RunHarmony(rna, "sample")
rna <- RunUMAP(object = rna, dims = 1:50)

#saveRDS(rna, "./rna_seurat.rds")


rna@meta.data$ct_asdas_tp <- paste0(rna@meta.data$cell_type_global, "_", rna@meta.data$ASDAS_response_CRP, "_", rna@meta.data$time_point)
rna@meta.data$ct_asdas_tp <- gsub("_1_", "_12_", rna@meta.data$ct_asdas_tp)
rna@meta.data$ct_asdas_tp <- gsub("_2_", "_12_", rna@meta.data$ct_asdas_tp)
Idents(rna) <- rna@meta.data$ct_asdas_tp

#Run scRank to obtain GRNs

for (target in c('TNF', 'JAK1')){
    print(target)
    obj <- CreateScRank(input = rna,
                     species = 'human',
                     cell_type = 'ct_asdas_tp',
                     target = target)

    obj <- scRank::Constr_net(obj)

    obj <- scRank::rank_celltype(obj)
    saveRDS(obj, paste0("./scRank/scRank_ct_asdas_0vs12_tp1_", target, ".rds"))
}


#GRN metric

library(igraph)
library(Matrix)
library(tidyverse)

targets <- c("TNF", "JAK1")
metrics_to_compute <- c("degree_weighted", "eigenvector")

for (target in targets) {
  message("🔍 Target: ", target)
  scres <- readRDS(paste0("./scRank/scRank_ct_asdas_0vs12_tp1_", target, ".rds"))

  all_results <- list()

  for (metric in metrics_to_compute) {
    for (grn_name in names(scres@net)) {

      grn <- scres@net[[grn_name]]
      grn_sparse <- as(grn, "dgCMatrix")

      edge_df <- as.data.frame(summary(grn_sparse)) %>%
        rename(from = i, to = j, weight = x)

      rownames_map <- rownames(grn)
      edge_df$from <- rownames_map[edge_df$from]
      edge_df$to <- rownames_map[edge_df$to]

      edge_df <- edge_df %>%
        filter(from != to) %>%
        rowwise() %>%
        mutate(pair = paste(sort(c(from, to)), collapse = "_")) %>%
        distinct(pair, .keep_all = TRUE) %>%
        ungroup() %>%
        select(from, to, weight)

      edf <- edge_df %>% mutate(weight = abs(weight))
      if (nrow(edf) == 0) next

      g <- graph_from_data_frame(edf, directed = FALSE)

      metric_values <- switch(
        metric,
        degree_weighted = strength(g, weights = E(g)$weight),
        eigenvector = eigen_centrality(g, weights = E(g)$weight)$vector
      )

      result <- data.frame(
        GRN = grn_name,
        gene = names(metric_values),
        metric = metric,
        metric_value = metric_values
      )

      all_results[[length(all_results) + 1]] <- result
    }
  }


  combined_results <- bind_rows(all_results)
  outdir <- file.path("./results/scRank_metrics_results")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write_csv(combined_results, file.path(outdir, paste0("combined_", target, "_metrics.csv")))

}

