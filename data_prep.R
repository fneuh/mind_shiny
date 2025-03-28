## prep data for web app ##
library(Seurat)
library(data.table)
library(Matrix)


## -----------------------------------------------------------------------------
## Inhibitory YANA

load("/datastore_share/Users/ykotlyarenko/Inhibitory_cells_development/6_All_inhibitory_datasets/Datasets/Inhibitory_datasets.Rdata")

#inhibitory_datasets_cellID_vec <- colnames(Inhibitory_datasets)

## df for umap plotting:
umap_embedding <- as.data.frame(Inhibitory_datasets@reductions$umap2@cell.embeddings)
## switch axis:
umap_embedding$UMAP2_1 <- -umap_embedding$UMAP2_1
## add metadata (hard-coded)
umap_embedding$stage <- as.character(Inhibitory_datasets$Collection_stage)
umap_embedding$stage <- factor(umap_embedding$stage, levels = c("E12","E14","E16","P0"))
umap_embedding$experiment <- Inhibitory_datasets$Experiment2
umap_embedding$experiment <- factor(umap_embedding$experiment, levels = c("WT","CFSE","LINEAGE"))
umap_embedding$study <- "Kotylarenko et al. 2024"
umap_embedding$cluster <- Inhibitory_datasets$Fine_annotation
umap_embedding$cluster <- factor(umap_embedding$cluster, levels = c(
  "Fabp7","Top2a","Fabp7_Ccnd2","Ube2c","Nkx2_1","Abracl","Npy","Maf_Sst","Snhg11_Lhx8","Snhg11","Tcf4_Nr2f2",
  "Tshz1","Six3_Gucy1a3","Gucy1a3","Ebf1_Isl1"
))
umap_embedding$cellID <- rownames(umap_embedding)

# 
#umap_plot_ggplot(umap_embedding, col_attr = "cluster", split_attr = "experiment", point_size = 0.2)
write.table(umap_embedding, "data/inhibitory_datasets_umap2_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


## df for feature plot:
log_count_mtx <- as.data.table(t(GetAssayData(Inhibitory_datasets, assay = "RNA", slot = "data")))

variable_genes <- VariableFeatures(FindVariableFeatures(Inhibitory_datasets, nfeatures = 5000))
mean_expr <- apply(log_count_mtx, 2, mean)
th <- 0.5
sum(mean_expr > th)
expr_high_genes <- colnames(log_count_mtx)[mean_expr > th]
genes_sub <- unique(c(variable_genes, expr_high_genes))



fwrite(log_count_mtx, file = "data/inhibitory_datasets_log_count_mtx.csv", sep = ",")
#writeMM(log_count_mtx, file = "data/inhibitory_datasets_log_count_mtx.mtx")
#system("gzip data/inhibitory_datasets_log_count_mtx.mtx")

log_count_HVG_mtx <- log_count_mtx[, colnames(log_count_mtx) %in% VariableFeatures(Inhibitory_datasets), with = FALSE]
fwrite(log_count_HVG_mtx, file = "data/inhibitory_datasets_log_count_HVG_mtx.csv", sep = ",")


log_count_sub_mtx <- log_count_mtx[, colnames(log_count_mtx) %in% genes_sub, with = FALSE]
fwrite(log_count_sub_mtx, file = "data/inhibitory_datasets_log_count_sub_mtx.csv", sep = ",")

scaled_mtx <- as.data.table(t(GetAssayData(Inhibitory_datasets, assay = "RNA", slot = "scale.data")))
fwrite(scaled_mtx, file = "data/inhibitory_datasets_scaled_count_mtx.csv", sep = ",")



## -----------------------------------------------------------------------------
## MERGED INHIBITORY


## -----------------------------------------------------------------------------
## MERGED EI

EI_seurat <- readRDS("/datastore_share/Users/neuhaus/dorsal_ventral_comp/Paper/results/final_results/1_01_data_prep_dorsal_ventral_WT/EXCIT_INHIBIT_cleaned_sub.rds")

## df for umap plotting:
umap_embedding <- as.data.frame(EI_seurat@reductions$umap2@cell.embeddings)

## add metadata (hard-coded)
umap_embedding$stage <- as.character(EI_seurat$Stage_DV1)
umap_embedding$stage <- factor(umap_embedding$stage, levels = c("E12","E13","E14","E15","E16"))
umap_embedding$experiment <- "WT"
umap_embedding$study <- EI_seurat$Dataset_big
db_map <- c("dorsal_Arlotta" = "Di Bella et al. 2021", "ventral_CA_Rachel" = "Bandler et al. 2022", "ventral_Yana" = "Kotylarenko et al. 2024")
umap_embedding$study <- db_map[umap_embedding$study]
umap_embedding$study <- factor(umap_embedding$study, levels = c("Kotylarenko et al. 2024","Bandler et al. 2022","Di Bella et al. 2021"))

umap_embedding$cluster <- EI_seurat$Gene_Annotation
umap_embedding$cluster <- factor(umap_embedding$cluster, levels = c(
  "Gas1_Ldha","Hes1_Fabp7","Fabp7_Mt3","Hist1h1b_Top2a","Ccnd2_Nudt4","Nkx2-1_Lhx8","Npy_Nxph1","Sst_Maf","Nr2f2_Nr2f1","Isl1_Zfp503","Foxp1_Gucy1a3","Ebf1_Foxp1",
  "Neurog2_Rrm2","Neurog2_Eomes","Neurod2_Neurod6","Neurod6_Mef2c"
))
umap_embedding$cellID <- rownames(umap_embedding)

# 
#umap_plot_ggplot(umap_embedding, col_attr = "cluster", split_attr = "study", point_size = 0.2)
write.table(umap_embedding, "data/EI_merged_umap2_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

## df for feature plot:

log_count_mtx <- as.data.table(t(GetAssayData(EI_seurat, assay = "RNA", slot = "data")))

variable_genes <- VariableFeatures(FindVariableFeatures(EI_seurat, nfeatures = 5000))
mean_expr <- apply(log_count_mtx, 2, mean)
th <- 0.5
sum(mean_expr > th)
expr_high_genes <- colnames(log_count_mtx)[mean_expr > th]
genes_sub <- unique(c(variable_genes, expr_high_genes))

fwrite(log_count_mtx, file = "data/EI_merged_log_count_mtx.csv", sep = ",")
#writeMM(log_count_mtx, file = "data/EI_merged_log_count_mtx.mtx")
#system("gzip data/inhibitory_datasets_log_count_mtx.mtx")

log_count_HVG_mtx <- log_count_mtx[, colnames(log_count_mtx) %in% VariableFeatures(EI_seurat), with = FALSE]
fwrite(log_count_HVG_mtx, file = "data/EI_merged_log_count_HVG_mtx.csv", sep = ",")

log_count_sub_mtx <- log_count_mtx[, colnames(log_count_mtx) %in% genes_sub, with = FALSE]
fwrite(log_count_sub_mtx, file = "data/EI_merged_log_count_sub_mtx.csv", sep = ",")

scaled_mtx <- as.data.table(t(GetAssayData(EI_seurat, assay = "RNA", slot = "scale.data")))
fwrite(scaled_mtx, file = "data/EI_merged_scaled_count_mtx.csv", sep = ",")
