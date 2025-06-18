## prep data for web app ##
library(Seurat)
library(data.table)
library(Matrix)
library(rhdf5)


## -----------------------------------------------------------------------------
## Inhibitory YANA DEVELOPMENT

load("/datastore_share/Users/ykotlyarenko/Inhibitory_cells_development/6_All_inhibitory_datasets/Datasets/Inhibitory_datasets.Rdata")

#inhibitory_datasets_cellID_vec <- colnames(Inhibitory_datasets)

## df for umap plotting:
umap_embedding <- as.data.frame(Inhibitory_datasets@reductions$umap2@cell.embeddings)
## invert axis:
umap_embedding$UMAP2_1 <- -umap_embedding$UMAP2_1
## add metadata (hard-coded)
umap_embedding$stage <- as.character(Inhibitory_datasets$Collection_stage[rownames(umap_embedding)])
umap_embedding$stage <- factor(umap_embedding$stage, levels = c("E12","E14","E16","P0"))
umap_embedding$experiment <- Inhibitory_datasets$Experiment2[rownames(umap_embedding)]
umap_embedding$experiment <- factor(umap_embedding$experiment, levels = c("WT","CFSE","LINEAGE"))
umap_embedding$study <- "Bright et al. 2025"
umap_embedding$cluster <- Inhibitory_datasets$Fine_annotation[rownames(umap_embedding)]
umap_embedding$cluster <- factor(umap_embedding$cluster, levels = c(
  "Fabp7","Top2a","Fabp7_Ccnd2","Ube2c","Nkx2_1","Abracl","Npy","Maf_Sst","Snhg11_Lhx8","Snhg11","Tcf4_Nr2f2",
  "Tshz1","Six3_Gucy1a3","Gucy1a3","Ebf1_Isl1"
))
umap_embedding$class <- Inhibitory_datasets$Main_splits[rownames(umap_embedding)]
umap_embedding$class[umap_embedding$class == "Interneurons"] <- "Interneuron Precursor"
umap_embedding$class[umap_embedding$class == "Projections"] <- "Inhibitory Projection Neuron Precursor"
umap_embedding$cellID <- rownames(umap_embedding)


write.table(umap_embedding, "data/inhibitory_datasets_umap2_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
## save expression values as hdf5:
log_count_mtx <- t(GetAssayData(Inhibitory_datasets, assay = "RNA", slot = "data"))


if(file.exists("data/inhibitory_datasets.h5")) {
  system("rm data/inhibitory_datasets.h5")
}

h5createFile("data/inhibitory_datasets.h5")
## save cell IDs and gene names:
h5createGroup("data/inhibitory_datasets.h5", "cell_IDs")
h5write(matrix(colnames(Inhibitory_datasets@assays$RNA@data), ncol = 1), "data/inhibitory_datasets.h5", "cell_IDs/cell_id_mtx")

h5createGroup("data/inhibitory_datasets.h5", "gene_names")
h5write(matrix(rownames(Inhibitory_datasets@assays$RNA@data), ncol = 1), "data/inhibitory_datasets.h5", "gene_names/gene_name_mtx")

## write log-normalized expression mtx:
h5createGroup("data/inhibitory_datasets.h5", "log_expression")

h5createDataset("data/inhibitory_datasets.h5", "log_expression/log_expr_mtx", c(nrow(log_count_mtx), ncol(log_count_mtx)), 
                storage.mode = "double", chunk = c(nrow(log_count_mtx), 100), level = 2)

h5write(matrix(log_count_mtx, ncol = ncol(log_count_mtx)), "data/inhibitory_datasets.h5", "log_expression/log_expr_mtx")
h5closeAll()




## -----------------------------------------------------------------------------
## INHIBITORY POSTNATAL

sticr_seurat <- readRDS("/data/mayerlab/neuhaus/CA_paper/STICR.basic.seuratobject.RDS")

## df for umap plotting:
umap_embedding <- as.data.frame(sticr_seurat@reductions$umap@cell.embeddings)

umap_embedding$UMAP2_1 <- umap_embedding$UMAP_1
umap_embedding$UMAP2_2 <- umap_embedding$UMAP_2

## add metadata (hard-coded)
umap_embedding$class <- sticr_seurat$refined_COUP_class
umap_embedding$class[umap_embedding$class == "Remove"] <- "Mixed"
umap_embedding$cluster <- sticr_seurat$refined_COUP_clust
umap_embedding$study <- "Bandler et al. 2022"
umap_embedding$dataset <- sticr_seurat$orig.ident
umap_embedding$cellID <- rownames(umap_embedding)
stage_vec <- c(
  "CA199" = "P8",
  "CA203" = "P7",
  "CA211" = "P9",
  "CA222" = "P7",
  "CA233" = "P5",
  "CA238" = "P7",
  "CA239" = "P7",
  "CA240" = "P7",
  "CA241" = "P7",
  "CA242" = "P7",
  "CA287" = "P13/P14",
  "CA289" = "P15",
  "CA304" = "P7",
  "LV1" = "P8",
  "LV3" = "P8",
  "LV4" = "P9",
  "LV7" = "P11",
  "LV12" = "P7"
)
umap_embedding$stage <- stage_vec[umap_embedding$dataset]
umap_embedding$stage <- factor(umap_embedding$stage, levels = c("P5","P7","P8","P9","P11","P13/P14","P15"))
exp_vec <- c(
  "CA199" = "STICR E12.5 - P8",
  "CA203" = "STICR E13.5 - P7",
  "CA211" = "STICR E12.5 - P9",
  "CA222" = "STICR E13.5 - P7",
  "CA233" = "STICR E13.5 - P5",
  "CA238" = "STICR E13.5 - P7",
  "CA239" = "STICR E14.5 - P7",
  "CA240" = "STICR E14.5 - P7",
  "CA241" = "STICR E13.5 - P7",
  "CA242" = "STICR E14.5 - P7",
  "CA287" = "STICR E13.5 - P13/P14",
  "CA289" = "STICR E12.5 - P15",
  "CA304" = "STICR E14.5 - P7",
  "LV1" = "STICR E10.5 - P8",
  "LV3" = "STICR E10.5 - P8",
  "LV4" = "STICR E10.5 - P9",
  "LV7" = "STICR E10.5 - P11",
  "LV12" = "STICR E12.5 - P7"
)
umap_embedding$experiment <- exp_vec[umap_embedding$dataset]
umap_embedding$experiment <- factor(umap_embedding$experiment, levels = c(
  "STICR E10.5 - P8","STICR E10.5 - P9","STICR E10.5 - P11",
  "STICR E12.5 - P7","STICR E12.5 - P8","STICR E12.5 - P9","STICR E12.5 - P15",
  "STICR E13.5 - P5","STICR E13.5 - P7","STICR E13.5 - P13/P14",
  "STICR E14.5 - P7"
))

write.table(umap_embedding, "data/STICR_umap_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


## save expression values as hdf5:
log_count_mtx <- t(GetAssayData(sticr_seurat, assay = "RNA", slot = "data"))

if(file.exists("data/STICR_datasets.h5")) {
  system("rm data/STICR_datasets.h5")
}

h5createFile("data/STICR_datasets.h5")
## save cell IDs and gene names:
h5createGroup("data/STICR_datasets.h5", "cell_IDs")
h5write(matrix(colnames(sticr_seurat@assays$RNA@data), ncol = 1), "data/STICR_datasets.h5", "cell_IDs/cell_id_mtx")

h5createGroup("data/STICR_datasets.h5", "gene_names")
h5write(matrix(rownames(sticr_seurat@assays$RNA@data), ncol = 1), "data/STICR_datasets.h5", "gene_names/gene_name_mtx")

## write log-normalized expression mtx:
h5createGroup("data/STICR_datasets.h5", "log_expression")

h5createDataset("data/STICR_datasets.h5", "log_expression/log_expr_mtx", c(nrow(log_count_mtx), ncol(log_count_mtx)), 
                storage.mode = "double", chunk = c(nrow(log_count_mtx), 100), level = 2)

h5write(matrix(log_count_mtx, ncol = ncol(log_count_mtx)), "data/STICR_datasets.h5", "log_expression/log_expr_mtx")
h5closeAll()




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
db_map <- c("dorsal_Arlotta" = "Di Bella et al. 2021", "ventral_CA_Rachel" = "Bandler et al. 2022", "ventral_Yana" = "Bright et al. 2025")
umap_embedding$study <- db_map[umap_embedding$study]
umap_embedding$study <- factor(umap_embedding$study, levels = c("Bright et al. 2025","Bandler et al. 2022","Di Bella et al. 2021"))

umap_embedding$cluster <- EI_seurat$Gene_Annotation
# umap_embedding$cluster <- factor(umap_embedding$cluster, levels = c(
#   "Gas1_Ldha","Hes1_Fabp7","Fabp7_Mt3","Hist1h1b_Top2a","Ccnd2_Nudt4","Nkx2-1_Lhx8","Npy_Nxph1","Sst_Maf","Nr2f2_Nr2f1","Isl1_Zfp503","Foxp1_Gucy1a3","Ebf1_Foxp1",
#   "Neurog2_Rrm2","Neurog2_Eomes","Neurod2_Neurod6","Neurod6_Mef2c"
# ))

class_vec <- c("Hes1_Fabp7" = "Mitotic", "Neurod2_Neurod6" = "Excitatory Neuron Precursor", "Neurog2_Rrm2" = "Mitotic", 
               "Gas1_Ldha" = "Mitotic", "Fabp7_Mt3" = "Mitotic", "Neurod6_Mef2c" = "Excitatory Neuron Precursor", 
               "Neurog2_Eomes" = "Mitotic", "Ccnd2_Nudt4" = "Inhibitory Neuron Precursor", "Npy_Nxph1" = "Inhibitory Neuron Precursor", 
               "Nr2f2_Nr2f1" = "Inhibitory Neuron Precursor", "Hist1h1b_Top2a" = "Mitotic", "Ebf1_Foxp1" = "Inhibitory Neuron Precursor", 
               "Nkx2-1_Lhx8" = "Inhibitory Neuron Precursor", "Isl1_Zfp503" = "Inhibitory Neuron Precursor", "Sst_Maf" = "Inhibitory Neuron Precursor",
               "Foxp1_Gucy1a3" = "Inhibitory Neuron Precursor")
umap_embedding$class <- class_vec[umap_embedding$cluster]
#umap_embedding$class <- factor(umap_embedding$class, levels = c("Mitotic", "Inhibitory Neuron Precursor", "Excitatory Neuron Precursor"))
umap_embedding$cellID <- rownames(umap_embedding)


write.table(umap_embedding, "data/EI_merged_umap2_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


## save expression values as hdf5:
log_count_mtx <- t(GetAssayData(EI_seurat, assay = "RNA", slot = "data"))

if(file.exists("data/EI_merged_datasets.h5")) {
  system("rm data/EI_merged_datasets.h5")
}

h5createFile("data/EI_merged_datasets.h5")
## save cell IDs and gene names:
h5createGroup("data/EI_merged_datasets.h5", "cell_IDs")
h5write(matrix(colnames(EI_seurat@assays$RNA@data), ncol = 1), "data/EI_merged_datasets.h5", "cell_IDs/cell_id_mtx")

h5createGroup("data/EI_merged_datasets.h5", "gene_names")
h5write(matrix(rownames(EI_seurat@assays$RNA@data), ncol = 1), "data/EI_merged_datasets.h5", "gene_names/gene_name_mtx")

## write log-normalized expression mtx:
h5createGroup("data/EI_merged_datasets.h5", "log_expression")

h5createDataset("data/EI_merged_datasets.h5", "log_expression/log_expr_mtx", c(nrow(log_count_mtx), ncol(log_count_mtx)), 
                storage.mode = "double", chunk = c(nrow(log_count_mtx), 100), level = 2)

h5write(matrix(log_count_mtx, ncol = ncol(log_count_mtx)), "data/EI_merged_datasets.h5", "log_expression/log_expr_mtx")
h5closeAll()
