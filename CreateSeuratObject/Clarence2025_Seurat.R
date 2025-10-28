library(Seurat)
library(tidyverse)
cla <- readRDS("a97cd3b7-8083-4d4b-b2d8-681078d7c94f.rds")
annot <- read.table("/home/lilab/wuyong/project/scRNA/data/Nano2025NatNeu/Ensemble_Symbol.txt", header =F)
colnames(annot) <- c("ensembl_id", "gene_symbol")
gene_map <- annot[!duplicated(annot$ensembl_id), ]
gene_rename <- setNames(gene_map$gene_symbol, gene_map$ensembl_id)
gene_rename <- gene_rename[rownames(cla)]
gene_rename[is.na(gene_rename)] <- rownames(cla)[is.na(gene_rename)]
gene_rename <- make.unique(gene_rename)

DefaultAssay(cla) <- "RNA"
assay_obj <- cla[["RNA"]]
counts_mat <- assay_obj@counts
rownames(counts_mat) <- gene_rename

data_mat <- assay_obj@data
rownames(data_mat) <- gene_rename

if (nrow(assay_obj@scale.data) > 0) {
  scale_mat <- assay_obj@scale.data
  rownames(scale_mat) <- gene_rename
  assay_obj@scale.data <- scale_mat
}

meta_feat <- assay_obj@meta.features
rownames(meta_feat) <- gene_rename

assay_obj@counts <- counts_mat
assay_obj@data <- data_mat
assay_obj@meta.features <- meta_feat

cla[["RNA"]] <- assay_obj

Clarence2025 <- subset(cla, tissue=="dorsolateral prefrontal cortex")

meta <- tibble(barcode = Cells(Clarence2025), donorid = Clarence2025@meta.data$donor_id, sex = Clarence2025@meta.data$sex, age = Clarence2025@meta.data$development_stage) %>% mutate(sex = ifelse(sex=="male", "M", "F"), age = ifelse(age == "infant stage", "0y", sub("-year-old stage", "y", age)))
Clarence2025@meta.data <- Clarence2025@meta.data[,21:23]
colnames(Clarence2025@meta.data) <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
saveRDS(Clarence2025, file = "Clarence2025.rds")
write.csv(meta, "meta_Clarence2025.csv")

