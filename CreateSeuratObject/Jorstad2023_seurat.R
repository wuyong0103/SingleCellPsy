library(Seurat)
library(tidyverse)

expr <- read.csv("matrix.csv", row.names = 1, check.names = FALSE)
expr <- t(expr)

library(Matrix)
sparse_expr <- Matrix(as.matrix(expr), sparse = TRUE)

Jorstad2023_SmartSeq <- CreateSeuratObject(counts = sparse_expr, project = "Jorstad2023_SmartSeq")

saveRDS(Jorstad2023, file = "Jorstad2023_SmartSeq.rds")

rds_files <- list.files(".", pattern = "^[ADMSV].*\\.RDS$", full.names = TRUE)
seurat_list <- lapply(seq_along(rds_files), function(i) {
        mat <- readRDS(rds_files[i])
        CreateSeuratObject(counts = mat)
})
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], project = "Jorstad2023")
annot <- read_tsv("/home/lilab/wuyong/project/scRNA/data/Nano2025NatNeu/Ensemble_Symbol.txt", col_names = c("ensembl", "name"))
genes_to_keep <- intersect(rownames(merged_seurat), annot$name)

meta <- read_csv("meta_10x_2_11_22.csv") %>% select(sample_id, donor) %>% rename(barcode = sample_id, donorid = donor)
cell_to_keep <- intersect(colnames(merged_seurat), meta$barcode)

seurat_filtered <- subset(merged_seurat, features = genes_to_keep, cells = cell_to_keep)
dim(seurat_filtered)
saveRDS(seurat_filtered, file = "Jorstad2023.rds")

meta <- meta %>% filter(barcode %in% cell_to_keep)
meta <- read_tsv("sample_info.txt") %>% right_join(meta, by = "donorid") %>% select(barcode, donorid, sex, age)
meta <- meta[match(rownames(seurat_filtered@meta.data), meta$barcode), ]
write.csv(meta, "meta_Jorstad2023.csv")
