library(Matrix)
library(Seurat)
library(tidyverse)
counts <- Matrix::readMM("GSE213982_combined_counts_matrix.mtx.gz")
gene <- read.table("GSE213982_combined_counts_matrix_genes_rows.csv.gz", header = T)
cell <- read.table("GSE213982_combined_counts_matrix_cells_columns.csv.gz", header = T)
rownames(counts) <- gene$x
colnames(counts) <- cell$x
meta <- read.table("meta_Maitra2023.txt", header = T)
counts <- counts[, colnames(counts) %in% meta$barcode]
Maitra2023 <- CreateSeuratObject(counts = counts, project = "Maitra2023")
dim(Maitra2023)
saveRDS(Maitra2023, file = "Maitra2023.rds")

meta <- tibble(barcode = rownames(Maitra2023@meta.data), donorid=sub("^([FM][0-9]+).*", "\\1", rownames(Maitra2023@meta.data)))
meta <- read_tsv("sampleinfo.txt") %>% 
    right_join(meta, by = "donorid") %>% 
    select(barcode, donorid, sex, age) %>% 
    mutate(age = paste0(age, "y"))
write.csv(meta, "meta_Maitra2023.csv")
