library(Matrix)
library(Seurat)
library(tidyverse)
data_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
genes <- read_tsv("genes.txt")
seurat_list <- list()
for (sample_dir in data_dirs) {
    counts <- readMM(paste0(sample_dir, "/counts_fil.mtx"))
    barcode <- read_tsv(paste0(sample_dir, "/col_metadata.tsv"))
    colnames(counts) <- barcode$Barcode
    rownames(counts) <- make.unique(genes$name)
    sample_name <- sub(".*_PN_([0-9]+)_snRNA-[A-Z0-9]+", "PN\\1", sample_dir)
    seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
    seurat_list[[sample_name]] <- seurat_obj
}

barcode <- tibble(barcode = Cells(Pineda2024), donorid = Pineda2024@meta.data$orig.ident)
meta <- read_tsv("sampleinfo.txt") %>% right_join(barcode, by = "donorid") %>% select(barcode, donorid, sex, age)
meta <- meta[match(Cells(Pineda2024), meta$barcode), ]
dim(Pineda2024)
dim(meta)
saveRDS(Pineda2024, file = "Pineda2024.rds")
write.csv(meta, "meta_Pineda2024.csv")
