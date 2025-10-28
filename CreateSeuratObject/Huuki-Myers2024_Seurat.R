library(Seurat)
library(tidyverse)
h5fs <- list.files(".", pattern = "matrix\\.h5$")

seurat_list <- list()

for (sample in h5fs) {
    data <- Read10X_h5(sample)
    name <- sub("_filtered_feature_bc_matrix.h5", "", sample)
    colnames(data) <- paste0(name, "-", colnames(data))
    seurat_obj <- CreateSeuratObject(counts = data, project = name)
    seurat_list[[name]] <- seurat_obj
}

Huuki-Myers2024 <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    project = "Huuki-Myers2024"
)

dim(Huuki-Myers2024)

barcode <- tibble(barcode = Cells(seurat_merged), file = seurat_merged@meta.data$orig.ident)
meta <- read_tsv("sampleinfo.txt") %>% right_join(barcode, by = "file") %>% select(barcode, donorid, sex, age)
meta <- meta[match(Cells(seurat_merged), meta$barcode), ]
dim(meta)

saveRDS(seurat_merged, file = "Huuki-Myers2024.rds")
write.csv(meta, "meta_Huuki-Myers2024.csv")
