library(Matrix)
library(Seurat)
Steyn2024 <- readRDS("GSE280569_Steyn_etal_seurat_object.rds")
dim(Steyn2024)
#3000 144438
raw_counts <- GetAssayData(Steyn2024, assay = "RNA", layer = "counts")
dim(raw_counts)
#32037 144438
unique(Steyn2024@meta.data$biological_replicate)
#12 individuals

meta <- as_tibble(Steyn2024@meta.data) %>% select(...1, biological_replicate, sex) %>% 
    rename(barcode=...1, donorid=biological_replicate) %>% 
    mutate(donorid=sub("_yr_old", "y", donorid), age=sub("_yr_old.*", "y", donorid))
write.csv(meta, "meta_Steyn2024.csv")

Steyn2024 <- CreateSeuratObject(counts = raw_counts, project = "Steyn2024")
saveRDS(Steyn2024, file = "Steyn2024.rds")
