library(Matrix)
library(tidyverse)
library(Seurat)
counts <- Read10X_h5("GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
meta <- read_csv("GSE174367_snRNA-seq_cell_meta.csv.gz") %>% filter(Diagnosis == "Control") %>% 
    select(Barcode, SampleID, Sex, Age) %>% 
    rename(barcode=Barcode, donorid=SampleID, sex=Sex, age=Age) %>% 
    mutate(age=paste0(age, "y"))
cells <- intersect(colnames(counts), meta$barcode)
counts <- counts[, colnames(counts) %in% cells]
Morabito2021 <- CreateSeuratObject(counts = counts, , project = "Morabito2021")
saveRDS(Morabito2021, file = "Morabito2021.rds")
write.csv(meta, "meta_Morabito2021.csv")
