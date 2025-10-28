library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
sce <- load("SCE_DLPFC-n3_tran-etal.rda")
sce
counts_dense <- as.matrix(counts(sce.dlpfc.tran))
Tran2021 <- CreateSeuratObject(counts = counts_dense, project = "Tran2021")
dim(Tran2021)
meta <- tibble(barcode = Cells(Tran2021), donorid = Tran2021@meta.data$orig.ident)
meta <- read_tsv("sampleinfo.txt") %>% right_join(meta, by = "donorid") %>% select(barcode, donorid, sex, age)
dim(meta)
meta <- meta[match(Cells(Tran2021), meta$barcode), ]
write.csv(meta, "meta_Tran2021.csv")
saveRDS(Tran2021, file = "Tran2021.rds")
