library(Seurat)
library(Matrix)
library(tidyverse)

samples <- c("SRR17380403", "SRR17380404", "SRR17380405", "SRR17380406")
paths <- file.path(samples, "outs", "filtered_feature_bc_matrix")

seurat_list <- list()

for (i in seq_along(samples)) {
  sample <- samples[i]
  path <- paths[i]
  counts <- Read10X(data.dir = path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)  
  seurat_obj$sample <- sample
  seurat_list[[sample]] <- seurat_obj
}

Caglayan2023 <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = samples,
  project = "Caglayan2023"
)

barcode <- tibble(barcode = Cells(Caglayan2023), donorid = Caglayan2023@meta.data$orig.ident)
meta <- read_tsv("SampleInfo.txt") %>% right_join(barcode, by = "donorid") %>% select(barcode, donorid, sex, age) %>% mutate(sex = ifelse(sex=="male", "M", "F"))
meta <- meta[match(Cells(Caglayan2023), meta$barcode), ]
Caglayan2023@meta.data <- Caglayan2023@meta.data[,1:3]
saveRDS(Caglayan2023, file = "Caglayan2023.rds")
write.csv(meta, "meta_Caglayan2023.csv")
