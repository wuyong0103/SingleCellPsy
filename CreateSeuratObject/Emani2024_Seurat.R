library(Seurat)
library(Matrix)
library(tidyverse)
library(data.table)

samples <- c("RT00372N", "RT00373N", "RT00374N", "RT00375N", "RT00376N", "RT00377N", "RT00378N", "RT00379N", "RT00380N", "RT00381N", "RT00382N", "RT00383N", "RT00384N", "RT00385N", "RT00386N", "RT00387N", "RT00388N", "RT00389N", "RT00390N", "RT00391N")
paths <- paste0("~/project/scRNA/data/Emani2024Science/SRA/", samples, "/outs/filtered_feature_bc_matrix")

seurat_list <- list()

for (i in seq_along(samples)) {
  sample <- samples[i]
  path <- paths[i]
  
  matrix_file <- paste0(path, "/matrix.mtx.gz")
  feature_file <- paste0(path, "/features.tsv.gz")
  barcode_file <- paste0(path, "/barcodes.tsv.gz")

  features <- fread(feature_file, header = FALSE)
  rna_features <- features[V3 == "Gene Expression" & !duplicated(V2)]
  rna_indices <- which(features$V3 == "Gene Expression" & !duplicated(features$V2))

  barcodes <- fread(barcode_file, header = FALSE)$V1
  expr_matrix <- readMM(matrix_file)
  expr_matrix <- as(expr_matrix, "dgCMatrix")

  rna_matrix <- expr_matrix[rna_indices, ]
  rownames(rna_matrix) <- rna_features$V2
  colnames(rna_matrix) <- barcodes

  seurat_obj <- CreateSeuratObject(counts = rna_matrix, project = sample)
  seurat_obj$sample <- sample
  seurat_list[[sample]] <- seurat_obj
}

Emani2024 <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = samples,
  project = "Emani2024"
)

barcode <- tibble(barcode = Cells(Emani2024), donorid = Emani2024@meta.data$orig.ident)
meta <- read_tsv("SampleInfo.txt") %>% right_join(barcode, by = "donorid") %>% select(barcode, donorid, sex, age) %>% mutate(sex = ifelse(sex=="male", "M", "F"))
meta <- meta[match(Cells(Emani2024), meta$barcode), ]
Emani2024@meta.data <- Emani2024@meta.data[,1:3]
saveRDS(Emani2024, file = "Emani2024.rds")
write.csv(meta, "meta_Emani2024.csv")

