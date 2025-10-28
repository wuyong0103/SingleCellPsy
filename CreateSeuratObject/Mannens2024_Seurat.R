library(Seurat)
library(Matrix)
library(tidyverse)

samples <- c("10X346_2_ABCD_2", "10X370_1_ABCD_2", "10X402_1_ABCD_1", "10X406_7_ABCD_2")
path <- "/home/lilab/wuyong/project/scRNA/data/Mannens2024Nature"
seurat_list <- list()

for (sample in samples) {
  message("Processing sample: ", sample)

  feature_file <- paste0(path, "/", sample, "/filtered_feature_bc_matrix/features.tsv.gz")
  barcode_file <- paste0(path, "/", sample, "/filtered_feature_bc_matrix/barcodes.tsv.gz")
  matrix_file  <- paste0(path, "/", sample, "/filtered_feature_bc_matrix/matrix.mtx.gz")

  features <- read.delim(feature_file, header = FALSE, sep="\t", stringsAsFactors = FALSE)
  colnames(features) <- c("id", "name", "type", "chr", "start", "end")

  gene_rows <- which(features$type == "Gene Expression" & grepl("^ENSG", features$id))
  gene_features <- features[gene_rows, ]

  mat <- readMM(matrix_file)
  barcodes <- readLines(barcode_file)
  sample_name <- gsub(".*X(\\d+)_.*", "\\1", sample)
  barcodes <- paste0(sample_name, "-", barcodes)  # 避免 barcode 冲突

  gene_mat <- mat[gene_rows, ]
  rownames(gene_mat) <- make.unique(gene_features$name)
  colnames(gene_mat) <- barcodes

  seu <- CreateSeuratObject(counts = gene_mat, project = sample_name)
  seurat_list[[sample_name]] <- seu
}

Mannens2024 <- merge(seurat_list[[1]], y = seurat_list[-1], project = "Mannens2024", add.cell.ids = names(seurat_list))
dim(Mannens2024)

saveRDS(Mannens2024, file = "Mannens2024.rds")

meta <- tibble(barcode = Cells(Mannens2024), sample=as.numeric(Mannens2024@meta.data$orig.ident))
meta <- read_tsv("sampleinfo.txt") %>% 
    right_join(meta, by = "sample") %>% 
    select(barcode, sample, sex, age) %>% 
    rename(donorid = sample)
write.csv(meta, "meta_Mannens2024.csv")
