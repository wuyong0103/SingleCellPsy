library(Seurat)
library(Matrix)
library(dplyr)

sample_paths <- list.dirs("/home/lilab/wuyong/project/scRNA/data/Schirmer2019Nature", recursive = FALSE)
samples <- basename(sample_paths)
paths <- file.path(samples, "outs", "filtered_feature_bc_matrix")

seurat_list <- list()

for (i in seq_along(samples)) {
  sample <- samples[i]
  path <- paths[i]
  counts <- Read10X(data.dir = path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  seurat_list[[sample]] <- seurat_obj
}

Schirmer2019 <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = samples,
  project = "Schirmer2019"
)

saveRDS(Schirmer2019, file = "Schirmer2019.rds")
dim(Schirmer2019)

meta <- Schirmer2019@meta.data
meta$barcode <- rownames(meta)
meta <- meta %>% rename(id=orig.ident) %>% 
    left_join(sampleinfo, by="id") %>% 
    select(barcode, sample, sex, age) %>% 
    mutate(age=paste0(age, "y"), sex=ifelse(sex=="male", "M", "F")) %>% 
    rename(donorid=sample)
write.csv(meta, "meta_Schirmer2019.csv")
