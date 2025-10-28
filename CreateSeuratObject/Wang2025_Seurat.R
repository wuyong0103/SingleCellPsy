library(Seurat)
library(tidyverse)
Wang2025 <- readRDS("snMultiome_atlas_Seurat_object.rds")
meta <- Wang2025@meta.data
meta$barcode <- rownames(meta)
meta <- as_tibble(meta) %>% select(barcode, Ident, sex, Estimated_postconceptional_age_in_days) %>% 
    rename(donorid=Ident, age=Estimated_postconceptional_age_in_days) %>%
    mutate(age = ifelse(age < 280, paste0(ceiling(age/7), "pcw"), paste0(ceiling((age - 280) / 365), "y"))) %>%
    mutate(donorid = sub("^([A-Za-z]+-[M0-9]+).*$", '\\1', donorid))
write.csv(meta, "meta_Wang2025.csv")
system("ln -s snMultiome_atlas_Seurat_object.rds Wang2025.rds")
