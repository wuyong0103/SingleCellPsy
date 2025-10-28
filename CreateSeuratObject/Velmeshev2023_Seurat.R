library(Seurat)
library(tidyverse)
expr <- Read10X("./", gene.column = 1)
Velmeshev2023 <- CreateSeuratObject(counts = expr, project = "Velmeshev2023")
saveRDS(Velmeshev2023, file = "Velmeshev2023.rds")

meta <- read_tsv("meta.tsv") %>% select(cell, individual, sex, `age(days)`) %>% 
    rename(barcode=cell, donorid=individual, age=`age(days)`) %>% 
    mutate(age=ifelse(age<280, paste0(ceiling(age / 7), "pcw"), paste0(ceiling((age-280)/365), "y")), sex=sub("[emal]+", "", sex))
write.csv(meta, "meta_Velmeshev2023.csv")
