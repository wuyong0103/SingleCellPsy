library(Matrix)
library(Seurat)
exp <- readRDS("count_matrices.rds")
ctrl_samples <- c("Biopsy", "c_p1", "c_p2", "c_p3", "c_p4", "GTS217_2", "GTS219", "GTS233", "CTR215", "CTR240")

seurat_list <- lapply(ctrl_samples, function(sample_name) {
  mat <- exp[[sample_name]]
  colnames(mat) <- paste0(sample_name, "-", colnames(mat))
  rownames(mat) <- make.unique(rownames(mat))
  CreateSeuratObject(counts = mat)
})

Pfisterer2020 <- merge(seurat_list[[1]], y = seurat_list[-1], project = "Pfisterer2020")
saveRDS(Pfisterer2020, file = "Pfisterer2020.rds")

barcode <- tibble(barcode = Cells(Pfisterer2020)) %>% mutate(Name = sub("(^.*?)-.*", "\\1", barcode))
meta <- read_csv("sample_info.csv") %>% right_join(barcode, by = "Name") %>% 
    select(barcode, Name, Sex, Age) %>% 
    rename(donorid = Name, sex = Sex, age = Age) %>% 
    mutate(sex = ifelse(sex=="male", "M", "F"), age = paste0(age, "y"))
write.csv(meta, "meta_Pfisterer2020.csv")
