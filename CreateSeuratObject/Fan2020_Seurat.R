library(tidyverse)
library(Seurat)
expr <- read.table("GSE120046_brain_all_UMIcounts.txt", row.names = 1, header = T)
meta <- read_tsv("meta.txt") %>% filter(tissue != "PONS")
common_cells <- intersect(colnames(expr), meta$cell)
expr_sub <- expr[, common_cells]
meta_sub <- meta[common_cells, ]
Fan2020 <- CreateSeuratObject(counts = expr_sub, meta.data = meta_sub, project = "Fan2020")
p <- VlnPlot(Fan2020, features = "XIST", group.by = "sample") + ggtitle("XIST Expression by Sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("XIST_expression_by_sample.png", plot = p, width = 8, height = 6, dpi = 300)
y <- VlnPlot(Fan2020, features = c("RPS4Y1", "KDM5D", "DDX3Y", "UTY", "ZFY", "EIF1AY", "TSPY1"), group.by = "sample") + ggtitle("Ygenes Expression by Sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Ygenes_expression_by_sample.png", plot = y, width = 8, height = 6, dpi = 300)

barcode <- tibble(barcode = Cells(Fan2020), donorid = sub("_.*", "", Cells(Fan2020)), age = sub("^HE([0-9]+)W_.*", "\\1pcw", Cells(Fan2020)))
meta <- read_tsv("sample.txt") %>% right_join(barcode, by = "donorid") %>% select(barcode, donorid, sex, age)
meta <- meta[match(colnames(Fan2020), meta$barcode), ]
Fan2020@meta.data <- Fan2020@meta.data[,1:3]
saveRDS(Fan2020, file = "Fan2020.rds")
write.csv(meta, "meta_Fan2020.csv")
