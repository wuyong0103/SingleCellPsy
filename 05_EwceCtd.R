library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(tidyverse)

# Find all counts layers
data <- readRDS("Integrate_RunUMAP.rds")
meta <- read.table("Integrate-Meta.tsv", header = TRUE, sep = "\t")
data@meta.data$lineage <- meta$lineage
data@meta.data$sub_lineage <- meta$sub_lineage
data <- subset(data, sub_lineage != "Unknown")
count_layers <- Layers(data[["RNA"]])
count_layers <- count_layers[grepl("^counts", count_layers)]

# output to hdf5 file
out_file <- "Integrate_counts.h5"

# write to HDF5
for (i in seq_along(count_layers)) {
    cat("Processing", count_layers[i], "...\n")
  
    # extract raw count dgmatrix
    m <- GetAssayData(data, assay = "RNA", layer = count_layers[i])
    
    if (i == 1) {
        writeHDF5Array(m, filepath = out_file, name = "counts")
    } else {
        writeHDF5Array(m, filepath = out_file, name = "counts", append = TRUE)
    }
}

merged_counts <- HDF5Array(out_file, "counts")
dim(merged_counts)
class(merged_counts)

meta <- read_tsv("Integrate-Meta.tsv")
meta <- meta[meta$sub_lineage != "Unknown", ]
m <- match(colnames(merged_counts), meta$cell)
meta <- meta[m, ]

exp <- EWCE::sct_normalize(merged_counts)
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = exp, input_species = "human", output_species = "human", level2annot = meta$sub_lineage, no_cores=40)
annotLevels <- list(level1class=meta$lineage,level2class=meta$sub_lineage)
ctd <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = "Cortex_All", savePath=getwd(), no_cores=40)

for (age in c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder")) {
    meta_age <- meta[meta$age == age,]
    m <- match(meta_age$cell, colnames(exp))
    exp_age <- exp[,m]
    exp_age_dropped <- EWCE::drop_uninformative_genes(exp = exp_age, input_species = "human", output_species = "human", level2annot = meta_age$sub_lineage, no_cores=40)
    annotLevels_age <- list(level1class=meta_age$lineage, level2class=meta_age$sub_lineage)
    gn <- paste("Cortex_", age, sep = "")
    ctd_age <- EWCE::generate_celltype_data(exp = exp_age_dropped, annotLevels = annotLevels_age, groupName = gn, savePath=getwd(), no_cores=40)
}
