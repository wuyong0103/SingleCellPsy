library(Seurat)
library(Matrix)
library(harmony)
options(future.globals.maxSize = 10 * 1024^3)
data <- readRDS("integrate_RunUMAP.rds")

Idents(data) <- "CellType"

exneu_cells <- colnames(data)[data$CellType %in% c("ExNeu", "Ex-Progenitor")]
cat(paste0("[INFO]: Total Ex Neuron cell number, ", length(exneu_cells), "\n"))

layer_names <- Layers(data[["RNA"]])
count_layer <- grep("count", layer_names, value = T)

seurat_list <- list()

for (layer in count_layer) {
    temp_count <- GetAssayData(data, assay = "RNA", layer = layer)
    common_cells <- intersect(colnames(temp_count), exneu_cells)
    cat(paste0("[INFO]:", length(common_cells), "\n"))
    if (length(common_cells) > 1){
        temp_count <- temp_count[, common_cells]
        seurat_list[[layer]] <- CreateSeuratObject(counts = temp_count, meta.data = data@meta.data[common_cells, 1:12])
        cat(paste0("[INFO]: ", layer, " ", dim(seurat_list[[layer]]), "\n"))
    }
}

ExNeu <- merge(seurat_list[[1]], y = seurat_list[-1])
cat(paste0("[INFO]: Seurat object created successfully, ", dim(ExNeu), "\n"))

ExNeu <- NormalizeData(ExNeu)
ExNeu <- FindVariableFeatures(ExNeu, selection.method='mean.var.plot')
ExNeu <- ScaleData(ExNeu, features=VariableFeatures(ExNeu))
ExNeu <- RunPCA(ExNeu, features = VariableFeatures(ExNeu), verbose=F)


png(filename = "ElbowPlot-ExNeu.png", width = 20, height = 10, units = "in", res = 300)
ElbowPlot(object = ExNeu, ndims = 40)
dev.off()

saveRDS(ExNeu, file = "ExNeu-RunPCA.rds")

ExNeu <- readRDS("ExNeu-RunPCA.rds")

N = 24
ExNeu <- RunHarmony(ExNeu, group.by.vars = "orig.ident", theta = 2, max.iter.harmony = 20, dims.use = 1:N)
ExNeu <- FindNeighbors(ExNeu, reduction = "harmony", dims = 1:N)
ExNeu <- FindClusters(ExNeu, resolution = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3))
ExNeu <- RunUMAP(ExNeu, reduction = "harmony", dims = 1:N)
#saveRDS(ExNeu, "ExNeu-RunUMAP.rds")


#====================DimPlot of different resolutions=======================================================
Idents(ExNeu) <- "RNA_snn_res.0.02"
png(filename = "DimPlot_ExNeu_0.02.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()

Idents(ExNeu) <- "RNA_snn_res.0.04"
png(filename = "DimPlot_ExNeu_0.04.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()

Idents(ExNeu) <- "RNA_snn_res.0.06"
png(filename = "DimPlot_ExNeu_0.06.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()

Idents(ExNeu) <- "RNA_snn_res.0.08"
png(filename = "DimPlot_ExNeu_0.08.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()

Idents(ExNeu) <- "RNA_snn_res.0.1"
png(filename = "DimPlot_ExNeu_0.1.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()

Idents(ExNeu) <- "RNA_snn_res.0.2"
png(filename = "DimPlot_ExNeu_0.2.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()

Idents(ExNeu) <- "RNA_snn_res.0.3"
png(filename = "DimPlot_ExNeu_0.3.png", width = 50, height = 40, units = "in", res = 300)
DimPlot(object = ExNeu, reduction = "umap")
dev.off()


#====================FeaturePlot of different markers=======================================================
png(filename = "FeaturePlot_ExET.png", width = 20, height = 22, units = "in", res = 300)
FeaturePlot(ExNeu, features = c("FEZF2", "BCL11B", "CTIP2", "TCF4", "EGR1", "RBP4", "SIM1", "ADCYAP1", "TCERG1L", "SEMA3E"), ncol = 4)
dev.off()

png(filename = "FeaturePlot_ExIT.png", width = 20, height = 18, units = "in", res = 300)
FeaturePlot(ExNeu, features = c("LHX2", "SATB2", "CUX1", "CUX2", "RORB", "TOX", "TBR2", "PLXND1", "TLX3"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_ExCT.png", width = 20, height = 12, units = "in", res = 300)
FeaturePlot(ExNeu, features = c("NTSR1", "TLE4", "TBR1", "FOXP2", "SAMD3"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_ExNP.png", width = 20, height = 18, units = "in", res = 300)
FeaturePlot(ExNeu, features = c("FOS", "BHLHE22", "LMO3", "LAMP5", "ADRA1A", "RPRM", "C1QL2"), ncol = 3)
dev.off()

Idents(Ex) <- "RNA_snn_res.0.2"
png(filename = "Dimplot_0.2_label.png", width = 20, height = 20, units = "in", res = 300, type = "cairo")
DimPlot(Ex, reduction = "umap", label = T)
dev.off()

Ex@meta.data$CellType <- Ex@meta.data$RNA_snn_res.0.2
Ex@meta.data$CellType <- as.character(Ex$CellType)
Ex@meta.data$CellType[Ex@meta.data$CellType %in% c("0", "7", "11")] <- "Ex-L23IT"
Ex@meta.data$CellType[Ex@meta.data$CellType == "5"] <- "Ex-L4IT"
Ex@meta.data$CellType[Ex@meta.data$CellType %in% c("1", "14")] <- "Ex-L5IT"
Ex@meta.data$CellType[Ex@meta.data$CellType == "6"] <- "Ex-L6IT"
Ex@meta.data$CellType[Ex@meta.data$CellType == "12"] <- "Ex-L6IT-Car3"
Ex@meta.data$CellType[Ex@meta.data$CellType == "15"] <- "Ex-L5ET"
Ex@meta.data$CellType[Ex@meta.data$CellType == "10"] <- "Ex-L6b"
Ex@meta.data$CellType[Ex@meta.data$CellType == "8"] <- "Ex-L6CT"
Ex@meta.data$CellType[Ex@meta.data$CellType == "13"] <- "Ex-L56NP"
Ex@meta.data$CellType[Ex@meta.data$CellType == "3"] <- "Ex-Progenitor"
Ex@meta.data$CellType[Ex@meta.data$CellType %in% c("2", "4", "9")] <- "Ex-Inter"
saveRDS(ExNeu, "ExNeu-RunUMAP.rds")
