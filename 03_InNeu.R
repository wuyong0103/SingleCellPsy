library(Seurat)
library(Matrix)
library(harmony)
options(future.globals.maxSize = 10 * 1024^3)
data <- readRDS("Integrate_RunUMAP.rds")

Idents(data) <- "CellType"


inneu_cells <- colnames(data)[data$CellType %in% c("In-PVALB", "In-Progenitor", "In-SST", "In-VIP", "In-RELN", "In-LAMP5")]
cat(paste0("[INFO]: Total In Neuron cell number, ", length(inneu_cells), "\n"))

layer_names <- Layers(data[["RNA"]])
count_layer <- grep("count", layer_names, value = T)

seurat_list <- list()

for (layer in count_layer) {
    temp_count <- GetAssayData(data, assay = "RNA", layer = layer)
    common_cells <- intersect(colnames(temp_count), inneu_cells)
    cat(paste0("[INFO]:", length(common_cells), "\n"))
    if (length(common_cells) > 1){
        temp_count <- temp_count[, common_cells]
        seurat_list[[layer]] <- CreateSeuratObject(counts = temp_count, meta.data = data@meta.data[common_cells, 1:12])
        cat(paste0("[INFO]: ", layer, " ", dim(seurat_list[[layer]]), "\n"))
    }
}

InNeu <- merge(seurat_list[[1]], y = seurat_list[-1])
cat(paste0("[INFO]: Seurat object created successfully, ", dim(InNeu), "\n"))

InNeu <- NormalizeData(InNeu)
InNeu <- FindVariableFeatures(InNeu, selection.method='mean.var.plot')
InNeu <- ScaleData(InNeu, features=VariableFeatures(InNeu))
InNeu <- RunPCA(InNeu, features = VariableFeatures(InNeu), verbose=F)


png(filename = "ElbowPlot-InNeu.png", width = 20, height = 10, units = "in", res = 300)
ElbowPlot(object = InNeu, ndims = 40)
dev.off()

saveRDS(InNeu, file = "InNeu-RunPCA.rds")

InNeu <- readRDS("InNeu-RunPCA.rds")

N = 26
InNeu <- RunHarmony(InNeu, group.by.vars = "orig.ident", theta = 2, max.iter.harmony = 20, dims.use = 1:N)
InNeu <- FindNeighbors(InNeu, reduction = "harmony", dims = 1:N)
InNeu <- FindClusters(InNeu, resolution = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4))
InNeu <- RunUMAP(InNeu, reduction = "harmony", dims = 1:N)
#saveRDS(InNeu, "InNeu-RunUMAP.rds")


#====================DimPlot of different resolutions=======================================================
Idents(InNeu) <- "RNA_snn_res.0.02"
png(filename = "DimPlot_InNeu_0.02.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.04"
png(filename = "DimPlot_InNeu_0.04.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.06"
png(filename = "DimPlot_InNeu_0.06.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.08"
png(filename = "DimPlot_InNeu_0.08.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.1"
png(filename = "DimPlot_InNeu_0.1.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.2"
png(filename = "DimPlot_InNeu_0.2.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.3"
png(filename = "DimPlot_InNeu_0.3.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

Idents(InNeu) <- "RNA_snn_res.0.4"
png(filename = "DimPlot_InNeu_0.4.png", width = 40, height = 40, units = "in", res = 300)
DimPlot(object = InNeu, reduction = "umap")
dev.off()

#====================FeaturePlot of different markers=======================================================
png(filename = "FeaturePlot_Chandelier.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("THSD7A", "CA8", "ANK1", "GULP1", "FAM19A4", "NPNT", "PVALB", "CRH", "PLCXD3", "GPR149", "UNC5B", "PLEKHH2"), ncol = 4)
dev.off()

png(filename = "FeaturePlot_LAMP5-LHX6.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("LAMP5", "LHX6", "FGF13", "EYA4", "PDGFD", "NTNG1", "SV2C", "CHST9", "MYO16", "CA1", "SPHKAP", "KIT"), ncol = 4)
dev.off()

png(filename = "FeaturePlot_LAMP5.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("LAMP5", "KIT", "EYA4", "SV2C", "TRPC3", "FREM1", "EGFR", "PKP2", "CPLX3", "CXCL14", "BMP6", "TPD52L1"), ncol = 4)
dev.off()

png(filename = "FeaturePlot_SNCG.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("SNCG", "CXCL14", "THSD7B", "COL21A1", "NR2F2", "PROX1", "CHRNA7", "NDST4", "ADAM33", "GNG12", "CNR1", "NPAS1", "ADARB2"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_VIP.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("VIP", "TAC3", "SHISA8", "PROX1", "CALB2", "SLC22A3", "GPR149", "CHRNA2", "SYNPR", "LAMA3", "SLC24A3", "IQGAP2", "ADARB2"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_PAX6.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("RELN", "CNR1", "CXCL14", "ADARB2", "CHRNA7", "SLC35F4", "CRH", "WIF1", "SORCS3", "PLS3", "AP1S2", "DDR2", "RBMS3"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_PVALB.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("BTBD11", "ADAMTS17", "MYO5B", "SULF1", "TAC1", "RNF144B", "FGF12", "ANK1", "CNTNAP3", "DOCK11", "SOX6", "SLIT2", "ERBB4"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_SST.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("SST", "FLT3", "GRIK1", "PLCH1", "NXPH1", "STXBP6", "GRIK3", "SPOCK3", "KIF26B", "XKR4", "SOX6", "SYNPR", "PAWR"), ncol = 3)
dev.off()

png(filename = "FeaturePlot_SST-CHODL.png", width = 40, height = 30, units = "in", res = 300)
FeaturePlot(InNeu, features = c("SST", "CHODL", "NPY", "CRHBP", "TRPC6", "NOS1", "CDHR1", "CCDC109B", "TAC1", "NXPH2", "SOX6", "KLF5", "NPY2R"), ncol = 3)
dev.off()


Idents(In) <- "RNA_snn_res.0.2"
png(filename = "Dimplot_0.2_label.png", width = 20, height = 20, units = "in", res = 300, type = "cairo")
DimPlot(In, reduction = "umap", label = T)
dev.off()

In@meta.data$CellType <- In@meta.data$RNA_snn_res.0.2
In@meta.data$CellType <- as.character(In$CellType)
In@meta.data$CellType[In@meta.data$CellType %in% c("0", "11")] <- "In-PVALB"
In@meta.data$CellType[In@meta.data$CellType %in% c("1", "4")] <- "In-SST"
In@meta.data$CellType[In@meta.data$CellType == "10"] <- "In-Chandelier"
In@meta.data$CellType[In@meta.data$CellType == "8"] <- "In-LAMP5-LHX6"
In@meta.data$CellType[In@meta.data$CellType == "5"] <- "In-LAMP5"
In@meta.data$CellType[In@meta.data$CellType == "7"] <- "In-SNCG"
In@meta.data$CellType[In@meta.data$CellType == "9"] <- "In-PAX6"
In@meta.data$CellType[In@meta.data$CellType %in% c("6", "2")] <- "In-VIP"
In@meta.data$CellType[In@meta.data$CellType == "13"] <- "In-SST-CHODL"
In@meta.data$CellType[In@meta.data$CellType %in% c("3", "12")] <- "In-Progenitor"
In@meta.data$CellType[In@meta.data$CellType == "14"] <- "Unknown"

Idents(In) <- "CellType"
imeta <- In@meta.data
pdf("Dimplot_In0.2_Lable.pdf", width = 13, height = 10)
DimPlot(In, reduction = "umap", label = T)
dev.off()
saveRDS(InNeu, "InNeu-RunUMAP.rds")