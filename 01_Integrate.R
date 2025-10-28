library(Seurat)
library(harmony)
datadir <- "/gpfs/home/wuyong/project/scRNA/data"
seurat_list <- list()
setwd("/gpfs/home/wuyong/project/scRNA/integrate")

seurat_list[[1]] <- readRDS(paste0(datadir, "/Fan2020SciAdv/Fan2020.rds"))
#raw    8324
#filter 2247
seurat_list[[2]] <- readRDS(paste0(datadir, "/Velmeshev2023Science/Velmeshev2023.rds"))
#raw    709372
#filter 586360
seurat_list[[3]] <- readRDS(paste0(datadir, "/Bakken2021Nature/Bakken2021.rds"))
#raw    76533
#filter 76450
seurat_list[[4]] <- readRDS(paste0(datadir, "/Caglayan2023Nature/Caglayan2023.rds"))
#raw    62699
#filter 55889
seurat_list[[5]] <- readRDS(paste0(datadir, "/Jorstad2023Science/Jorstad2023.rds"))
#raw    1155822
#filter 1153673
seurat_list[[6]] <- readRDS(paste0(datadir, "/Emani2024Science/Emani2024.rds"))
#raw    96673
#filter 96673
seurat_list[[7]] <- readRDS(paste0(datadir, "/Morabito2021NatNeu/Morabito2021.rds"))
#raw    22796
#filter 22191
seurat_list[[8]] <- readRDS(paste0(datadir, "/Pfisterer2020NC/Pfisterer2020.rds"))
#raw    62671
#filter 60266
seurat_list[[9]] <- readRDS(paste0(datadir, "/Schirmer2019Nature/Schirmer2019.rds"))
#raw    21609 
#filter 16391
seurat_list[[10]] <- readRDS(paste0(datadir, "/Clarence2025NatGenet/Clarence2025.rds"))
#raw    33422
#filter 28297
seurat_list[[11]] <- readRDS(paste0(datadir, "/Zhu2023SciAdv/Zhu2023.rds"))
#raw    45549
#filter 45535
seurat_list[[12]] <- readRDS(paste0(datadir, "/Batiuk2022SciAdv/Batiuk2022.rds"))
#raw    209053
#filter 198573
seurat_list[[13]] <- readRDS(paste0(datadir, "/Steyn2024NatGenet/Steyn2024.rds"))
#raw    144438
#filter 138124
seurat_list[[14]] <- readRDS(paste0(datadir, "/Wang2025Nature/Wang2025.rds"))
seurat_list[[14]] <- UpdateSeuratObject(seurat_list[[14]])
seurat_list[[14]]@meta.data <- seurat_list[[14]]@meta.data[,2:4]
colnames(seurat_list[[14]]@meta.data) <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
#raw    232328
#filter 226840
seurat_list[[15]] <- readRDS(paste0(datadir, "/Maitra2023NatCommun/Maitra2023.rds"))
#raw    67438
#filter 52750
seurat_list[[16]] <- readRDS(paste0(datadir, "/Tran2021Neuron/Tran2021.rds"))
#raw    11202 
#filter 11185
seurat_list[[17]] <- readRDS(paste0(datadir, "/Pineda2024Cell/Pineda2024.rds"))
#raw    78142
#filter 76775
seurat_list[[18]] <- readRDS(paste0(datadir, "/Huuki-Myers2024Science/Huuki-Myers2024.rds"))
#raw    106086
#filter 87420
seurat_list[[19]] <- readRDS(paste0(datadir, "/Mannens2024Nature/Mannens2024.rds"))
#raw    43083
#filter 33870

seurat_list[[1]]@meta.data$orig.ident <- "Fan2020"
Idents(seurat_list[[1]]) <- seurat_list[[1]]@"meta.data"$orig.ident
Idents(seurat_list[[2]]) <- seurat_list[[2]]@"meta.data"$orig.ident
Idents(seurat_list[[3]]) <- seurat_list[[3]]@"meta.data"$orig.ident
seurat_list[[4]]@meta.data$orig.ident <- "Caglayan2023"
Idents(seurat_list[[4]]) <- seurat_list[[4]]@"meta.data"$orig.ident
seurat_list[[5]]@meta.data$orig.ident <- "Jorstad2023"
Idents(seurat_list[[5]]) <- seurat_list[[5]]@"meta.data"$orig.ident
seurat_list[[6]]@meta.data$orig.ident <- "Emani2024"
Idents(seurat_list[[6]]) <- seurat_list[[6]]@"meta.data"$orig.ident
Idents(seurat_list[[7]]) <- seurat_list[[7]]@"meta.data"$orig.ident
seurat_list[[8]]@meta.data$orig.ident <- "Pfisterer2020"
Idents(seurat_list[[8]]) <- seurat_list[[8]]@"meta.data"$orig.ident
seurat_list[[9]]@meta.data$orig.ident <- "Schirmer2019"
Idents(seurat_list[[9]]) <- seurat_list[[9]]@"meta.data"$orig.ident
seurat_list[[10]]@meta.data$orig.ident <- "Clarence2025"
Idents(seurat_list[[10]]) <- seurat_list[[10]]@"meta.data"$orig.ident
seurat_list[[11]]@meta.data$orig.ident <- "Zhu2023"
Idents(seurat_list[[11]]) <- seurat_list[[11]]@"meta.data"$orig.ident
seurat_list[[12]]@meta.data$orig.ident <- "Batiuk2022"
Idents(seurat_list[[12]]) <- seurat_list[[12]]@"meta.data"$orig.ident
seurat_list[[13]]@meta.data$orig.ident <- "Steyn2024"
Idents(seurat_list[[13]]) <- seurat_list[[13]]@"meta.data"$orig.ident
seurat_list[[14]]@meta.data$orig.ident <- "Wang2025"
Idents(seurat_list[[14]]) <- seurat_list[[14]]@"meta.data"$orig.ident
Idents(seurat_list[[15]]) <- seurat_list[[15]]@"meta.data"$orig.ident
seurat_list[[16]]@meta.data$orig.ident <- "Tran2021"
Idents(seurat_list[[16]]) <- seurat_list[[16]]@"meta.data"$orig.ident
seurat_list[[17]]@meta.data$orig.ident <- "Pineda2024"
Idents(seurat_list[[17]]) <- seurat_list[[17]]@"meta.data"$orig.ident
seurat_list[[18]]@meta.data$orig.ident <- "Huuki-Myers2024"
Idents(seurat_list[[18]]) <- seurat_list[[18]]@"meta.data"$orig.ident
seurat_list[[19]]@meta.data$orig.ident <- "Mannens2024"
Idents(seurat_list[[19]]) <- seurat_list[[19]]@"meta.data"$orig.ident



#get the intersect features of all the seurat objects
feature_lists <- lapply(seurat_list, function(obj) Features(obj))
common_features <- Reduce(intersect, feature_lists)

seurat_list_filtered <- lapply(seurat_list, function(obj) {
                                   obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
                                    obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
                                    obj <- subset(obj, subset = nFeature_RNA >= 400 & nCount_RNA >= 1000 & percent.mt < 10 & percent.rb < 10)
                                    obj <- subset(obj, features = common_features)
                                    return(obj)
})
rm(seurat_list)



#Merge all the seurat objects
studies <- c("Fan2020", "Velmeshev2023", "Bakken2021", "Caglayan2023", "Jorstad2023", "Emani2024", "Morabito2021", "Pfisterer2020", "Schirmer2019", "Clarence2025", "Zhu2023", "Batiuk2022", "Steyn2024", "Wang2025", "Maitra2023", "Tran2021", "Pineda2024", "Huuki-Myers2024", "Mannens2024")
data <- merge(seurat_list_filtered[[1]], y = seurat_list_filtered[-1], add.cell.ids = studies)
data <- NormalizeData(data)
data <-FindVariableFeatures(data,selection.method='mean.var.plot')
data <-ScaleData(data,features=VariableFeatures(data))
#run PCA. Select significant PCs based on a scree plot. Look for the last point before the plot becomes flat
data <-RunPCA(data, features = VariableFeatures(data), verbose=F)

png( filename = "Integrate-ElbowPlot.png", width = 20, height = 10, units = "in", res = 300)
ElbowPlot(object = data, ndims = 40)
dev.off()

meta_df <- read.table("meta_AllStudy_Included_stage.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(meta_df) <- meta_df$barcode
meta_df_subset <- meta_df[colnames(data), ]
data <- AddMetaData(data, metadata = meta_df_subset)

saveRDS(data, "Integrate_RunPCA.rds")


N = 23
data <- readRDS("Integrate_RunPCA.rds")
data@reductions$pca2 <- data@reductions$pca
data@reductions$pca2@"cell.embeddings" <- data@reductions$pca2@"cell.embeddings"[,1:N]
data <- RunHarmony(data, group.by.vars = "orig.ident", theta = 2, max.iter.harmony = 20, dims.use = 1:N)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:N)
data <- FindClusters(data, resolution = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
data <- RunUMAP(data, dims = 1:N, reduction = 'harmony')
#saveRDS(data, "Integrate_RunUMAP.rds")

Idents(data) <- "RNA_snn_res.0.5"
data@meta.data$CellType <- data@meta.data$RNA_snn_res.0.5
data@meta.data$CellType <- as.character(data$CellType)

umap_coords <- Embeddings(data, "umap")

cells_in_cluster32 <- WhichCells(data, idents = "32")
cells_to_reassign <- intersect(cells_in_cluster32, rownames(umap_coords)[umap_coords[, "umap_1"] > 7])
data@meta.data$CellType[colnames(data) %in% cells_to_reassign] <- "Oligo"
cells_to_reassign <- intersect(cells_in_cluster32, rownames(umap_coords)[umap_coords[, "umap_1"] <= 7])
data@meta.data$CellType[colnames(data) %in% cells_to_reassign] <- "ExNeu"


cells_in_cluster39 <- WhichCells(data, idents = "39")
cells_to_reassign <- intersect(cells_in_cluster39, rownames(umap_coords)[umap_coords[, "umap_2"] > 5])
data@meta.data$CellType[colnames(data) %in% cells_to_reassign] <- "Astro"
cells_to_reassign <- intersect(cells_in_cluster39, rownames(umap_coords)[umap_coords[, "umap_2"] <= 5])
data@meta.data$CellType[colnames(data) %in% cells_to_reassign] <- "Oligo"

data@meta.data$CellType[data@meta.data$CellType %in% c("4", "41")] <- "Astro"
data@meta.data$CellType[data@meta.data$CellType == "12"] <- "OPC"
data@meta.data$CellType[data@meta.data$CellType %in% c("1", "45")] <- "Oligo"
data@meta.data$CellType[data@meta.data$CellType == "18"] <- "Micro"
data@meta.data$CellType[data@meta.data$CellType == "20"] <- "Unknown"
data@meta.data$CellType[data@meta.data$CellType %in% c("5", "29", "3", "7", "47", "17", "16")] <- "InNeu"
data@meta.data$CellType[data@meta.data$CellType %in% c("0", "2", "6", "8", "9", "11", "13", "14", "19", "21", "22", "23", "24", "26", "27", "28", "30", "31", "37", "38", "40", "42", "46")] <- "ExNeu"
data@meta.data$CellType[data@meta.data$CellType == "25"] <- "Glia-Progenitor"
data@meta.data$CellType[data@meta.data$CellType %in% c("15", "10")] <- "Neu-Progenitor"

data@meta.data$CellType[data@meta.data$CellType %in% c("33", "34", "35", "36", "43", "44")] <- "Unknown"
gmeta <- read.csv("Glia-Meta.csv", header = T, row.names = 1)
data@meta.data$CellType[data@meta.data$barcode %in% gmeta$barcode[gmeta$CellType == "VLMC"]] <- "VLMC"
data@meta.data$CellType[data@meta.data$barcode %in% gmeta$barcode[gmeta$CellType %in% c("Endo1", "Endo2")]] <- "Endo"

data@meta.data$CellType <- factor(data@meta.data$CellType)

Idents(data) <- "CellType"
pdf("Dimplot_0.5_CellType.pdf", width = 13, height = 10)
DimPlot(data, reduction = "umap", label = T)
dev.off()

write.csv(data@meta.data, "Integrate-Meta.csv")

data$CellType <- factor(Idents(data), levels = rev(c("Neu-Progenitor", "InNeu", "ExNeu", "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "VLMC", "Endo", "Unknown")))
Idents(data) <- "CellType"
pdf("DotPlot_All.pdf", width = 15, height = 10)
DotPlot(data, features = c("SOX4", "SNAP25", "GAD1", "SLC17A7", "AQP4", "PDGFRA", "MOG", "P2RY12", "COL1A2", "FLT1"),
                idents = rev(c("Neu-Progenitor", "InNeu", "ExNeu", "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "VLMC", "Endo")),
                        dot.scale = 6)
dev.off()
saveRDS(data, "Integrate_RunUMAP.rds")