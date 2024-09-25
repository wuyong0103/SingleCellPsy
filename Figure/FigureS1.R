library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
library("ggplot2")
setwd("F:/CellType/Figure")

# 读取数据
ldsc_type <- read.table("../LDSC-CellType/All_Disorder_Type.txt", row.names = 1, header = T)
magma_type <- read.table("../MAGMA-CellType/All_Disorder_Type.txt", row.names = 1, header = T)
ldsc_subtype <- read.table("../LDSC-SubCellType/All_Disorder_Subtype.txt", row.names = 1, header = T)
magma_subtype <- read.table("../MAGMA-SubCellType/All_Disorder_Subtype.txt", row.names = 1, header = T)
loci <- read.table("../Others/Disorder_loci.txt", row.names = 1, header = T)

#定义热图疾病绘制顺序
nc <- rownames(loci)

# 创建样本数量柱状图
loci <- anno_barplot(loci, border = FALSE, add_numbers = TRUE, axis = FALSE, bar_width = 0.8,
                     numbers_gp = gpar(fontsize = 10, col = "red"), numbers_rot = 90, numbers_offset = unit(0.3, "mm"))

#定义热图细胞类型绘制顺序
nrt <- c("ExNeu", "IN", "MG", "OL", "OPC", "AST", "GLIALPROG", "VASC")
nrs <- c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
         "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
         "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc")

####################################################################################
#Figure S1A
#MAGMA Cell Type
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
magma_type <- magma_type[nrt, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- magma_type[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}
mt <- Heatmap(as.matrix(-log10(magma_type)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "MAGMA cell type", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'GnBu'))(4)),
			  width = ncol(magma_type)*unit(6, "mm"), height = nrow(magma_type)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), 
			  heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(mt)

####################################################################################
#Figure S1C
#MAGMA Cell Subtype
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
magma_subtype <- magma_subtype[nrs, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- magma_subtype[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

# 创建基因功能注释数据框
cell_type <- data.frame(
  subtype = rownames(magma_subtype),
  type = c(rep("ExNeu", 8), rep("IN", 12), "MG", "OL", "OPC", "AST", "AST", "GLIALPROG", "VASC", "VASC")
)
rownames(cell_type) <- cell_type$subtype

# 创建左侧注释
left_annotation <- rowAnnotation(
  CellType = cell_type$type,
  col = list(CellType = c("ExNeu" = "#1b9e77", "IN" = "#d95f02", "MG" = "#7570b3", "OL"="#e7298a", "OPC" = "#66a61e",
                          "AST"="#e6ab02", "GLIALPROG"="#a6761d", "VASC"="#377eb8"))
)

ms <- Heatmap(as.matrix(-log10(magma_subtype)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "MAGMA cell subtype", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'GnBu'))(4)),
			  width = ncol(magma_subtype)*unit(6, "mm"), height = nrow(magma_subtype)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), left_annotation = left_annotation,
			  heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(ms)


####################################################################################
#Figure S1B
#LDSC Cell Type
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
ldsc_type <- ldsc_type[nrt, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- ldsc_type[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

lt <- Heatmap(as.matrix(-log10(ldsc_type)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "LDSC cell type", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'PuRd'))(4)),
			  width = ncol(ldsc_type)*unit(6, "mm"), height = nrow(ldsc_type)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), 
			  heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(lt)

####################################################################################
#Figure S1D
#LDSC Cell Subtype
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
ldsc_subtype <- ldsc_subtype[nrs, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- ldsc_subtype[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

ls <- Heatmap(as.matrix(-log10(ldsc_subtype)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "LDSC cell subtype", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'PuRd'))(4)),
			  width = ncol(ldsc_subtype)*unit(6, "mm"), height = nrow(ldsc_subtype)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), left_annotation = left_annotation,
			  heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(ls)
