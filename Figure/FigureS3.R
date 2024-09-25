library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
library("ggplot2")
setwd("F:/CellType/Figure")

#Figure S3A-S3D
# 读取数据
ldsc_type_male <- read.table("../LDSC-CellType/All_Disorder_Type_Male.txt", row.names = 1, header = T)
ldsc_type_female <- read.table("../LDSC-CellType/All_Disorder_Type_Female.txt", row.names = 1, header = T)
magma_type_male <- read.table("../MAGMA-CellType/All_Disorder_Type_Male.txt", row.names = 1, header = T)
magma_type_female <- read.table("../MAGMA-CellType/All_Disorder_Type_Female.txt", row.names = 1, header = T)

loci <- read.table("../Others/Disorder_loci.txt", row.names = 1, header = T)

#定义热图疾病绘制顺序
nc <- rownames(loci)

# 创建样本数量柱状图
loci <- anno_barplot(loci, border = FALSE, add_numbers = TRUE, axis = FALSE, bar_width = 0.8,
                     numbers_gp = gpar(fontsize = 10, col = "red"), numbers_rot = 90, numbers_offset = unit(0.3, "mm"))

#定义热图细胞类型绘制顺序
nrt <- c("ExNeu", "IN", "MG", "OL", "OPC", "AST", "GLIALPROG", "VASC")

####################################################################################
#Figure S3A
#MAGMA Cell Type Male
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
magma_type_male <- magma_type_male[nrt, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- magma_type_male[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}
mtm <- Heatmap(as.matrix(-log10(magma_type_male)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "MAGMA cell type", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'GnBu'))(4)), 
              width = ncol(magma_type_male)*unit(6, "mm"), height = nrow(magma_type_male)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), 
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(mtm)

####################################################################################
#Figure S3B
#LDSC Cell Type Male
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
ldsc_type_male <- ldsc_type_male[nrt, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- ldsc_type_male[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}
ltm <- Heatmap(as.matrix(-log10(ldsc_type_male)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "LDSC cell type", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'PuRd'))(4)),
              width = ncol(ldsc_type_male)*unit(6, "mm"), height = nrow(ldsc_type_male)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci),
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(ltm)

####################################################################################
#Figure S3C
#MAGMA Cell Type Female
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
magma_type_female <- magma_type_female[nrt, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- magma_type_female[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

mtf <- Heatmap(as.matrix(-log10(magma_type_female)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "MAGMA cell type", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'GnBu'))(4)), 
              width = ncol(magma_type_female)*unit(6, "mm"), height = nrow(magma_type_female)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), 
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(mtf)

####################################################################################
#Figure S3D
#LDSC Cell Type Female
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
ldsc_type_female <- ldsc_type_female[nrt, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- ldsc_type_female[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

ltf <- Heatmap(as.matrix(-log10(ldsc_type_female)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "LDSC cell type", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'PuRd'))(4)),
              width = ncol(ldsc_type_female)*unit(6, "mm"), height = nrow(ldsc_type_female)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), 
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(ltf)



#Figure S3E-S3H
# 读取数据
ldsc_subtype_male <- read.table("../LDSC-SubCellType/All_Disorder_Subtype_Male.txt", row.names = 1, header = T)
ldsc_subtype_female <- read.table("../LDSC-SubCellType/All_Disorder_Subtype_Female.txt", row.names = 1, header = T)
magma_subtype_male <- read.table("../MAGMA-SubCellType/All_Disorder_Subtype_Male.txt", row.names = 1, header = T)
magma_subtype_female <- read.table("../MAGMA-SubCellType/All_Disorder_Subtype_Female.txt", row.names = 1, header = T)

loci <- read.table("../Others/Disorder_loci.txt", row.names = 1, header = T)

#定义热图疾病绘制顺序
nc <- rownames(loci)

# 创建样本数量柱状图
loci <- anno_barplot(loci, border = FALSE, add_numbers = TRUE, axis = FALSE, bar_width = 0.8,
                     numbers_gp = gpar(fontsize = 10, col = "red"), numbers_rot = 90, numbers_offset = unit(0.3, "mm"))

#定义热图细胞类型绘制顺序
nrs <- c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
         "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
         "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc")

####################################################################################
#Figure S3E
#MAGMA Cell Subtype Male
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
magma_subtype_male <- magma_subtype_male[nrs, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- magma_subtype_male[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

# 创建基因功能注释数据框
cell_type <- data.frame(
  subtype = rownames(magma_subtype_male),
  type = c(rep("ExNeu", 8), rep("IN", 12), "MG", "OL", "OPC", "AST", "AST", "GLIALPROG", "VASC", "VASC")
)
rownames(cell_type) <- cell_type$subtype

# 创建左侧注释
left_annotation <- rowAnnotation(
  CellType = cell_type$type,
  col = list(CellType = c("ExNeu" = "#1b9e77", "IN" = "#d95f02", "MG" = "#7570b3", "OL"="#e7298a", "OPC" = "#66a61e",
                          "AST"="#e6ab02", "GLIALPROG"="#a6761d", "VASC"="#377eb8"))
)

msm <- Heatmap(as.matrix(-log10(magma_subtype_male)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "MAGMA cell subtype", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'GnBu'))(4)),
              width = ncol(magma_subtype_male)*unit(6, "mm"), height = nrow(magma_subtype_male)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), rect_gp = gpar(col = "grey", lwd = 1), left_annotation = left_annotation,
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(msm)



####################################################################################
#Figure S1F
#LDSC Cell Subtype Male
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
ldsc_subtype_male <- ldsc_subtype_male[nrs, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- ldsc_subtype_male[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

lsm <- Heatmap(as.matrix(-log10(ldsc_subtype_male)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "LDSC cell subtype", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'PuRd'))(4)),
              width = ncol(ldsc_subtype_male)*unit(6, "mm"), height = nrow(ldsc_subtype_male)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), rect_gp = gpar(col = "grey", lwd = 1), left_annotation = left_annotation,
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(lsm)

####################################################################################
#Figure S3G
#MAGMA Cell Subtype Female
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
magma_subtype_female <- magma_subtype_female[nrs, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- magma_subtype_female[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

# 创建基因功能注释数据框
cell_type <- data.frame(
  subtype = rownames(magma_subtype_female),
  type = c(rep("ExNeu", 8), rep("IN", 12), "MG", "OL", "OPC", "AST", "AST", "GLIALPROG", "VASC", "VASC")
)
rownames(cell_type) <- cell_type$subtype

# 创建左侧注释
left_annotation <- rowAnnotation(
  CellType = cell_type$type,
  col = list(CellType = c("ExNeu" = "#1b9e77", "IN" = "#d95f02", "MG" = "#7570b3", "OL"="#e7298a", "OPC" = "#66a61e",
                          "AST"="#e6ab02", "GLIALPROG"="#a6761d", "VASC"="#377eb8"))
)

msf <- Heatmap(as.matrix(-log10(magma_subtype_female)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "MAGMA cell subtype", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'GnBu'))(4)),
              width = ncol(magma_subtype_female)*unit(6, "mm"), height = nrow(magma_subtype_female)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), rect_gp = gpar(col = "grey", lwd = 1), left_annotation = left_annotation,
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(msf)

####################################################################################
#Figure S3H
#LDSC Cell Subtype Female
# 创建热图并使用cell_fun
# 自定义cell_fun以添加星号
ldsc_subtype_female <- ldsc_subtype_female[nrs, nc]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = fill, col = NA, fontsize = 16))
  expr_value <- ldsc_subtype_female[i, j]
  if (expr_value < 9.92e-5) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 16))
  } else if (expr_value < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 16))
  }
}

lsf <- Heatmap(as.matrix(-log10(ldsc_subtype_female)), name = "-log10(P)", cluster_rows = FALSE, cluster_columns = FALSE, 
              column_title = "LDSC cell subtype", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
              col = colorRamp2(c(0, 3, 6, 9), colorRampPalette(brewer.pal(9,'PuRd'))(4)),
              width = ncol(ldsc_subtype_female)*unit(6, "mm"), height = nrow(ldsc_subtype_female)*unit(6, "mm"),
              top_annotation = HeatmapAnnotation(GWASLoci = loci), rect_gp = gpar(col = "grey", lwd = 1), left_annotation = left_annotation,
              heatmap_legend_param = list(at = c(0, 3, 6, 9), labels = c(0, 3, 6, 9), legend_height = unit(3.5, "cm")),
              cell_fun = cell_fun)

# 绘制热图
draw(lsf)