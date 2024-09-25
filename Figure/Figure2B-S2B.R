library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
setwd("F:/CellType/Figure")

####################################################################################
#Figure 2A
#读取数据
ldsc <- read.table("../LDSC-SubCellType/LDSC_Subtype.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("../MAGMA-SubCellType/MAGMA_Subtype.txt", header = TRUE, row.names = 1, sep = "\t")
rnames <- rownames(ldsc)
cnames <- colnames(ldsc)
magma <- magma[rnames, cnames]

# 计算LDSC与MAGMA的P值平均值的对数
mean_log_p <- (-log10(ldsc) -log10(magma))/2
mean_log_p[is.na(ldsc) | is.na(ldsc)] <- NA

row_type_name <- c("ExNeu","IN","MG","OL","OPC","AST","GLIALPROG","VASC")
row_subtype_name <- c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
                      "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
                      "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc")

# 创建基因功能注释数据框
cell_type <- data.frame(
  subtype = row_subtype_name,
  type = c(rep("ExNeu", 8), rep("IN", 12), "MG", "OL", "OPC", "AST", "AST", "GLIALPROG", "VASC", "VASC")
)
rownames(cell_type) <- cell_type$subtype

# 创建左侧注释
left_annotation <- rowAnnotation(
  CellType = cell_type$type,
  col = list(CellType = c("ExNeu" = "#1b9e77", "IN" = "#d95f02", "MG" = "#7570b3", "OL"="#e7298a", "OPC" = "#66a61e",
                          "AST"="#e6ab02", "GLIALPROG"="#a6761d", "VASC"="#377eb8"))
)
####################################################################################
#SCZ
scz_mean <- mean_log_p[grep("^SCZ", rownames(mean_log_p)), 2:9]
scz_ldsc <- ldsc[grep("^SCZ", rownames(ldsc)), 2:9]
scz_magma <- magma[grep("^SCZ", rownames(magma)), 2:9]
rownames(scz_mean) <- gsub("[A-Za-z0-9]+-","", rownames(scz_mean))
colnames(scz_mean) <- gsub("Mean_SubCellType_", "", colnames(scz_mean))
rownames(scz_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(scz_ldsc))
colnames(scz_ldsc) <- gsub("Mean_SubCellType_", "", colnames(scz_ldsc))
rownames(scz_magma) <- gsub("[A-Za-z0-9]+-","", rownames(scz_magma))
colnames(scz_magma) <- gsub("Mean_SubCellType_", "", colnames(scz_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(scz_ldsc), ncol = ncol(scz_ldsc))
rownames(marker_matrix) <- rownames(scz_ldsc)
colnames(marker_matrix) <- colnames(scz_ldsc)

# 添加标识
for (i in 1:nrow(scz_ldsc)) {
  for (j in 1:ncol(scz_ldsc)) {
    if (!is.na(scz_ldsc[i, j]) & !is.na(scz_magma[i, j])) {
      if (scz_ldsc[i, j] < 1.736E-4 & scz_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (scz_ldsc[i, j] < 1.736E-4 | scz_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
scz <- Heatmap(scz_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(scz_mean)*unit(6, "mm"), 
               height = nrow(scz_mean)*unit(6, "mm"),column_title = "SCZ cell subtype",
               left_annotation = left_annotation, cell_fun = text_fun)

####################################################################################
#BD
bd_mean <- mean_log_p[grep("^BD", rownames(mean_log_p)), 2:9]
bd_ldsc <- ldsc[grep("^BD", rownames(ldsc)), 2:9]
bd_magma <- magma[grep("^BD", rownames(magma)), 2:9]
rownames(bd_mean) <- gsub("[A-Za-z0-9]+-","", rownames(bd_mean))
colnames(bd_mean) <- gsub("Mean_SubCellType_", "", colnames(bd_mean))
rownames(bd_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(bd_ldsc))
colnames(bd_ldsc) <- gsub("Mean_SubCellType_", "", colnames(bd_ldsc))
rownames(bd_magma) <- gsub("[A-Za-z0-9]+-","", rownames(bd_magma))
colnames(bd_magma) <- gsub("Mean_SubCellType_", "", colnames(bd_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(bd_ldsc), ncol = ncol(bd_ldsc))
rownames(marker_matrix) <- rownames(bd_ldsc)
colnames(marker_matrix) <- colnames(bd_ldsc)

# 添加标识
for (i in 1:nrow(bd_ldsc)) {
  for (j in 1:ncol(bd_ldsc)) {
    if (!is.na(bd_ldsc[i, j]) & !is.na(bd_magma[i, j])) {
      if (bd_ldsc[i, j] < 1.736E-4 & bd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (bd_ldsc[i, j] < 1.736E-4 | bd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
bd <- Heatmap(bd_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(bd_mean)*unit(6, "mm"), 
              height = nrow(bd_mean)*unit(6, "mm"),column_title = "BD cell subtype",
              cell_fun = text_fun)


####################################################################################
#MDD
mdd_mean <- mean_log_p[grep("^MDD", rownames(mean_log_p)), 2:9]
mdd_ldsc <- ldsc[grep("^MDD", rownames(ldsc)), 2:9]
mdd_magma <- magma[grep("^MDD", rownames(magma)), 2:9]
rownames(mdd_mean) <- gsub("[A-Za-z0-9]+-","", rownames(mdd_mean))
colnames(mdd_mean) <- gsub("Mean_SubCellType_", "", colnames(mdd_mean))
rownames(mdd_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(mdd_ldsc))
colnames(mdd_ldsc) <- gsub("Mean_SubCellType_", "", colnames(mdd_ldsc))
rownames(mdd_magma) <- gsub("[A-Za-z0-9]+-","", rownames(mdd_magma))
colnames(mdd_magma) <- gsub("Mean_SubCellType_", "", colnames(mdd_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(mdd_ldsc), ncol = ncol(mdd_ldsc))
rownames(marker_matrix) <- rownames(mdd_ldsc)
colnames(marker_matrix) <- colnames(mdd_ldsc)

# 添加标识
for (i in 1:nrow(mdd_ldsc)) {
  for (j in 1:ncol(mdd_ldsc)) {
    if (!is.na(mdd_ldsc[i, j]) & !is.na(mdd_magma[i, j])) {
      if (mdd_ldsc[i, j] < 1.736E-4 & mdd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (mdd_ldsc[i, j] < 1.736E-4 | mdd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
mdd <- Heatmap(mdd_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(mdd_mean)*unit(6, "mm"), 
               height = nrow(mdd_mean)*unit(6, "mm"),column_title = "MDD cell subtype",
               cell_fun = text_fun)


####################################################################################
#INT
int_mean <- mean_log_p[grep("^Inte", rownames(mean_log_p)), 2:9]
int_ldsc <- ldsc[grep("^Inte", rownames(ldsc)), 2:9]
int_magma <- magma[grep("^Inte", rownames(magma)), 2:9]
rownames(int_mean) <- gsub("[A-Za-z0-9]+-","", rownames(int_mean))
colnames(int_mean) <- gsub("Mean_SubCellType_", "", colnames(int_mean))
rownames(int_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(int_ldsc))
colnames(int_ldsc) <- gsub("Mean_SubCellType_", "", colnames(int_ldsc))
rownames(int_magma) <- gsub("[A-Za-z0-9]+-","", rownames(int_magma))
colnames(int_magma) <- gsub("Mean_SubCellType_", "", colnames(int_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(int_ldsc), ncol = ncol(int_ldsc))
rownames(marker_matrix) <- rownames(int_ldsc)
colnames(marker_matrix) <- colnames(int_ldsc)

# 添加标识
for (i in 1:nrow(int_ldsc)) {
  for (j in 1:ncol(int_ldsc)) {
    if (!is.na(int_ldsc[i, j]) & !is.na(int_magma[i, j])) {
      if (int_ldsc[i, j] < 1.736E-4 & int_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (int_ldsc[i, j] < 1.736E-4 | int_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
int <- Heatmap(int_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(int_mean)*unit(6, "mm"), 
               height = nrow(int_mean)*unit(6, "mm"),column_title = "INT cell subtype",
               cell_fun = text_fun)


####################################################################################
#AD
ad_mean <- mean_log_p[grep("^AD2022", rownames(mean_log_p)), 2:9]
ad_ldsc <- ldsc[grep("^AD2022", rownames(ldsc)), 2:9]
ad_magma <- magma[grep("^AD2022", rownames(magma)), 2:9]
rownames(ad_mean) <- gsub("[A-Za-z0-9]+-","", rownames(ad_mean))
colnames(ad_mean) <- gsub("Mean_SubCellType_", "", colnames(ad_mean))
rownames(ad_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(ad_ldsc))
colnames(ad_ldsc) <- gsub("Mean_SubCellType_", "", colnames(ad_ldsc))
rownames(ad_magma) <- gsub("[A-Za-z0-9]+-","", rownames(ad_magma))
colnames(ad_magma) <- gsub("Mean_SubCellType_", "", colnames(ad_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(ad_ldsc), ncol = ncol(ad_ldsc))
rownames(marker_matrix) <- rownames(ad_ldsc)
colnames(marker_matrix) <- colnames(ad_ldsc)

# 添加标识
for (i in 1:nrow(ad_ldsc)) {
  for (j in 1:ncol(ad_ldsc)) {
    if (!is.na(ad_ldsc[i, j]) & !is.na(ad_magma[i, j])) {
      if (ad_ldsc[i, j] < 1.736E-4 & ad_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (ad_ldsc[i, j] < 1.736E-4 | ad_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
ad <- Heatmap(ad_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(ad_mean)*unit(6, "mm"), 
              height = nrow(ad_mean)*unit(6, "mm"),column_title = "AD cell subtype",
              left_annotation = left_annotation, cell_fun = text_fun)


####################################################################################
#oct
oct_mean <- mean_log_p[grep("^OCT", rownames(mean_log_p)), 2:9]
oct_ldsc <- ldsc[grep("^OCT", rownames(ldsc)), 2:9]
oct_magma <- magma[grep("^OCT", rownames(magma)), 2:9]
rownames(oct_mean) <- gsub("[A-Za-z0-9]+-","", rownames(oct_mean))
colnames(oct_mean) <- gsub("Mean_SubCellType_", "", colnames(oct_mean))
rownames(oct_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(oct_ldsc))
colnames(oct_ldsc) <- gsub("Mean_SubCellType_", "", colnames(oct_ldsc))
rownames(oct_magma) <- gsub("[A-Za-z0-9]+-","", rownames(oct_magma))
colnames(oct_magma) <- gsub("Mean_SubCellType_", "", colnames(oct_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(oct_ldsc), ncol = ncol(oct_ldsc))
rownames(marker_matrix) <- rownames(oct_ldsc)
colnames(marker_matrix) <- colnames(oct_ldsc)

# 添加标识
for (i in 1:nrow(oct_ldsc)) {
  for (j in 1:ncol(oct_ldsc)) {
    if (!is.na(oct_ldsc[i, j]) & !is.na(oct_magma[i, j])) {
      if (oct_ldsc[i, j] < 1.736E-4 & oct_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (oct_ldsc[i, j] < 1.736E-4 | oct_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
oct <- Heatmap(oct_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(oct_mean)*unit(6, "mm"), 
               height = nrow(oct_mean)*unit(6, "mm"),column_title = "OCT cell subtype",
               cell_fun = text_fun)


####################################################################################
#an
an_mean <- mean_log_p[grep("^Anxiety", rownames(mean_log_p)), 2:9]
an_ldsc <- ldsc[grep("^Anxiety", rownames(ldsc)), 2:9]
an_magma <- magma[grep("^Anxiety", rownames(magma)), 2:9]
rownames(an_mean) <- gsub("[A-Za-z0-9]+-","", rownames(an_mean))
colnames(an_mean) <- gsub("Mean_SubCellType_", "", colnames(an_mean))
rownames(an_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(an_ldsc))
colnames(an_ldsc) <- gsub("Mean_SubCellType_", "", colnames(an_ldsc))
rownames(an_magma) <- gsub("[A-Za-z0-9]+-","", rownames(an_magma))
colnames(an_magma) <- gsub("Mean_SubCellType_", "", colnames(an_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(an_ldsc), ncol = ncol(an_ldsc))
rownames(marker_matrix) <- rownames(an_ldsc)
colnames(marker_matrix) <- colnames(an_ldsc)

# 添加标识
for (i in 1:nrow(an_ldsc)) {
  for (j in 1:ncol(an_ldsc)) {
    if (!is.na(an_ldsc[i, j]) & !is.na(an_magma[i, j])) {
      if (an_ldsc[i, j] < 1.736E-4 & an_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (an_ldsc[i, j] < 1.736E-4 | an_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
an <- Heatmap(an_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(an_mean)*unit(6, "mm"), 
              height = nrow(an_mean)*unit(6, "mm"),column_title = "AN cell subtype",
              cell_fun = text_fun)


####################################################################################
#HS
hs_mean <- mean_log_p[grep("^Hoarding", rownames(mean_log_p)), 2:9]
hs_ldsc <- ldsc[grep("^Hoarding", rownames(ldsc)), 2:9]
hs_magma <- magma[grep("^Hoarding", rownames(magma)), 2:9]
rownames(hs_mean) <- gsub("[A-Za-z0-9]+-","", rownames(hs_mean))
colnames(hs_mean) <- gsub("Mean_SubCellType_", "", colnames(hs_mean))
rownames(hs_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(hs_ldsc))
colnames(hs_ldsc) <- gsub("Mean_SubCellType_", "", colnames(hs_ldsc))
rownames(hs_magma) <- gsub("[A-Za-z0-9]+-","", rownames(hs_magma))
colnames(hs_magma) <- gsub("Mean_SubCellType_", "", colnames(hs_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(hs_ldsc), ncol = ncol(hs_ldsc))
rownames(marker_matrix) <- rownames(hs_ldsc)
colnames(marker_matrix) <- colnames(hs_ldsc)

# 添加标识
for (i in 1:nrow(hs_ldsc)) {
  for (j in 1:ncol(hs_ldsc)) {
    if (!is.na(hs_ldsc[i, j]) & !is.na(hs_magma[i, j])) {
      if (hs_ldsc[i, j] < 1.736E-4 & hs_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (hs_ldsc[i, j] < 1.736E-4 | hs_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
hs <- Heatmap(hs_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(hs_mean)*unit(6, "mm"), 
              height = nrow(hs_mean)*unit(6, "mm"),column_title = "HS cell subtype",
              cell_fun = text_fun)

####################################################################################
#TS
ts_mean <- mean_log_p[grep("^Tourette", rownames(mean_log_p)), 2:9]
ts_ldsc <- ldsc[grep("^Tourette", rownames(ldsc)), 2:9]
ts_magma <- magma[grep("^Tourette", rownames(magma)), 2:9]
rownames(ts_mean) <- gsub("[A-Za-z0-9]+-","", rownames(ts_mean))
colnames(ts_mean) <- gsub("Mean_SubCellType_", "", colnames(ts_mean))
rownames(ts_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(ts_ldsc))
colnames(ts_ldsc) <- gsub("Mean_SubCellType_", "", colnames(ts_ldsc))
rownames(ts_magma) <- gsub("[A-Za-z0-9]+-","", rownames(ts_magma))
colnames(ts_magma) <- gsub("Mean_SubCellType_", "", colnames(ts_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(ts_ldsc), ncol = ncol(ts_ldsc))
rownames(marker_matrix) <- rownames(ts_ldsc)
colnames(marker_matrix) <- colnames(ts_ldsc)

# 添加标识
for (i in 1:nrow(ts_ldsc)) {
  for (j in 1:ncol(ts_ldsc)) {
    if (!is.na(ts_ldsc[i, j]) & !is.na(ts_magma[i, j])) {
      if (ts_ldsc[i, j] < 1.736E-4 & ts_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (ts_ldsc[i, j] < 1.736E-4 | ts_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
ts <- Heatmap(ts_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(ts_mean)*unit(6, "mm"), 
              height = nrow(ts_mean)*unit(6, "mm"),column_title = "TS cell subtype",
              cell_fun = text_fun)


####################################################################################
#INS
ins_mean <- mean_log_p[grep("^Insomina", rownames(mean_log_p)), 2:9]
ins_ldsc <- ldsc[grep("^Insomina", rownames(ldsc)), 2:9]
ins_magma <- magma[grep("^Insomina", rownames(magma)), 2:9]
rownames(ins_mean) <- gsub("[A-Za-z0-9]+-","", rownames(ins_mean))
colnames(ins_mean) <- gsub("Mean_SubCellType_", "", colnames(ins_mean))
rownames(ins_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(ins_ldsc))
colnames(ins_ldsc) <- gsub("Mean_SubCellType_", "", colnames(ins_ldsc))
rownames(ins_magma) <- gsub("[A-Za-z0-9]+-","", rownames(ins_magma))
colnames(ins_magma) <- gsub("Mean_SubCellType_", "", colnames(ins_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(ins_ldsc), ncol = ncol(ins_ldsc))
rownames(marker_matrix) <- rownames(ins_ldsc)
colnames(marker_matrix) <- colnames(ins_ldsc)

# 添加标识
for (i in 1:nrow(ins_ldsc)) {
  for (j in 1:ncol(ins_ldsc)) {
    if (!is.na(ins_ldsc[i, j]) & !is.na(ins_magma[i, j])) {
      if (ins_ldsc[i, j] < 1.736E-4 & ins_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (ins_ldsc[i, j] < 1.736E-4 | ins_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
ins <- Heatmap(ins_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(ins_mean)*unit(6, "mm"), 
               height = nrow(ins_mean)*unit(6, "mm"),column_title = "INS cell subtype",
               cell_fun = text_fun)


####################################################################################
#su
su_mean <- mean_log_p[grep("^Suicide", rownames(mean_log_p)), 2:9]
su_ldsc <- ldsc[grep("^Suicide", rownames(ldsc)), 2:9]
su_magma <- magma[grep("^Suicide", rownames(magma)), 2:9]
rownames(su_mean) <- gsub("[A-Za-z0-9]+-","", rownames(su_mean))
colnames(su_mean) <- gsub("Mean_SubCellType_", "", colnames(su_mean))
rownames(su_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(su_ldsc))
colnames(su_ldsc) <- gsub("Mean_SubCellType_", "", colnames(su_ldsc))
rownames(su_magma) <- gsub("[A-Za-z0-9]+-","", rownames(su_magma))
colnames(su_magma) <- gsub("Mean_SubCellType_", "", colnames(su_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(su_ldsc), ncol = ncol(su_ldsc))
rownames(marker_matrix) <- rownames(su_ldsc)
colnames(marker_matrix) <- colnames(su_ldsc)

# 添加标识
for (i in 1:nrow(su_ldsc)) {
  for (j in 1:ncol(su_ldsc)) {
    if (!is.na(su_ldsc[i, j]) & !is.na(su_magma[i, j])) {
      if (su_ldsc[i, j] < 1.736E-4 & su_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (su_ldsc[i, j] < 1.736E-4 | su_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
su <- Heatmap(su_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(su_mean)*unit(6, "mm"), 
              height = nrow(su_mean)*unit(6, "mm"),column_title = "SU cell subtype",
              cell_fun = text_fun)


####################################################################################
#ED
ed_mean <- mean_log_p[grep("^Eating", rownames(mean_log_p)), 2:9]
ed_ldsc <- ldsc[grep("^Eating", rownames(ldsc)), 2:9]
ed_magma <- magma[grep("^Eating", rownames(magma)), 2:9]
rownames(ed_mean) <- gsub("[A-Za-z0-9]+-","", rownames(ed_mean))
colnames(ed_mean) <- gsub("Mean_SubCellType_", "", colnames(ed_mean))
rownames(ed_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(ed_ldsc))
colnames(ed_ldsc) <- gsub("Mean_SubCellType_", "", colnames(ed_ldsc))
rownames(ed_magma) <- gsub("[A-Za-z0-9]+-","", rownames(ed_magma))
colnames(ed_magma) <- gsub("Mean_SubCellType_", "", colnames(ed_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(ed_ldsc), ncol = ncol(ed_ldsc))
rownames(marker_matrix) <- rownames(ed_ldsc)
colnames(marker_matrix) <- colnames(ed_ldsc)

# 添加标识
for (i in 1:nrow(ed_ldsc)) {
  for (j in 1:ncol(ed_ldsc)) {
    if (!is.na(ed_ldsc[i, j]) & !is.na(ed_magma[i, j])) {
      if (ed_ldsc[i, j] < 1.736E-4 & ed_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (ed_ldsc[i, j] < 1.736E-4 | ed_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
ed <- Heatmap(ed_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
              name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
              width = ncol(ed_mean)*unit(6, "mm"), 
              height = nrow(ed_mean)*unit(6, "mm"),column_title = "ED cell subtype",
              cell_fun = text_fun)


####################################################################################
#ADHD
adhd_mean <- mean_log_p[grep("^ADHD", rownames(mean_log_p)), 2:9]
adhd_ldsc <- ldsc[grep("^ADHD", rownames(ldsc)), 2:9]
adhd_magma <- magma[grep("^ADHD", rownames(magma)), 2:9]
rownames(adhd_mean) <- gsub("[A-Za-z0-9]+-","", rownames(adhd_mean))
colnames(adhd_mean) <- gsub("Mean_SubCellType_", "", colnames(adhd_mean))
rownames(adhd_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(adhd_ldsc))
colnames(adhd_ldsc) <- gsub("Mean_SubCellType_", "", colnames(adhd_ldsc))
rownames(adhd_magma) <- gsub("[A-Za-z0-9]+-","", rownames(adhd_magma))
colnames(adhd_magma) <- gsub("Mean_SubCellType_", "", colnames(adhd_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(adhd_ldsc), ncol = ncol(adhd_ldsc))
rownames(marker_matrix) <- rownames(adhd_ldsc)
colnames(marker_matrix) <- colnames(adhd_ldsc)

# 添加标识
for (i in 1:nrow(adhd_ldsc)) {
  for (j in 1:ncol(adhd_ldsc)) {
    if (!is.na(adhd_ldsc[i, j]) & !is.na(adhd_magma[i, j])) {
      if (adhd_ldsc[i, j] < 1.736E-4 & adhd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (adhd_ldsc[i, j] < 1.736E-4 | adhd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
adhd <- Heatmap(adhd_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                width = ncol(adhd_mean)*unit(6, "mm"), 
                height = nrow(adhd_mean)*unit(6, "mm"),column_title = "ADHD cell subtype",
                cell_fun = text_fun)

####################################################################################
#ASD
asd_mean <- mean_log_p[grep("^ASD", rownames(mean_log_p)), 2:9]
asd_ldsc <- ldsc[grep("^ASD", rownames(ldsc)), 2:9]
asd_magma <- magma[grep("^ASD", rownames(magma)), 2:9]
rownames(asd_mean) <- gsub("[A-Za-z0-9]+-","", rownames(asd_mean))
colnames(asd_mean) <- gsub("Mean_SubCellType_", "", colnames(asd_mean))
rownames(asd_ldsc) <- gsub("[A-Za-z0-9]+-","", rownames(asd_ldsc))
colnames(asd_ldsc) <- gsub("Mean_SubCellType_", "", colnames(asd_ldsc))
rownames(asd_magma) <- gsub("[A-Za-z0-9]+-","", rownames(asd_magma))
colnames(asd_magma) <- gsub("Mean_SubCellType_", "", colnames(asd_magma))

# 创建一个标识矩阵
marker_matrix <- matrix("", nrow = nrow(asd_ldsc), ncol = ncol(asd_ldsc))
rownames(marker_matrix) <- rownames(asd_ldsc)
colnames(marker_matrix) <- colnames(asd_ldsc)

# 添加标识
for (i in 1:nrow(asd_ldsc)) {
  for (j in 1:ncol(asd_ldsc)) {
    if (!is.na(asd_ldsc[i, j]) & !is.na(asd_magma[i, j])) {
      if (asd_ldsc[i, j] < 1.736E-4 & asd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "**"
      } else if (asd_ldsc[i, j] < 1.736E-4 | asd_magma[i, j] < 1.736E-4) {
        marker_matrix[i, j] <- "*"
      }
    }
  }
}
marker_matrix <- marker_matrix[row_subtype_name,]

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(marker_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

# 绘制热图
asd <- Heatmap(asd_mean[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "Mean(-log10(P))",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(asd_mean)*unit(6, "mm"), 
               height = nrow(asd_mean)*unit(6, "mm"),column_title = "ASD cell subtype",
               cell_fun = text_fun)


#Figure2B
scz+bd+mdd+int

#FigureS2B
ad+oct+an+hs+ts+ins+su+ed+adhd+asd