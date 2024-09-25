library(ggplot2)
library(reshape2)
setwd("F:/CellType/Figure")

####################################################################################
#Figure 3A
# 读取数据
ldsc <- read.table("../LDSC-CellType/All_Disorder_Type_Male.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("../MAGMA-CellType/All_Disorder_Type_Male.txt", header = TRUE, row.names = 1, sep = "\t")

# 计算LDSC与MAGMA的P值平均值的对数
mean_log_p <- (-log10(ldsc) -log10(magma))/2

# 转换数据格式
mean_log_p_melted <- melt(as.matrix(mean_log_p))

# 设置显著性阈值
#Bonferroni correction, 8 major cell types, 28 Subtypes, 14 disorder
#P = 0.05/((8+28)*14) = 9.92e-5
threshold <- 9.92e-5

# 合并数据并添加显著性信息
ldsc_melted <- melt(as.matrix(ldsc))
magma_melted <- melt(as.matrix(magma))
combined <- merge(mean_log_p_melted, ldsc_melted, by = c("Var1", "Var2"))
combined <- merge(combined, magma_melted, by = c("Var1", "Var2"))
colnames(combined) <- c("CellType", "Disease", "Mean_Log_P", "LDSC_P", "MAGMA_P")

combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))

# 设置颜色
combined$Color <- with(combined, ifelse(Significance == "both", "#0291BA",
                                        ifelse(Significance == "ldsc", "#F4931A",
                                               ifelse(Significance == "magma", "#BD56A2", "#C8AF81"))))


#调整顺序
combined$CellType <- factor(combined$CellType,levels=c("VASC","GLIALPROG","AST","OPC","OL","MG","IN","ExNeu"))
combined$Disease <- factor(combined$Disease,levels=c("BD","MDD","Intelligence","SCZ","AD","OCT","Anxiety","HoardingSymptoms",
                                                     "TouretteSyndrome","Insomina","Suicide","EatingDisorder","ADHD","ASD"))

# 绘制柱状图
ggplot(combined, aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#377eb8", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ Disease, scales = "fixed", nrow = 1) +
  labs(x = "Mean(-log10(P))", y = "Cell Types") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))

ggsave(
  filename = "Figure3A.pdf", # 名字为文件名称；后缀代表什么格式的图片。
  width = 11.38,
  height = 2.57, 
  units = "in",
  device = "pdf"
)


####################################################################################
#Figure 3B
# 读取数据
ldsc <- read.table("../LDSC-CellType/All_Disorder_Type_Female.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("../MAGMA-CellType/All_Disorder_Type_Female.txt", header = TRUE, row.names = 1, sep = "\t")

# 计算LDSC与MAGMA的P值平均值的对数
mean_log_p <- (-log10(ldsc) -log10(magma))/2

# 转换数据格式
mean_log_p_melted <- melt(as.matrix(mean_log_p))

# 设置显著性阈值
#Bonferroni correction, 8 major cell types, 28 Subtypes, 14 disorder
#P = 0.05/((8+28)*14) = 9.92e-5
threshold <- 9.92e-5

# 合并数据并添加显著性信息
ldsc_melted <- melt(as.matrix(ldsc))
magma_melted <- melt(as.matrix(magma))
combined <- merge(mean_log_p_melted, ldsc_melted, by = c("Var1", "Var2"))
combined <- merge(combined, magma_melted, by = c("Var1", "Var2"))
colnames(combined) <- c("CellType", "Disease", "Mean_Log_P", "LDSC_P", "MAGMA_P")

combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))

# 设置颜色
combined$Color <- with(combined, ifelse(Significance == "both", "#0291BA",
                                        ifelse(Significance == "ldsc", "#F4931A",
                                               ifelse(Significance == "magma", "#BD56A2", "#C8AF81"))))


#调整顺序
combined$CellType <- factor(combined$CellType,levels=c("VASC","GLIALPROG","AST","OPC","OL","MG","IN","ExNeu"))
combined$Disease <- factor(combined$Disease,levels=c("BD","MDD","Intelligence","SCZ","AD","OCT","Anxiety","HoardingSymptoms",
                                                     "TouretteSyndrome","Insomina","Suicide","EatingDisorder","ADHD","ASD"))

# 绘制柱状图
ggplot(combined, aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#DB4D47", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ Disease, scales = "fixed", nrow = 1) +
  labs(x = "Mean(-log10(P))", y = "Cell Types") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))

ggsave(
  filename = "Figure3B.pdf", # 名字为文件名称；后缀代表什么格式的图片。
  width = 11.38,
  height = 2.57, 
  units = "in",
  device = "pdf"
)


####################################################################################
#Figure 3C
# 读取数据
ldsc <- read.table("../LDSC-SubCellType/All_Disorder_Subtype_Male.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("../MAGMA-SubCellType/All_Disorder_Subtype_Male.txt", header = TRUE, row.names = 1, sep = "\t")

# 计算LDSC与MAGMA的P值平均值的对数
mean_log_p <- (-log10(ldsc) -log10(magma))/2

# 转换数据格式
mean_log_p_melted <- melt(as.matrix(mean_log_p))

# 设置显著性阈值
#Bonferroni correction, 8 major cell types, 28 Subtypes, 14 disorder
#P = 0.05/((8+28)*14) = 9.92e-5
threshold <- 9.92e-5

# 合并数据并添加显著性信息
ldsc_melted <- melt(as.matrix(ldsc))
magma_melted <- melt(as.matrix(magma))
combined <- merge(mean_log_p_melted, ldsc_melted, by = c("Var1", "Var2"))
combined <- merge(combined, magma_melted, by = c("Var1", "Var2"))
colnames(combined) <- c("CellType", "Disease", "Mean_Log_P", "LDSC_P", "MAGMA_P")

combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))

# 设置颜色
combined$Color <- with(combined, ifelse(Significance == "both", "#0291BA",
                                        ifelse(Significance == "ldsc", "#F4931A",
                                               ifelse(Significance == "magma", "#BD56A2", "#C8AF81"))))


#调整顺序
combined$CellType <- factor(combined$CellType,levels=rev(c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
                                                           "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
                                                           "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc")))
combined$Disease <- factor(combined$Disease,levels=c("BD","MDD","Intelligence","SCZ","AD","OCT","Anxiety","HoardingSymptoms",
                                                     "TouretteSyndrome","Insomina","Suicide","EatingDisorder","ADHD","ASD"))

# 绘制柱状图
ggplot(combined, aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#377eb8", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ Disease, scales = "fixed", nrow = 1) +
  labs(x = "Mean(-log10(P))", y = "Cell Types") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))

ggsave(
  filename = "Figure3C.pdf", # 名字为文件名称；后缀代表什么格式的图片。
  width = 11.93,
  height = 4.40, 
  units = "in",
  device = "pdf"
)



####################################################################################
#Figure 3D
# 读取数据
ldsc <- read.table("../LDSC-SubCellType/All_Disorder_Subtype_Female.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("../MAGMA-SubCellType/All_Disorder_Subtype_Female.txt", header = TRUE, row.names = 1, sep = "\t")

# 计算LDSC与MAGMA的P值平均值的对数
mean_log_p <- (-log10(ldsc) -log10(magma))/2

# 转换数据格式
mean_log_p_melted <- melt(as.matrix(mean_log_p))

# 设置显著性阈值
#Bonferroni correction, 8 major cell types, 28 Subtypes, 14 disorder
#P = 0.05/((8+28)*14) = 9.92e-5
threshold <- 9.92e-5

# 合并数据并添加显著性信息
ldsc_melted <- melt(as.matrix(ldsc))
magma_melted <- melt(as.matrix(magma))
combined <- merge(mean_log_p_melted, ldsc_melted, by = c("Var1", "Var2"))
combined <- merge(combined, magma_melted, by = c("Var1", "Var2"))
colnames(combined) <- c("CellType", "Disease", "Mean_Log_P", "LDSC_P", "MAGMA_P")

combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))

# 设置颜色
combined$Color <- with(combined, ifelse(Significance == "both", "#0291BA",
                                        ifelse(Significance == "ldsc", "#F4931A",
                                               ifelse(Significance == "magma", "#BD56A2", "#C8AF81"))))


#调整顺序
combined$CellType <- factor(combined$CellType,levels=rev(c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
                                                           "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
                                                           "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc")))
combined$Disease <- factor(combined$Disease,levels=c("BD","MDD","Intelligence","SCZ","AD","OCT","Anxiety","HoardingSymptoms",
                                                     "TouretteSyndrome","Insomina","Suicide","EatingDisorder","ADHD","ASD"))

# 绘制柱状图
ggplot(combined, aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#DB4D47", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ Disease, scales = "fixed", nrow = 1) +
  labs(x = "Mean(-log10(P))", y = "Cell Types") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))

ggsave(
  filename = "Figure3D.pdf", # 名字为文件名称；后缀代表什么格式的图片。
  width = 11.93,
  height = 4.40, 
  units = "in",
  device = "pdf"
)