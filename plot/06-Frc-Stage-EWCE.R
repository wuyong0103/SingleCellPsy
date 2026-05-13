library(ggplot2)
library(tidyverse)
library(gridExtra)
setwd("D:/CellType/Figures-20260429")

#====================================================================================================================
#SubType level
data <- read_tsv('06-Frc_Sublineage-Stage_05.tsv')
data$CellType <- gsub("_", "-", data$CellType)
# 3. 数据预处理
df <- data %>% 
  mutate(logP = -log10(q)) %>%
  mutate(logP = ifelse(is.infinite(logP), -log10(0.0001), logP))
df$Tissue <- sub("Migration", "NeuronalMigration", df$Tissue)
df$Tissue <- sub("Grundtypus", "6-layer Grundtypus", df$Tissue)

df$CellType <- factor(df$CellType,levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                 "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                 "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))
df$Tissue <- factor(df$Tissue, levels=c("Neurogenesis", "NeuronalMigration", "Astrogliogenesis", "Synaptogenesis", "6-layer Grundtypus", "Oligodendrogenesis", "Myelination", "SynapticPruning"))

# 4. 绘制叠加图
pss5 <- df %>%
  filter(Trait == "SCZ") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )

psb5 <- df %>%
  filter(Trait == "BD") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
  
pdf("06-Frc-Subtype-Stage-EWCE-05.pdf", width = 10.5, height = 11.00)
grid.arrange(psb5, pss5, ncol = 1)
dev.off()



#====================================================================================================================
#Type level
data <- read_tsv('06-Frc_Lineage-Stage_05.tsv')
data$CellType <- gsub("_", "-", data$CellType)
# 3. 数据预处理
df <- data %>% 
  mutate(logP = -log10(q)) %>%
  mutate(logP = ifelse(is.infinite(logP), -log10(0.0001), logP))
df$Tissue <- sub("Migration", "NeuronalMigration", df$Tissue)
df$Tissue <- sub("Grundtypus", "6-layer Grundtypus", df$Tissue)

df$CellType <- factor(df$CellType,levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))
df$Tissue <- factor(df$Tissue, levels=c("Neurogenesis", "NeuronalMigration", "Astrogliogenesis", "Synaptogenesis", "6-layer Grundtypus", "Oligodendrogenesis", "Myelination", "SynapticPruning"))

# 4. 绘制叠加图
pts5 <- df %>%
  filter(Trait == "SCZ") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
  
  
ptb5 <- df %>%
  filter(Trait == "BD") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
  
pdf("06-Frc-Type-Stage-EWCE-05.pdf", width = 10.5, height = 7.17)
grid.arrange(ptb5, pts5, ncol = 1)
dev.off()






#====================================================================================================================
#SubType level
data <- read_tsv('06-Frc_Sublineage-Stage_01.tsv')
data$CellType <- gsub("_", "-", data$CellType)
# 3. 数据预处理
df <- data %>% 
  mutate(logP = -log10(q)) %>%
  mutate(logP = ifelse(is.infinite(logP), -log10(0.0001), logP))
df$Tissue <- sub("Migration", "NeuronalMigration", df$Tissue)
df$Tissue <- sub("Grundtypus", "6-layer Grundtypus", df$Tissue)

df$CellType <- factor(df$CellType,levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                 "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                 "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))
df$Tissue <- factor(df$Tissue, levels=c("Neurogenesis", "NeuronalMigration", "Astrogliogenesis", "Synaptogenesis", "6-layer Grundtypus", "Oligodendrogenesis", "Myelination", "SynapticPruning"))

# 4. 绘制叠加图
pss1 <- df %>%
  filter(Trait == "SCZ") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )

psb1 <- df %>%
  filter(Trait == "BD") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
  
pdf("06-Frc-Subtype-Stage-EWCE-01.pdf", width = 10.5, height = 11.00)
grid.arrange(psb1, pss1, ncol = 1)
dev.off()



#====================================================================================================================
#Type level
data <- read_tsv('06-Frc_Lineage-Stage_01.tsv')
data$CellType <- gsub("_", "-", data$CellType)
# 3. 数据预处理
df <- data %>% 
  mutate(logP = -log10(q)) %>%
  mutate(logP = ifelse(is.infinite(logP), -log10(0.0001), logP))
df$Tissue <- sub("Migration", "NeuronalMigration", df$Tissue)
df$Tissue <- sub("Grundtypus", "6-layer Grundtypus", df$Tissue)

df$CellType <- factor(df$CellType,levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))
df$Tissue <- factor(df$Tissue, levels=c("Neurogenesis", "NeuronalMigration", "Astrogliogenesis", "Synaptogenesis", "6-layer Grundtypus", "Oligodendrogenesis", "Myelination", "SynapticPruning"))

# 4. 绘制叠加图
pts1 <- df %>%
  filter(Trait == "SCZ") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
  
  
ptb1 <- df %>%
  filter(Trait == "BD") %>%
  ggplot(aes(x = CellType, y = logP)) +
  # 柱状图：表示均值
  geom_bar(stat = "identity", fill = "#FAB95B", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 4)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](q-value)),
  ) +
  scale_y_continuous(
    breaks = seq(0, 4, by = 1), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#547792", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
  
pdf("06-Frc-Type-Stage-EWCE-01.pdf", width = 10.5, height = 7.17)
grid.arrange(ptb1, pts1, ncol = 1)
dev.off()