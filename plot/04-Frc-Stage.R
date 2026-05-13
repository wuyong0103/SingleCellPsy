library(ggplot2)
library(tidyverse)
library(gridExtra)
setwd("D:/CellType/Figures-20260429")

#====================================================================================================================
#SubType level
data <- read_tsv('04-Frc-Subtype-Stage.txt')
# 3. 数据预处理
df <- data %>%
  filter(Trait %in% c("BD2025OConnel", "SCZ2022PGC3")) %>%
  mutate(logP = -log10(pvalue)) %>%
  # 关键步骤：按 Trait 和 cell_type 分组，计算不同 method 的均值
  group_by(Trait, stage, cell_type) %>%
  mutate(mean_logP = mean(logP)) %>%
  ungroup() %>%
  # 按均值排序，使图表更有序
  mutate(cell_type = reorder(cell_type, -mean_logP))
df$Trait <- sub("(\\d{4}).*", "\\1", df$Trait)
df$stage <- sub("Migration", "NeuronalMigration", df$stage)
df$stage <- sub("Grundtypus", "6-layer Grundtypus", df$stage)

df$cell_type <- factor(df$cell_type,levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                 "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                 "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))
df$stage <- factor(df$stage, levels=c("Neurogenesis", "NeuronalMigration", "Astrogliogenesis", "Synaptogenesis", "6-layer Grundtypus", "Oligodendrogenesis", "Myelination", "SynapticPruning"))

# 4. 绘制叠加图
p1 <- df %>% filter(Trait == "BD2025") %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.00179), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 9)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~stage, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  scale_y_continuous(
    breaks = seq(0, 9, by = 3), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )

p2 <- df %>% filter(Trait == "SCZ2022") %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.00179), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 9)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~stage, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  scale_y_continuous(
    breaks = seq(0, 9, by = 3), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
pdf("04-Frc-Subtype-Stage.pdf", width = 10.5, height = 11.00)
grid.arrange(p1, p2, ncol = 1)
dev.off()



#====================================================================================================================
#Type level
data <- read_tsv('04-Frc-Type-Stage.txt')
# 3. 数据预处理
df <- data %>%
  filter(Trait %in% c("BD2025OConnel", "SCZ2022PGC3")) %>%
  mutate(logP = -log10(pvalue)) %>%
  # 关键步骤：按 Trait 和 cell_type 分组，计算不同 method 的均值
  group_by(Trait, stage, cell_type) %>%
  mutate(mean_logP = mean(logP)) %>%
  ungroup() %>%
  # 按均值排序，使图表更有序
  mutate(cell_type = reorder(cell_type, -mean_logP))
df$Trait <- sub("(\\d{4}).*", "\\1", df$Trait)
df$stage <- sub("Migration", "NeuronalMigration", df$stage)
df$stage <- sub("Grundtypus", "6-layer Grundtypus", df$stage)

df$cell_type <- factor(df$cell_type,levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))
df$stage <- factor(df$stage, levels=c("Neurogenesis", "NeuronalMigration", "Astrogliogenesis", "Synaptogenesis", "6-layer Grundtypus", "Oligodendrogenesis", "Myelination", "SynapticPruning"))

# 4. 绘制叠加图
p1 <- df %>% filter(Trait == "BD2025") %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 12)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~stage, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  scale_y_continuous(
    breaks = seq(0, 12, by = 3), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )

p2 <- df %>% filter(Trait == "SCZ2022") %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  # 对调 X 轴和 Y 轴
  coord_flip(ylim = c(0, 12)) +
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~stage, nrow = 1, strip.position = "top") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  scale_y_continuous(
    breaks = seq(0, 12, by = 3), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 8),        # 标签文字加粗
    #    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
pdf("04-Frc-Type-Stage.pdf", width = 10.50, height = 7.17)
grid.arrange(p1, p2, ncol = 1)
dev.off()