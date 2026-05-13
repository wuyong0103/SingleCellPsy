library(ggplot2)
library(tidyverse)
library(patchwork)
setwd("D:/CellType/Figures-20260429")

#====================================================================================================================
#SubType level
data <- read_tsv('01-All-subtype.txt')
# 3. 数据预处理
df <- data %>%
  mutate(logP = -log10(pvalue)) %>%
  # 关键步骤：按 Trait 和 cell_type 分组，计算不同 method 的均值
  group_by(Trait, cell_type) %>%
  mutate(mean_logP = mean(logP)) %>%
  ungroup() %>%
  # 按均值排序，使图表更有序
  mutate(cell_type = reorder(cell_type, -mean_logP))
df$Trait <- sub("(\\d{4}).*", "\\1", df$Trait)

df$cell_type <- factor(df$cell_type,levels=c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                               "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                               "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC"))

# 4. 绘制叠加图
ps4 <- df %>% filter(str_starts(Trait, "SCZ")) %>%
  filter(Trait != "SCZ2026") %>%
  ggplot(aes(x = cell_type, y = logP)) +
    # 柱状图：表示均值
    stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
    # 绘制显著线 (P = 0.05)
    # linetype = "dashed" 设置为虚线，color 设置为红色
    geom_hline(yintercept = -log10(0.00179), linetype = "dashed", color = "red", linewidth = 0.8) +
  
    # 散点图：展示不同 method
    geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  
    # 分面设置：
    # ncol = 1 确保每行一个 Trait
    # strip.position = "right" 将标签放在右侧
    facet_wrap(~Trait, ncol = 1, strip.position = "right") +
    scale_y_continuous(
      breaks = seq(0, 6, by = 2), # 设置刻度间隔为 2 (0, 2, 4, ...)
      expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
    ) + 
  
    # 标签和主题
    labs(
      x = "Cell Type",
      y = expression(-log[10](p-value)),
      color = "Method"
    ) +
    theme_bw() +
    theme(
      # --- 去掉所有网格线 ---
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      # 自定义分面标签（Strip）
      strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
      strip.text = element_text(color = "white", face = "bold", size = 11),        # 标签文字加粗
      panel.spacing = unit(1, "lines"),                                  # 增加分面间距
      legend.position = "none"
    )
#ggsave(filename = "01-SCZ-DiffGWAS-Subtype.pdf", width = 7.00, height = 6.25, units = "in",device = "pdf")

ps2 <- df %>% filter(str_starts(Trait, "BD")) %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.00179), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Trait, ncol = 1, strip.position = "right") +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 11),        # 标签文字加粗
    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
#ggsave(filename = "01-BD-DiffGWAS-Subtype.pdf", width = 7.00, height = 6.25, units = "in",device = "pdf")


ps1 <- df %>% filter(Trait %in% c("SCZ2022", "BD2025")) %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.00179), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Trait, ncol = 1, strip.position = "right") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 11),        # 标签文字加粗
    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )



#====================================================================================================================
#Type level
data <- read_tsv('01-All-type.txt')
# 3. 数据预处理
df <- data %>%
  mutate(logP = -log10(pvalue)) %>%
  # 关键步骤：按 Trait 和 cell_type 分组，计算不同 method 的均值
  group_by(Trait, cell_type) %>%
  mutate(mean_logP = mean(logP)) %>%
  ungroup() %>%
  # 按均值排序，使图表更有序
  mutate(cell_type = reorder(cell_type, -mean_logP))
df$Trait <- sub("(\\d{4}).*", "\\1", df$Trait)

df$cell_type <- factor(df$cell_type, levels=c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC"))

# 4. 绘制叠加图
pt4 <- df %>% filter(str_starts(Trait, "SCZ")) %>%
  filter(Trait != "SCZ2026") %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Trait, ncol = 1, strip.position = "right") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 11),        # 标签文字加粗
    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
#ggsave(filename = "01-SCZ-DiffGWAS-Type.pdf", width = 3.76, height = 6.25, units = "in",device = "pdf")


# 4. 绘制叠加图
pt2 <- df %>% filter(str_starts(Trait, "BD")) %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Trait, ncol = 1, strip.position = "right") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 11),        # 标签文字加粗
    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
#ggsave(filename = "01-BD-DiffGWAS-Type.pdf", width = 3.76, height = 6.25, units = "in",device = "pdf")

pt1 <- df %>% filter(Trait %in% c("SCZ2022", "BD2025")) %>%
  ggplot(aes(x = cell_type, y = logP)) +
  # 柱状图：表示均值
  stat_summary(fun = mean, geom = "bar", fill = "gray90", color = "gray70", width = 0.7) +
  
  # 绘制显著线 (P = 0.05)
  # linetype = "dashed" 设置为虚线，color 设置为红色
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # 散点图：展示不同 method
  geom_point(aes(color = method), size = 1.8, alpha = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2), # 设置刻度间隔为 2 (0, 2, 4, ...)
    expand = expansion(mult = c(0, 0.1))    # 去除坐标轴与绘图区之间的空隙
  ) + 
  
  # 分面设置：
  # ncol = 1 确保每行一个 Trait
  # strip.position = "right" 将标签放在右侧
  facet_wrap(~Trait, ncol = 1, strip.position = "right") +
  
  # 标签和主题
  labs(
    x = "Cell Type",
    y = expression(-log[10](p-value)),
    color = "Method"
  ) +
  theme_bw() +
  theme(
    # --- 去掉所有网格线 ---
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # 自定义分面标签（Strip）
    strip.background = element_rect(fill = "#1C9E78", color = "black"), # 背景绿色
    strip.text = element_text(color = "white", face = "bold", size = 11),        # 标签文字加粗
    panel.spacing = unit(1, "lines"),                                  # 增加分面间距
    legend.position = "none"
  )
#ggsave(filename = "01-SCZ2022-Type.pdf", width = 3.76, height = 3.41, units = "in",device = "pdf")


#组图
combined_plot <- pt1 | ps1
final_plot <- combined_plot + 
  plot_layout(widths = c(3, 7))

# 保存为 PDF
ggsave("01-NewGWAS.pdf", final_plot, width = 10.76, height = 5.3)


#组图
combined_plot <- (pt2 / pt4) | (ps2 / ps4)
final_plot <- combined_plot + 
  plot_layout(widths = c(3, 7))

# 保存为 PDF
ggsave("01-AllGWAS.pdf", final_plot, width = 10.76, height = 11.5)