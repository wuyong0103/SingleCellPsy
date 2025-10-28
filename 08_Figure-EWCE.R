library(ggplot2)
library(reshape2)
library(tidyverse)
setwd("D:/CellType/downsample500/EWCE")

#Type all level
# read the sd data
sd <- read_tsv("All_Disorder_Type_sd.txt")
sd_long <- sd %>% 
  select(CellType, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "sd", -CellType)

# read the qvalue data
q <- read_tsv("All_Disorder_Type_q.txt")
q_long <- q %>% 
  select(CellType, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "q", -CellType)

# merge sd and qvalue
merged_df <- merge(sd_long, q_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
merged_df <- merged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#1879B6", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
merged_df$CellType <- gsub("_", "-", merged_df$CellType)
merged_df$CellType <- factor(merged_df$CellType,levels=c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC"))
merged_df$group <- gsub("_Exome_All", "-Exome", merged_df$group)
merged_df$group <- factor(merged_df$group, levels = c("BD-Exome", "SCZ-Exome"))

# plot the barplot
ggplot(merged_df, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance), nudge_y = 0.2, size = 5) +
  scale_fill_identity() +
  facet_grid(group ~ ., scales = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  )+
  labs(x = "", y = "s.d. from mean") +
  scale_y_continuous(limits = c(0,3), breaks = c(0,1,2,3))

ggsave(
  filename = "EWCE_All_Type.pdf",
  width = 4.83,
  height = 4.49, 
  units = "in",
  device = "pdf"
)


#Subtype all level
# read the sd data
sd <- read_tsv("All_Disorder_Subtype_sd.txt")
sd_long <- sd %>% 
  select(CellType, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "sd", -CellType)

# read the qvalue data
q <- read_tsv("All_Disorder_Subtype_q.txt")
q_long <- q %>% 
  select(CellType, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "q", -CellType)

# merge sd and qvalue
merged_df <- merge(sd_long, q_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
merged_df <- merged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#1879B6", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
merged_df$CellType <- gsub("_", "-", merged_df$CellType)
merged_df$CellType <- factor(merged_df$CellType,levels = c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                           "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                           "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC"))
merged_df$group <- gsub("_Exome_All", "-Exome", merged_df$group)
merged_df$group <- factor(merged_df$group, levels = c("BD-Exome", "SCZ-Exome"))

# plot the barplot
ggplot(merged_df, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance), nudge_y = 0.2, size = 5) +
  scale_fill_identity() +
  facet_grid(group ~ ., scales = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  )+
  labs(x = "", y = "s.d. from mean") +
  scale_y_continuous(limits = c(0,4), breaks = c(0,2,4))

ggsave(
  filename = "EWCE_All_Subtype.pdf",
  width = 8.81,
  height = 4.49, 
  units = "in",
  device = "pdf"
)



#Type all level
# read the sd data
sd <- read_tsv("EWCE_Type_sd.txt")
sd_long <- sd %>% 
  select(-All) %>% 
  gather(key = "group", value = "sd", -CellType)

# read the qvalue data
q <- read_tsv("EWCE_Type_q.txt")
q_long <- q %>% 
  select(-All) %>% 
  gather(key = "group", value = "q", -CellType)

# merge sd and qvalue
merged_df <- merge(sd_long, q_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
merged_df <- merged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#1879B6", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
merged_bd <- merged_df %>% filter(str_detect(CellType, pattern = "^BD"))
merged_bd$CellType <- gsub("BD_Exome-", "", merged_bd$CellType)
merged_bd$CellType <- gsub("_", "-", merged_bd$CellType)
merged_bd$CellType <- factor(merged_bd$CellType,levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))
merged_bd$group <- factor(merged_bd$group, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))

merged_sz <- merged_df %>% filter(str_detect(CellType, pattern = "^SCZ"))
merged_sz$CellType <- gsub("SCZ_Exome-", "", merged_sz$CellType)
merged_sz$CellType <- gsub("_", "-", merged_sz$CellType)
merged_sz$CellType <- factor(merged_sz$CellType,levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))
merged_sz$group <- factor(merged_sz$group, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))

# 绘制柱状图
ggplot(merged_bd, aes(y = CellType, x = sd_abs, fill = Color)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance), nudge_x = 0.2, size = 5) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
  labs(x = "s.d. from mean", y = "Cell Types") +
  scale_x_continuous(limits = c(0,4), breaks = c(0,2,4))

ggsave(
  filename = "EWCE-BD-stage-type.pdf", 
  width = 9.75,
  height = 3.33, 
  units = "in",
  device = "pdf"
)

ggplot(merged_sz, aes(y = CellType, x = sd_abs, fill = Color)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance), nudge_x = 0.2, size = 5) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
  labs(x = "s.d. from mean", y = "Cell Types") +
  scale_x_continuous(limits = c(0,4), breaks = c(0,2,4))

ggsave(
  filename = "EWCE-SCZ-stage-type.pdf",
  width = 9.75,
  height = 3.33, 
  units = "in",
  device = "pdf"
)



#Subtype all level
# read the sd data
sd <- read_tsv("EWCE_Subtype_sd.txt")
sd_long <- sd %>% 
  select(-All) %>% 
  gather(key = "group", value = "sd", -CellType)

# read the qvalue data
q <- read_tsv("EWCE_Subtype_q.txt")
q_long <- q %>% 
  select(-All) %>% 
  gather(key = "group", value = "q", -CellType)

# merge sd and qvalue
merged_df <- merge(sd_long, q_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
merged_df <- merged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#1879B6", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
merged_bd <- merged_df %>% filter(str_detect(CellType, pattern = "^BD"))
merged_bd$CellType <- gsub("BD_Exome-", "", merged_bd$CellType)
merged_bd$CellType <- gsub("_", "-", merged_bd$CellType)
merged_bd$CellType <- factor(merged_bd$CellType,levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                             "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                             "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))
merged_bd$group <- factor(merged_bd$group, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))

merged_sz <- merged_df %>% filter(str_detect(CellType, pattern = "^SCZ"))
merged_sz$CellType <- gsub("SCZ_Exome-", "", merged_sz$CellType)
merged_sz$CellType <- gsub("_", "-", merged_sz$CellType)
merged_sz$CellType <- factor(merged_sz$CellType,levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                             "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                             "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))
merged_sz$group <- factor(merged_sz$group, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))

ggplot(merged_bd, aes(y = CellType, x = sd_abs, fill = Color)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance), nudge_x = 0.3, size = 5) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
  labs(x = "s.d. from mean", y = "Cell Types") +
  scale_x_continuous(limits = c(0,6), breaks = c(0,2,4,6))

ggsave(
  filename = "EWCE-BD-stage-subtype.pdf", 
  width = 9.75,
  height = 5.66, 
  units = "in",
  device = "pdf"
)

ggplot(merged_sz, aes(y = CellType, x = sd_abs, fill = Color)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance), nudge_x = 0.3, size = 5) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
  labs(x = "s.d. from mean", y = "Cell Types") +
  scale_x_continuous(limits = c(0,6), breaks = c(0,2,4,6))

ggsave(
  filename = "EWCE-SCZ-stage-subtype.pdf", 
  width = 9.75,
  height = 5.66, 
  units = "in",
  device = "pdf"
)