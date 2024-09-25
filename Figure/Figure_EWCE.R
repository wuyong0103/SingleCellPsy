library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
setwd("F:/CellType/Figure")

####################################################################################################
#Subtype all level
# read the sd data
sd <- read_tsv("../EWCE/All_Disorder_Subtype_sd.txt")
sd_long <- sd %>% 
  select(CellType, ASD_Diff_All, BD_Diff_All, SCZ_Diff_All, ASD_SFARI_All, ASD_Exome_All, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "sd", ASD_Diff_All:SCZ_Exome_All)

# read the qvalue data
q <- read_tsv("../EWCE/All_Disorder_Subtype_q.txt")
q_long <- q %>% 
  select(CellType, ASD_Diff_All, BD_Diff_All, SCZ_Diff_All, ASD_SFARI_All, ASD_Exome_All, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "q", ASD_Diff_All:SCZ_Exome_All)

# merge sd and qvalue
merged_df <- merge(sd_long, q_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
merged_df <- merged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#1879B6", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
merged_df$CellType <- factor(merged_df$CellType,levels = c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
                                                           "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
                                                           "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc"))
merged_df$group <- factor(merged_df$group, levels = c("ASD_SFARI_All", "ASD_Exome_All", "BD_Exome_All", "SCZ_Exome_All", 
                                                      "ASD_Diff_All", "BD_Diff_All", "SCZ_Diff_All"))

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
  scale_y_continuous(limits = c(0,15), breaks = c(0,5,10,15))

ggsave(
  filename = "EWCE_All_Subtype.pdf",
  width = 7.56,
  height = 6.32, 
  units = "in",
  device = "pdf"
)

#Subtype Exome
merged_ex <- merged_df %>% filter(group %in% c("ASD_SFARI_All", "ASD_Exome_All", "BD_Exome_All", "SCZ_Exome_All"))
# plot the barplot
ggplot(merged_ex, aes(x = CellType, y = sd_abs, fill = Color)) +
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
  scale_y_continuous(limits = c(0,12), breaks = c(0,4,8,12))
ggsave(
  filename = "EWCE_All_Exome_Subtype.pdf",
  width = 7.56,
  height = 4.89, 
  units = "in",
  device = "pdf"
)

#Subtype Differentail expression
merged_diff <- merged_df %>% filter(group %in% c("ASD_Diff_All", "BD_Diff_All", "SCZ_Diff_All"))
# plot the barplot
ggplot(merged_diff, aes(x = CellType, y = sd_abs, fill = Color)) +
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
  scale_y_continuous(limits = c(0,15), breaks = c(0,5,10,15))
ggsave(
  filename = "EWCE_All_Diff_Subtype.pdf",
  width = 7.56,
  height = 4.39, 
  units = "in",
  device = "pdf"
)


####################################################################################################
#Type all level
# read the sd data
sd <- read_tsv("../EWCE/All_Disorder_Type_sd.txt")
sd_long <- sd %>% 
  select(CellType, ASD_Diff_All, BD_Diff_All, SCZ_Diff_All, ASD_SFARI_All, ASD_Exome_All, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "sd", ASD_Diff_All:SCZ_Exome_All)

# read the qvalue data
q <- read_tsv("../EWCE/All_Disorder_Type_q.txt")
q_long <- q %>% 
  select(CellType, ASD_Diff_All, BD_Diff_All, SCZ_Diff_All, ASD_SFARI_All, ASD_Exome_All, BD_Exome_All, SCZ_Exome_All) %>% 
  gather(key = "group", value = "q", ASD_Diff_All:SCZ_Exome_All)

# merge sd and qvalue
merged_df <- merge(sd_long, q_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
merged_df <- merged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#1879B6", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
merged_df$CellType <- factor(merged_df$CellType,levels=rev(c("VASC","GLIALPROG","AST","OPC","OL","MG","IN","ExNeu")))
merged_df$group <- factor(merged_df$group, levels = c("ASD_SFARI_All", "ASD_Exome_All", "BD_Exome_All", "SCZ_Exome_All", 
                                                      "ASD_Diff_All", "BD_Diff_All", "SCZ_Diff_All"))

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
  scale_y_continuous(limits = c(0,12), breaks = c(0,4,8,12))

ggsave(
  filename = "EWCE_All_Type.pdf",
  width = 5.50,
  height = 6.32, 
  units = "in",
  device = "pdf"
)


#Type Exome
merged_ex <- merged_df %>% filter(group %in% c("ASD_SFARI_All", "ASD_Exome_All", "BD_Exome_All", "SCZ_Exome_All"))
# plot the barplot
ggplot(merged_ex, aes(x = CellType, y = sd_abs, fill = Color)) +
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
  scale_y_continuous(limits = c(0,8), breaks = c(0,4,8))
ggsave(
  filename = "EWCE_All_Exome_Type.pdf",
  width = 5.50,
  height = 4.89, 
  units = "in",
  device = "pdf"
)


#Type Differentail expression
merged_diff <- merged_df %>% filter(group %in% c("ASD_Diff_All", "BD_Diff_All", "SCZ_Diff_All"))
# plot the barplot
ggplot(merged_diff, aes(x = CellType, y = sd_abs, fill = Color)) +
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
  scale_y_continuous(limits = c(0,12), breaks = c(0,4,8,12))
ggsave(
  filename = "EWCE_All_Diff_Type.pdf",
  width = 7.56,
  height = 4.39, 
  units = "in",
  device = "pdf"
)



####################################################################################################
#Subtype Gender level
# read the sd data
msd <- read_tsv("../EWCE/All_Disorder_Subtype_Male_sd.txt")
msd_long <- msd %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "sd", ASD_Diff:SCZ_Exome)

# read the qvalue data
mq <- read_tsv("../EWCE/All_Disorder_Subtype_Male_q.txt")
mq_long <- mq %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "q", ASD_Diff:SCZ_Exome)

# merge sd and qvalue
mmerged_df <- merge(msd_long, mq_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
mmerged_df <- mmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#377eb8", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
mmerged_df$sex <- "Male"

fsd <- read_tsv("../EWCE/All_Disorder_Subtype_Female_sd.txt")
fsd_long <- fsd %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "sd", ASD_Diff:SCZ_Exome)

# read the qvalue data
fq <- read_tsv("../EWCE/All_Disorder_Subtype_Female_q.txt")
fq_long <- fq %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "q", ASD_Diff:SCZ_Exome)

# merge sd and qvalue
fmerged_df <- merge(fsd_long, fq_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
fmerged_df <- fmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#DB4D47", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
fmerged_df$sex <- "Female"

merged_df <- rbind(mmerged_df,fmerged_df)

merged_df$CellType <- factor(merged_df$CellType,levels = c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
                                                           "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
                                                           "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc"))
merged_df$group <- factor(merged_df$group, levels = c("ASD_SFARI", "ASD_Exome", "BD_Exome", "SCZ_Exome", 
                                                      "ASD_Diff", "BD_Diff", "SCZ_Diff"))

# plot the barplot
ggplot(merged_df, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance, y = sd + 0.01), position = position_dodge(width = 0.9), vjust = 0) +
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
  scale_y_continuous(limits = c(0,15), breaks = c(0,5,10,15))

ggsave(
  filename = "EWCE_All_Subtype_Sex.pdf",
  width = 7.56,
  height = 6.32, 
  units = "in",
  device = "pdf"
)

#Subtype Exome
merged_ex <- merged_df %>% filter(group %in% c("ASD_SFARI", "ASD_Exome", "BD_Exome", "SCZ_Exome"))
# plot the barplot
ggplot(merged_ex, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance, y = sd + 0.01), position = position_dodge(width = 0.9), vjust = 0) +
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
  scale_y_continuous(limits = c(0,12), breaks = c(0,4,8,12))
ggsave(
  filename = "EWCE_All_Exome_Subtype_Sex.pdf",
  width = 7.56,
  height = 4.89, 
  units = "in",
  device = "pdf"
)

#Subtype Differentail expression
merged_diff <- merged_df %>% filter(group %in% c("ASD_Diff", "BD_Diff", "SCZ_Diff"))
# plot the barplot
ggplot(merged_diff, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance, y = sd + 0.01), position = position_dodge(width = 0.9), vjust = 0) +
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
  scale_y_continuous(limits = c(0,15), breaks = c(0,5,10,15))
ggsave(
  filename = "EWCE_All_Diff_Subtype_Sex.pdf",
  width = 7.56,
  height = 4.39, 
  units = "in",
  device = "pdf"
)



####################################################################################################
#Type Gender level
# read the sd data
msd <- read_tsv("../EWCE/All_Disorder_Type_Male_sd.txt")
msd_long <- msd %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "sd", ASD_Diff:SCZ_Exome)

# read the qvalue data
mq <- read_tsv("../EWCE/All_Disorder_Type_Male_q.txt")
mq_long <- mq %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "q", ASD_Diff:SCZ_Exome)

# merge sd and qvalue
mmerged_df <- merge(msd_long, mq_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
mmerged_df <- mmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#377eb8", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
mmerged_df$sex <- "Male"

fsd <- read_tsv("../EWCE/All_Disorder_Type_Female_sd.txt")
fsd_long <- fsd %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "sd", ASD_Diff:SCZ_Exome)

# read the qvalue data
fq <- read_tsv("../EWCE/All_Disorder_Type_Female_q.txt")
fq_long <- fq %>% 
  select(CellType, ASD_Diff, BD_Diff, SCZ_Diff, ASD_SFARI, ASD_Exome, BD_Exome, SCZ_Exome) %>% 
  gather(key = "group", value = "q", ASD_Diff:SCZ_Exome)

# merge sd and qvalue
fmerged_df <- merge(fsd_long, fq_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
fmerged_df <- fmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#DB4D47", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
fmerged_df$sex <- "Female"

merged_df <- rbind(mmerged_df,fmerged_df)

merged_df$CellType <- factor(merged_df$CellType,levels = rev(c("VASC","GLIALPROG","AST","OPC","OL","MG","IN","ExNeu")))
merged_df$group <- factor(merged_df$group, levels = c("ASD_SFARI", "ASD_Exome", "BD_Exome", "SCZ_Exome", 
                                                      "ASD_Diff", "BD_Diff", "SCZ_Diff"))

# plot the barplot
ggplot(merged_df, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance, y = sd + 0.01), position = position_dodge(width = 0.9), vjust = 0) +
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
  scale_y_continuous(limits = c(0,10), breaks = c(0,5,10))

ggsave(
  filename = "EWCE_All_Type_Sex.pdf",
  width = 7.56,
  height = 6.32, 
  units = "in",
  device = "pdf"
)

#Subtype Exome
merged_ex <- merged_df %>% filter(group %in% c("ASD_SFARI", "ASD_Exome", "BD_Exome", "SCZ_Exome"))
# plot the barplot
ggplot(merged_ex, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance, y = sd + 0.01), position = position_dodge(width = 0.9), vjust = 0) +
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
  scale_y_continuous(limits = c(0,8), breaks = c(0,4,8))
ggsave(
  filename = "EWCE_All_Exome_Type_Sex.pdf",
  width = 7.56,
  height = 4.89, 
  units = "in",
  device = "pdf"
)

#Subtype Differentail expression
merged_diff <- merged_df %>% filter(group %in% c("ASD_Diff", "BD_Diff", "SCZ_Diff"))
# plot the barplot
ggplot(merged_diff, aes(x = CellType, y = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  geom_text(aes(label = Significance, y = sd + 0.01), position = position_dodge(width = 0.9), vjust = 0) +
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
  scale_y_continuous(limits = c(0,10), breaks = c(0,5,10))
ggsave(
  filename = "EWCE_All_Diff_Type_Sex.pdf",
  width = 7.56,
  height = 4.39, 
  units = "in",
  device = "pdf"
)





####################################################################################################
#plot the complexheatmap across different developmental stages
ssd <- read.table("../EWCE/EWCE_Subtype_sd.txt", header = T)
sq <- read.table("../EWCE/EWCE_Subtype_q.txt", header = T)
rownames(ssd) <- ssd$CellType
rownames(sq) <- sq$CellType

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

# 定义自定义标记函数
text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(mark_matrix[i, j], x, y, gp = gpar(fontsize = 16))
}

##########################################################################################
#ASD Exome subtype
asd_ex_ssd <- ssd[grep("^ASD_Exome", ssd$CellType), 5:12]
rownames(asd_ex_ssd) <- gsub("ASD_Exome-","", rownames(asd_ex_ssd))
asd_ex_sq <- sq[grep("^ASD_Exome", sq$CellType), 5:12]
rownames(asd_ex_sq) <- gsub("ASD_Exome-","", rownames(asd_ex_sq))

# 添加标识
mark_positions <- asd_ex_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_ex_sq ), ncol = ncol(asd_ex_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_ex_sq)
colnames(mark_matrix) <- colnames(asd_ex_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_ex_subtype <- Heatmap(asd_ex_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
               width = ncol(asd_ex_ssd)*unit(6, "mm"), 
               height = nrow(asd_ex_ssd)*unit(6, "mm"),column_title = "ASD Exome subtype",
               cell_fun = text_fun)


##########################################################################################
#asd SFARI subtype
asd_sfari_ssd <- ssd[grep("^ASD_SFARI", ssd$CellType), 5:12]
rownames(asd_sfari_ssd) <- gsub("ASD_SFARI-","", rownames(asd_sfari_ssd))
asd_sfari_sq <- sq[grep("^ASD_SFARI", sq$CellType), 5:12]
rownames(asd_sfari_sq) <- gsub("ASD_SFARI-","", rownames(asd_sfari_sq))

# 添加标识
mark_positions <- asd_sfari_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_sfari_sq ), ncol = ncol(asd_sfari_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_sfari_sq)
colnames(mark_matrix) <- colnames(asd_sfari_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_sfari_subtype <- Heatmap(asd_sfari_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                       name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                       width = ncol(asd_sfari_ssd)*unit(6, "mm"), 
                       height = nrow(asd_sfari_ssd)*unit(6, "mm"),column_title = "ASD SFARI subtype",
                       left_annotation = left_annotation, cell_fun = text_fun)


##########################################################################################
#BD Exome subtype
bd_ex_ssd <- ssd[grep("^BD_Exome", ssd$CellType), 5:12]
rownames(bd_ex_ssd) <- gsub("BD_Exome-","", rownames(bd_ex_ssd))
bd_ex_sq <- sq[grep("^BD_Exome", sq$CellType), 5:12]
rownames(bd_ex_sq) <- gsub("BD_Exome-","", rownames(bd_ex_sq))

# 添加标识
mark_positions <- bd_ex_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_ex_sq ), ncol = ncol(bd_ex_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_ex_sq)
colnames(mark_matrix) <- colnames(bd_ex_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
bd_ex_subtype <- Heatmap(bd_ex_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(bd_ex_ssd)*unit(6, "mm"), 
                          height = nrow(bd_ex_ssd)*unit(6, "mm"),column_title = "BD Exome subtype",
                          cell_fun = text_fun)


##########################################################################################
#scz Exome subtype
scz_ex_ssd <- ssd[grep("^SCZ_Exome", ssd$CellType), 5:12]
rownames(scz_ex_ssd) <- gsub("SCZ_Exome-","", rownames(scz_ex_ssd))
scz_ex_sq <- sq[grep("^SCZ_Exome", sq$CellType), 5:12]
rownames(scz_ex_sq) <- gsub("SCZ_Exome-","", rownames(scz_ex_sq))

# 添加标识
mark_positions <- scz_ex_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_ex_sq ), ncol = ncol(scz_ex_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_ex_sq)
colnames(mark_matrix) <- colnames(scz_ex_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
scz_ex_subtype <- Heatmap(scz_ex_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(scz_ex_ssd)*unit(6, "mm"), 
                          height = nrow(scz_ex_ssd)*unit(6, "mm"),column_title = "SCZ Exome subtype",
                          cell_fun = text_fun)



#EWCE_All_Exome_Stages_Subtype
asd_sfari_subtype + asd_ex_subtype + bd_ex_subtype + scz_ex_subtype



##########################################################################################
#ASD Exome subtype Male
asd_ex_ssd_male <- ssd[grep("^ASD_Exome", ssd$CellType), 13:20]
rownames(asd_ex_ssd_male) <- gsub("ASD_Exome-","", rownames(asd_ex_ssd_male))
asd_ex_sq_male <- sq[grep("^ASD_Exome", sq$CellType), 13:20]
rownames(asd_ex_sq_male) <- gsub("ASD_Exome-","", rownames(asd_ex_sq_male))

# 添加标识
mark_positions <- asd_ex_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_ex_sq_male ), ncol = ncol(asd_ex_sq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_ex_sq_male)
colnames(mark_matrix) <- colnames(asd_ex_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_ex_subtype_male <- Heatmap(asd_ex_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(asd_ex_ssd_male)*unit(6, "mm"), 
                          height = nrow(asd_ex_ssd_male)*unit(6, "mm"),column_title = "ASD Exome subtype",
                          cell_fun = text_fun)


##########################################################################################
#ASD SFARI subtype Male
asd_sfari_ssd_male <- ssd[grep("^ASD_SFARI", ssd$CellType), 13:20]
rownames(asd_sfari_ssd_male) <- gsub("ASD_SFARI-","", rownames(asd_sfari_ssd_male))
asd_sfari_sq_male <- sq[grep("^ASD_SFARI", sq$CellType), 13:20]
rownames(asd_sfari_sq_male) <- gsub("ASD_SFARI-","", rownames(asd_sfari_sq_male))

# 添加标识
mark_positions <- asd_sfari_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_sfari_sq_male), ncol = ncol(asd_sfari_sq_male))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_sfari_sq_male)
colnames(mark_matrix) <- colnames(asd_sfari_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_sfari_subtype_male <- Heatmap(asd_sfari_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                             name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                             width = ncol(asd_sfari_ssd_male)*unit(6, "mm"), 
                             height = nrow(asd_sfari_ssd_male)*unit(6, "mm"),column_title = "ASD SFARI subtype",
                             left_annotation = left_annotation, cell_fun = text_fun)


##########################################################################################
#BD Exome subtype Male
bd_ex_ssd_male <- ssd[grep("^BD_Exome", ssd$CellType), 13:20]
rownames(bd_ex_ssd_male) <- gsub("BD_Exome-","", rownames(bd_ex_ssd_male))
bd_ex_sq_male <- sq[grep("^BD_Exome", sq$CellType), 13:20]
rownames(bd_ex_sq_male) <- gsub("BD_Exome-","", rownames(bd_ex_sq_male))

# 添加标识
mark_positions <- bd_ex_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_ex_sq_male ), ncol = ncol(bd_ex_sq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_ex_sq_male)
colnames(mark_matrix) <- colnames(bd_ex_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
bd_ex_subtype_male <- Heatmap(bd_ex_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                         name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                         width = ncol(bd_ex_ssd_male)*unit(6, "mm"), 
                         height = nrow(bd_ex_ssd_male)*unit(6, "mm"),column_title = "BD Exome subtype",
                         cell_fun = text_fun)


##########################################################################################
#scz Exome subtype Male
scz_ex_ssd_male <- ssd[grep("^SCZ_Exome", ssd$CellType), 13:20]
rownames(scz_ex_ssd_male) <- gsub("SCZ_Exome-","", rownames(scz_ex_ssd_male))
scz_ex_sq_male <- sq[grep("^SCZ_Exome", sq$CellType), 13:20]
rownames(scz_ex_sq_male) <- gsub("SCZ_Exome-","", rownames(scz_ex_sq_male))

# 添加标识
mark_positions <- scz_ex_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_ex_sq_male ), ncol = ncol(scz_ex_sq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_ex_sq_male)
colnames(mark_matrix) <- colnames(scz_ex_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
scz_ex_subtype_male <- Heatmap(scz_ex_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(scz_ex_ssd_male)*unit(6, "mm"), 
                          height = nrow(scz_ex_ssd_male)*unit(6, "mm"),column_title = "SCZ Exome subtype",
                          cell_fun = text_fun)



#EWCE_All_Exome_Stages_Subtype_Male
asd_sfari_subtype_male + asd_ex_subtype_male + bd_ex_subtype_male + scz_ex_subtype_male



##########################################################################################
#ASD Exome subtype female
asd_ex_ssd_female <- ssd[grep("^ASD_Exome", ssd$CellType), 21:28]
rownames(asd_ex_ssd_female) <- gsub("ASD_Exome-","", rownames(asd_ex_ssd_female))
asd_ex_sq_female <- sq[grep("^ASD_Exome", sq$CellType), 21:28]
rownames(asd_ex_sq_female) <- gsub("ASD_Exome-","", rownames(asd_ex_sq_female))

# 添加标识
mark_positions <- asd_ex_sq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_ex_sq_female ), ncol = ncol(asd_ex_sq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_ex_sq_female)
colnames(mark_matrix) <- colnames(asd_ex_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_ex_subtype_female <- Heatmap(asd_ex_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                               width = ncol(asd_ex_ssd_female)*unit(6, "mm"), 
                               height = nrow(asd_ex_ssd_female)*unit(6, "mm"),column_title = "ASD Exome subtype",
                               cell_fun = text_fun)


##########################################################################################
#ASD SFARI subtype female
asd_sfari_ssd_female <- ssd[grep("^ASD_SFARI", ssd$CellType), 21:28]
rownames(asd_sfari_ssd_female) <- gsub("ASD_SFARI-","", rownames(asd_sfari_ssd_female))
asd_sfari_sq_female <- sq[grep("^ASD_SFARI", sq$CellType), 21:28]
rownames(asd_sfari_sq_female) <- gsub("ASD_SFARI-","", rownames(asd_sfari_sq_female))

# 添加标识
mark_positions <- asd_sfari_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_sfari_sq_female), ncol = ncol(asd_sfari_sq_female))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_sfari_sq_female)
colnames(mark_matrix) <- colnames(asd_sfari_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_sfari_subtype_female <- Heatmap(asd_sfari_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                  name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                  width = ncol(asd_sfari_ssd_female)*unit(6, "mm"), 
                                  height = nrow(asd_sfari_ssd_female)*unit(6, "mm"),column_title = "ASD SFARI subtype",
                                  left_annotation = left_annotation, cell_fun = text_fun)


##########################################################################################
#BD Exome subtype female
bd_ex_ssd_female <- ssd[grep("^BD_Exome", ssd$CellType), 21:28]
rownames(bd_ex_ssd_female) <- gsub("BD_Exome-","", rownames(bd_ex_ssd_female))
bd_ex_sq_female <- sq[grep("^BD_Exome", sq$CellType), 21:28]
rownames(bd_ex_sq_female) <- gsub("BD_Exome-","", rownames(bd_ex_sq_female))

# 添加标识
mark_positions <- bd_ex_sq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_ex_sq_female ), ncol = ncol(bd_ex_sq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_ex_sq_female)
colnames(mark_matrix) <- colnames(bd_ex_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
bd_ex_subtype_female <- Heatmap(bd_ex_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                              name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                              width = ncol(bd_ex_ssd_female)*unit(6, "mm"), 
                              height = nrow(bd_ex_ssd_female)*unit(6, "mm"),column_title = "BD Exome subtype",
                              cell_fun = text_fun)


##########################################################################################
#SCZ Exome subtype female
scz_ex_ssd_female <- ssd[grep("^SCZ_Exome", ssd$CellType), 21:28]
rownames(scz_ex_ssd_female) <- gsub("SCZ_Exome-","", rownames(scz_ex_ssd_female))
scz_ex_sq_female <- sq[grep("^SCZ_Exome", sq$CellType), 21:28]
rownames(scz_ex_sq_female) <- gsub("SCZ_Exome-","", rownames(scz_ex_sq_female))

# 添加标识
mark_positions <- scz_ex_sq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_ex_sq_female ), ncol = ncol(scz_ex_sq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_ex_sq_female)
colnames(mark_matrix) <- colnames(scz_ex_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
scz_ex_subtype_female <- Heatmap(scz_ex_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                               width = ncol(scz_ex_ssd_female)*unit(6, "mm"), 
                               height = nrow(scz_ex_ssd_female)*unit(6, "mm"),column_title = "SCZ Exome subtype",
                               cell_fun = text_fun)



#EWCE_All_Exome_Stages_Subtype_female
asd_sfari_subtype_female + asd_ex_subtype_female + bd_ex_subtype_female + scz_ex_subtype_female





####################################################################################################
#plot the complexheatmap across different developmental stages
tsd <- read.table("../EWCE/EWCE_Type_sd.txt", header = T)
tq <- read.table("../EWCE/EWCE_Type_q.txt", header = T)
rownames(tsd) <- tsd$CellType
rownames(tq) <- tq$CellType

##########################################################################################
#ASD Exome type
asd_ex_tsd <- tsd[grep("^ASD_Exome", tsd$CellType), 5:12]
rownames(asd_ex_tsd) <- gsub("ASD_Exome-","", rownames(asd_ex_tsd))
asd_ex_tq <- tq[grep("^ASD_Exome", tq$CellType), 5:12]
rownames(asd_ex_tq) <- gsub("ASD_Exome-","", rownames(asd_ex_tq))

# 添加标识
mark_positions <- asd_ex_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_ex_tq ), ncol = ncol(asd_ex_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_ex_tq)
colnames(mark_matrix) <- colnames(asd_ex_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_ex_type <- Heatmap(asd_ex_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(asd_ex_tsd)*unit(6, "mm"), 
                          height = nrow(asd_ex_tsd)*unit(6, "mm"),column_title = "ASD Exome type",
                          cell_fun = text_fun)


##########################################################################################
#asd SFARI type
asd_sfari_tsd <- tsd[grep("^ASD_SFARI", tsd$CellType), 5:12]
rownames(asd_sfari_tsd) <- gsub("ASD_SFARI-","", rownames(asd_sfari_tsd))
asd_sfari_tq <- tq[grep("^ASD_SFARI", tq$CellType), 5:12]
rownames(asd_sfari_tq) <- gsub("ASD_SFARI-","", rownames(asd_sfari_tq))

# 添加标识
mark_positions <- asd_sfari_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_sfari_tq ), ncol = ncol(asd_sfari_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_sfari_tq)
colnames(mark_matrix) <- colnames(asd_sfari_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_sfari_type <- Heatmap(asd_sfari_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                             name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                             width = ncol(asd_sfari_tsd)*unit(6, "mm"), 
                             height = nrow(asd_sfari_tsd)*unit(6, "mm"),column_title = "ASD SFARI type",
                             cell_fun = text_fun)


##########################################################################################
#BD Exome type
bd_ex_tsd <- tsd[grep("^BD_Exome", tsd$CellType), 5:12]
rownames(bd_ex_tsd) <- gsub("BD_Exome-","", rownames(bd_ex_tsd))
bd_ex_tq <- tq[grep("^BD_Exome", tq$CellType), 5:12]
rownames(bd_ex_tq) <- gsub("BD_Exome-","", rownames(bd_ex_tq))

# 添加标识
mark_positions <- bd_ex_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_ex_tq ), ncol = ncol(bd_ex_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_ex_tq)
colnames(mark_matrix) <- colnames(bd_ex_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
bd_ex_type <- Heatmap(bd_ex_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                         name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                         width = ncol(bd_ex_tsd)*unit(6, "mm"), 
                         height = nrow(bd_ex_tsd)*unit(6, "mm"),column_title = "BD Exome type",
                         cell_fun = text_fun)


##########################################################################################
#scz Exome type
scz_ex_tsd <- tsd[grep("^SCZ_Exome", tsd$CellType), 5:12]
rownames(scz_ex_tsd) <- gsub("SCZ_Exome-","", rownames(scz_ex_tsd))
scz_ex_tq <- tq[grep("^SCZ_Exome", tq$CellType), 5:12]
rownames(scz_ex_tq) <- gsub("SCZ_Exome-","", rownames(scz_ex_tq))

# 添加标识
mark_positions <- scz_ex_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_ex_tq ), ncol = ncol(scz_ex_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_ex_tq)
colnames(mark_matrix) <- colnames(scz_ex_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
scz_ex_type <- Heatmap(scz_ex_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(scz_ex_tsd)*unit(6, "mm"), 
                          height = nrow(scz_ex_tsd)*unit(6, "mm"),column_title = "SCZ Exome type",
                          cell_fun = text_fun)



#EWCE_All_Exome_Stages_Type
asd_sfari_type + asd_ex_type + bd_ex_type + scz_ex_type



##########################################################################################
#ASD Exome type Male
asd_ex_tsd_male <- tsd[grep("^ASD_Exome", tsd$CellType), 13:20]
rownames(asd_ex_tsd_male) <- gsub("ASD_Exome-","", rownames(asd_ex_tsd_male))
asd_ex_tq_male <- tq[grep("^ASD_Exome", tq$CellType), 13:20]
rownames(asd_ex_tq_male) <- gsub("ASD_Exome-","", rownames(asd_ex_tq_male))

# 添加标识
mark_positions <- asd_ex_tq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_ex_tq_male ), ncol = ncol(asd_ex_tq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_ex_tq_male)
colnames(mark_matrix) <- colnames(asd_ex_tq_male)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_ex_type_male <- Heatmap(asd_ex_tsd_male[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                               width = ncol(asd_ex_tsd_male)*unit(6, "mm"), 
                               height = nrow(asd_ex_tsd_male)*unit(6, "mm"),column_title = "ASD Exome type",
                               cell_fun = text_fun)


##########################################################################################
#asd SFARI type
asd_sfari_tsd_male  <- tsd[grep("^ASD_SFARI", tsd$CellType), 13:20]
rownames(asd_sfari_tsd_male ) <- gsub("ASD_SFARI-","", rownames(asd_sfari_tsd_male ))
asd_sfari_tq_male  <- tq[grep("^ASD_SFARI", tq$CellType), 13:20]
rownames(asd_sfari_tq_male ) <- gsub("ASD_SFARI-","", rownames(asd_sfari_tq_male ))

# 添加标识
mark_positions <- asd_sfari_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_sfari_tq ), ncol = ncol(asd_sfari_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_sfari_tq)
colnames(mark_matrix) <- colnames(asd_sfari_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_sfari_type_male <- Heatmap(asd_sfari_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                  name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                  width = ncol(asd_sfari_tsd)*unit(6, "mm"), 
                                  height = nrow(asd_sfari_tsd)*unit(6, "mm"),column_title = "ASD SFARI type",
                                  cell_fun = text_fun)


##########################################################################################
#BD Exome type
bd_ex_tsd_male <- tsd[grep("^BD_Exome", tsd$CellType), 13:20]
rownames(bd_ex_tsd_male) <- gsub("BD_Exome-","", rownames(bd_ex_tsd_male))
bd_ex_tq_male <- tq[grep("^BD_Exome", tq$CellType), 13:20]
rownames(bd_ex_tq_male) <- gsub("BD_Exome-","", rownames(bd_ex_tq_male))

# 添加标识
mark_positions <- bd_ex_tq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_ex_tq_male ), ncol = ncol(bd_ex_tq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_ex_tq_male)
colnames(mark_matrix) <- colnames(bd_ex_tq_male)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
bd_ex_type_male <- Heatmap(bd_ex_tsd_male[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                              name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                              width = ncol(bd_ex_tsd_male)*unit(6, "mm"), 
                              height = nrow(bd_ex_tsd_male)*unit(6, "mm"),column_title = "BD Exome type",
                              cell_fun = text_fun)


##########################################################################################
#scz Exome type
scz_ex_tsd_male <- tsd[grep("^SCZ_Exome", tsd$CellType), 13:20]
rownames(scz_ex_tsd_male) <- gsub("SCZ_Exome-","", rownames(scz_ex_tsd_male))
scz_ex_tq_male <- tq[grep("^SCZ_Exome", tq$CellType), 13:20]
rownames(scz_ex_tq_male) <- gsub("SCZ_Exome-","", rownames(scz_ex_tq_male))

# 添加标识
mark_positions <- scz_ex_tq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_ex_tq_male ), ncol = ncol(scz_ex_tq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_ex_tq_male)
colnames(mark_matrix) <- colnames(scz_ex_tq_male)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
scz_ex_type_male <- Heatmap(scz_ex_tsd_male[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                               width = ncol(scz_ex_tsd_male)*unit(6, "mm"), 
                               height = nrow(scz_ex_tsd_male)*unit(6, "mm"),column_title = "SCZ Exome type",
                               cell_fun = text_fun)



#EWCE_All_Exome_Stages_Type_Male
asd_sfari_type_male + asd_ex_type_male + bd_ex_type_male + scz_ex_type_male



##########################################################################################
#ASD Exome type female
asd_ex_tsd_female <- tsd[grep("^ASD_Exome", tsd$CellType), 21:28]
rownames(asd_ex_tsd_female) <- gsub("ASD_Exome-","", rownames(asd_ex_tsd_female))
asd_ex_tq_female <- tq[grep("^ASD_Exome", tq$CellType), 21:28]
rownames(asd_ex_tq_female) <- gsub("ASD_Exome-","", rownames(asd_ex_tq_female))

# 添加标识
mark_positions <- asd_ex_tq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_ex_tq_female ), ncol = ncol(asd_ex_tq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_ex_tq_female)
colnames(mark_matrix) <- colnames(asd_ex_tq_female)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_ex_type_female <- Heatmap(asd_ex_tsd_female[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                 name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                 width = ncol(asd_ex_tsd_female)*unit(6, "mm"), 
                                 height = nrow(asd_ex_tsd_female)*unit(6, "mm"),column_title = "ASD Exome type",
                                 cell_fun = text_fun)


##########################################################################################
#asd SFARI type
asd_sfari_tsd_female <- tsd[grep("^ASD_SFARI", tsd$CellType), 21:28]
rownames(asd_sfari_tsd_female) <- gsub("ASD_SFARI-","", rownames(asd_sfari_tsd_female))
asd_sfari_tq_female <- tq[grep("^ASD_SFARI", tq$CellType), 21:28]
rownames(asd_sfari_tq_female) <- gsub("ASD_SFARI-","", rownames(asd_sfari_tq_female))

# 添加标识
mark_positions <- asd_sfari_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_sfari_tq ), ncol = ncol(asd_sfari_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_sfari_tq)
colnames(mark_matrix) <- colnames(asd_sfari_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_sfari_type_female <- Heatmap(asd_sfari_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                    name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                    width = ncol(asd_sfari_tsd)*unit(6, "mm"), 
                                    height = nrow(asd_sfari_tsd)*unit(6, "mm"),column_title = "ASD SFARI type",
                                    cell_fun = text_fun)


##########################################################################################
#BD Exome type
bd_ex_tsd_female <- tsd[grep("^BD_Exome", tsd$CellType), 21:28]
rownames(bd_ex_tsd_female) <- gsub("BD_Exome-","", rownames(bd_ex_tsd_female))
bd_ex_tq_female <- tq[grep("^BD_Exome", tq$CellType), 21:28]
rownames(bd_ex_tq_female) <- gsub("BD_Exome-","", rownames(bd_ex_tq_female))

# 添加标识
mark_positions <- bd_ex_tq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_ex_tq_female ), ncol = ncol(bd_ex_tq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_ex_tq_female)
colnames(mark_matrix) <- colnames(bd_ex_tq_female)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
bd_ex_type_female <- Heatmap(bd_ex_tsd_female[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                width = ncol(bd_ex_tsd_female)*unit(6, "mm"), 
                                height = nrow(bd_ex_tsd_female)*unit(6, "mm"),column_title = "BD Exome type",
                                cell_fun = text_fun)


##########################################################################################
#scz Exome type
scz_ex_tsd_female <- tsd[grep("^SCZ_Exome", tsd$CellType), 21:28]
rownames(scz_ex_tsd_female) <- gsub("SCZ_Exome-","", rownames(scz_ex_tsd_female))
scz_ex_tq_female <- tq[grep("^SCZ_Exome", tq$CellType), 21:28]
rownames(scz_ex_tq_female) <- gsub("SCZ_Exome-","", rownames(scz_ex_tq_female))

# 添加标识
mark_positions <- scz_ex_tq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_ex_tq_female ), ncol = ncol(scz_ex_tq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_ex_tq_female)
colnames(mark_matrix) <- colnames(scz_ex_tq_female)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
scz_ex_type_female <- Heatmap(scz_ex_tsd_female[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                 name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                 width = ncol(scz_ex_tsd_female)*unit(6, "mm"), 
                                 height = nrow(scz_ex_tsd_female)*unit(6, "mm"),column_title = "SCZ Exome type",
                                 cell_fun = text_fun)



#EWCE_All_Exome_Stages_Type_female
asd_sfari_type_female + asd_ex_type_female + bd_ex_type_female + scz_ex_type_female




##########################################################################################
#ASD Diff subtype
asd_dif_ssd <- ssd[grep("^ASD_Diff", ssd$CellType), 5:12]
rownames(asd_dif_ssd) <- gsub("ASD_Diff-","", rownames(asd_dif_ssd))
asd_dif_sq <- sq[grep("^ASD_Diff", sq$CellType), 5:12]
rownames(asd_dif_sq) <- gsub("ASD_Diff-","", rownames(asd_dif_sq))

# 添加标识
mark_positions <- asd_dif_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_dif_sq ), ncol = ncol(asd_dif_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_dif_sq)
colnames(mark_matrix) <- colnames(asd_dif_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_dif_subtype <- Heatmap(asd_dif_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                           name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                           width = ncol(asd_dif_ssd)*unit(6, "mm"), 
                           height = nrow(asd_dif_ssd)*unit(6, "mm"),column_title = "ASD Diff subtype",
                           left_annotation = left_annotation, cell_fun = text_fun)


##########################################################################################
#BD Diff subtype
bd_dif_ssd <- ssd[grep("^BD_Diff", ssd$CellType), 5:12]
rownames(bd_dif_ssd) <- gsub("BD_Diff-","", rownames(bd_dif_ssd))
bd_dif_sq <- sq[grep("^BD_Diff", sq$CellType), 5:12]
rownames(bd_dif_sq) <- gsub("BD_Diff-","", rownames(bd_dif_sq))

# 添加标识
mark_positions <- bd_dif_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_dif_sq ), ncol = ncol(bd_dif_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_dif_sq)
colnames(mark_matrix) <- colnames(bd_dif_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
bd_dif_subtype <- Heatmap(bd_dif_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                          name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                          width = ncol(bd_dif_ssd)*unit(6, "mm"), 
                          height = nrow(bd_dif_ssd)*unit(6, "mm"),column_title = "BD Diff subtype",
                          cell_fun = text_fun)


##########################################################################################
#scz Diff subtype
scz_dif_ssd <- ssd[grep("^SCZ_Diff", ssd$CellType), 5:12]
rownames(scz_dif_ssd) <- gsub("SCZ_Diff-","", rownames(scz_dif_ssd))
scz_dif_sq <- sq[grep("^SCZ_Diff", sq$CellType), 5:12]
rownames(scz_dif_sq) <- gsub("SCZ_Diff-","", rownames(scz_dif_sq))

# 添加标识
mark_positions <- scz_dif_sq < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_dif_sq ), ncol = ncol(scz_dif_sq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_dif_sq)
colnames(mark_matrix) <- colnames(scz_dif_sq)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
scz_dif_subtype <- Heatmap(scz_dif_ssd[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                           name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                           width = ncol(scz_dif_ssd)*unit(6, "mm"), 
                           height = nrow(scz_dif_ssd)*unit(6, "mm"),column_title = "SCZ Diff subtype",
                           cell_fun = text_fun)



#EWCE_All_Diff_Stages_Subtype
asd_dif_subtype + bd_dif_subtype + scz_dif_subtype


##########################################################################################
#ASD Diff subtype Male
asd_dif_ssd_male <- ssd[grep("^ASD_Diff", ssd$CellType), 13:20]
rownames(asd_dif_ssd_male) <- gsub("ASD_Diff-","", rownames(asd_dif_ssd_male))
asd_dif_sq_male <- sq[grep("^ASD_Diff", sq$CellType), 13:20]
rownames(asd_dif_sq_male) <- gsub("ASD_Diff-","", rownames(asd_dif_sq_male))

# 添加标识
mark_positions <- asd_dif_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_dif_sq_male ), ncol = ncol(asd_dif_sq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_dif_sq_male)
colnames(mark_matrix) <- colnames(asd_dif_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_dif_subtype_male <- Heatmap(asd_dif_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                               width = ncol(asd_dif_ssd_male)*unit(6, "mm"), 
                               height = nrow(asd_dif_ssd_male)*unit(6, "mm"),column_title = "ASD Diff subtype",
                               left_annotation = left_annotation, cell_fun = text_fun)


##########################################################################################
#BD Diff subtype
bd_dif_ssd_male <- ssd[grep("^BD_Diff", ssd$CellType), 13:20]
rownames(bd_dif_ssd_male) <- gsub("BD_Diff-","", rownames(bd_dif_ssd_male))
bd_dif_sq_male <- sq[grep("^BD_Diff", sq$CellType), 13:20]
rownames(bd_dif_sq_male) <- gsub("BD_Diff-","", rownames(bd_dif_sq_male))

# 添加标识
mark_positions <- bd_dif_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_dif_sq_male ), ncol = ncol(bd_dif_sq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_dif_sq_male)
colnames(mark_matrix) <- colnames(bd_dif_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
bd_dif_subtype_male <- Heatmap(bd_dif_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                              name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                              width = ncol(bd_dif_ssd_male)*unit(6, "mm"), 
                              height = nrow(bd_dif_ssd_male)*unit(6, "mm"),column_title = "BD Diff subtype",
                              cell_fun = text_fun)


##########################################################################################
#scz Diff subtype
scz_dif_ssd_male <- ssd[grep("^SCZ_Diff", ssd$CellType), 13:20]
rownames(scz_dif_ssd_male) <- gsub("SCZ_Diff-","", rownames(scz_dif_ssd_male))
scz_dif_sq_male <- sq[grep("^SCZ_Diff", sq$CellType), 13:20]
rownames(scz_dif_sq_male) <- gsub("SCZ_Diff-","", rownames(scz_dif_sq_male))

# 添加标识
mark_positions <- scz_dif_sq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_dif_sq_male ), ncol = ncol(scz_dif_sq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_dif_sq_male)
colnames(mark_matrix) <- colnames(scz_dif_sq_male)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
scz_dif_subtype_male <- Heatmap(scz_dif_ssd_male[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                               name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                               width = ncol(scz_dif_ssd_male)*unit(6, "mm"), 
                               height = nrow(scz_dif_ssd_male)*unit(6, "mm"),column_title = "SCZ Diff subtype",
                               cell_fun = text_fun)



#EWCE_All_Diff_Stages_Subtype_Male
asd_dif_subtype_male + bd_dif_subtype_male + scz_dif_subtype_male



##########################################################################################
#ASD Diff subtype female
asd_dif_ssd_female <- ssd[grep("^ASD_Diff", ssd$CellType), 21:28]
rownames(asd_dif_ssd_female) <- gsub("ASD_Diff-","", rownames(asd_dif_ssd_female))
asd_dif_sq_female <- sq[grep("^ASD_Diff", sq$CellType), 21:28]
rownames(asd_dif_sq_female) <- gsub("ASD_Diff-","", rownames(asd_dif_sq_female))

# 添加标识
mark_positions <- asd_dif_sq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_dif_sq_female ), ncol = ncol(asd_dif_sq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_dif_sq_female)
colnames(mark_matrix) <- colnames(asd_dif_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
asd_dif_subtype_female <- Heatmap(asd_dif_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                 name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                 width = ncol(asd_dif_ssd_female)*unit(6, "mm"), 
                                 height = nrow(asd_dif_ssd_female)*unit(6, "mm"),column_title = "ASD Diff subtype",
                                 left_annotation = left_annotation, cell_fun = text_fun)



##########################################################################################
#BD Diff subtype female
bd_dif_ssd_female <- ssd[grep("^BD_Diff", ssd$CellType), 21:28]
rownames(bd_dif_ssd_female) <- gsub("BD_Diff-","", rownames(bd_dif_ssd_female))
bd_dif_sq_female <- sq[grep("^BD_Diff", sq$CellType), 21:28]
rownames(bd_dif_sq_female) <- gsub("BD_Diff-","", rownames(bd_dif_sq_female))

# 添加标识
mark_positions <- bd_dif_sq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_dif_sq_female ), ncol = ncol(bd_dif_sq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_dif_sq_female)
colnames(mark_matrix) <- colnames(bd_dif_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
bd_dif_subtype_female <- Heatmap(bd_dif_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                width = ncol(bd_dif_ssd_female)*unit(6, "mm"), 
                                height = nrow(bd_dif_ssd_female)*unit(6, "mm"),column_title = "BD Diff subtype",
                                cell_fun = text_fun)


##########################################################################################
#SCZ Diff subtype female
scz_dif_ssd_female <- ssd[grep("^SCZ_Diff", ssd$CellType), 21:28]
rownames(scz_dif_ssd_female) <- gsub("SCZ_Diff-","", rownames(scz_dif_ssd_female))
scz_dif_sq_female <- sq[grep("^SCZ_Diff", sq$CellType), 21:28]
rownames(scz_dif_sq_female) <- gsub("SCZ_Diff-","", rownames(scz_dif_sq_female))

# 添加标识
mark_positions <- scz_dif_sq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_dif_sq_female ), ncol = ncol(scz_dif_sq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_dif_sq_female)
colnames(mark_matrix) <- colnames(scz_dif_sq_female)
mark_matrix <- mark_matrix[row_subtype_name,]

# 绘制热图
scz_dif_subtype_female <- Heatmap(scz_dif_ssd_female[row_subtype_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                                 name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                                 width = ncol(scz_dif_ssd_female)*unit(6, "mm"), 
                                 height = nrow(scz_dif_ssd_female)*unit(6, "mm"),column_title = "SCZ Diff subtype",
                                 cell_fun = text_fun)



#EWCE_All_Diff_Stages_Subtype_female
asd_dif_subtype_female + bd_dif_subtype_female + scz_dif_subtype_female


####################################################################################################
#plot the complexheatmap across different developmental stages
tsd <- read.table("../EWCE/EWCE_Type_sd.txt", header = T)
tq <- read.table("../EWCE/EWCE_Type_q.txt", header = T)
rownames(tsd) <- tsd$CellType
rownames(tq) <- tq$CellType

##########################################################################################
#ASD Diff type
asd_dif_tsd <- tsd[grep("^ASD_Diff", tsd$CellType), 5:12]
rownames(asd_dif_tsd) <- gsub("ASD_Diff-","", rownames(asd_dif_tsd))
asd_dif_tq <- tq[grep("^ASD_Diff", tq$CellType), 5:12]
rownames(asd_dif_tq) <- gsub("ASD_Diff-","", rownames(asd_dif_tq))

# 添加标识
mark_positions <- asd_dif_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_dif_tq ), ncol = ncol(asd_dif_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_dif_tq)
colnames(mark_matrix) <- colnames(asd_dif_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_dif_type <- Heatmap(asd_dif_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                       name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                       width = ncol(asd_dif_tsd)*unit(6, "mm"), 
                       height = nrow(asd_dif_tsd)*unit(6, "mm"),column_title = "ASD Diff type",
                       cell_fun = text_fun)


##########################################################################################
#BD Diff type
bd_dif_tsd <- tsd[grep("^BD_Diff", tsd$CellType), 5:12]
rownames(bd_dif_tsd) <- gsub("BD_Diff-","", rownames(bd_dif_tsd))
bd_dif_tq <- tq[grep("^BD_Diff", tq$CellType), 5:12]
rownames(bd_dif_tq) <- gsub("BD_Diff-","", rownames(bd_dif_tq))

# 添加标识
mark_positions <- bd_dif_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_dif_tq ), ncol = ncol(bd_dif_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_dif_tq)
colnames(mark_matrix) <- colnames(bd_dif_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
bd_dif_type <- Heatmap(bd_dif_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                      name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                      width = ncol(bd_dif_tsd)*unit(6, "mm"), 
                      height = nrow(bd_dif_tsd)*unit(6, "mm"),column_title = "BD Diff type",
                      cell_fun = text_fun)


##########################################################################################
#scz Diff type
scz_dif_tsd <- tsd[grep("^SCZ_Diff", tsd$CellType), 5:12]
rownames(scz_dif_tsd) <- gsub("SCZ_Diff-","", rownames(scz_dif_tsd))
scz_dif_tq <- tq[grep("^SCZ_Diff", tq$CellType), 5:12]
rownames(scz_dif_tq) <- gsub("SCZ_Diff-","", rownames(scz_dif_tq))

# 添加标识
mark_positions <- scz_dif_tq < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_dif_tq ), ncol = ncol(scz_dif_tq ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_dif_tq)
colnames(mark_matrix) <- colnames(scz_dif_tq)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
scz_dif_type <- Heatmap(scz_dif_tsd[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                       name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                       width = ncol(scz_dif_tsd)*unit(6, "mm"), 
                       height = nrow(scz_dif_tsd)*unit(6, "mm"),column_title = "SCZ Diff type",
                       cell_fun = text_fun)



#EWCE_All_Diff_Stages_Type
asd_dif_type + bd_dif_type + scz_dif_type



##########################################################################################
#ASD Diff type Male
asd_dif_tsd_male <- tsd[grep("^ASD_Diff", tsd$CellType), 13:20]
rownames(asd_dif_tsd_male) <- gsub("ASD_Diff-","", rownames(asd_dif_tsd_male))
asd_dif_tq_male <- tq[grep("^ASD_Diff", tq$CellType), 13:20]
rownames(asd_dif_tq_male) <- gsub("ASD_Diff-","", rownames(asd_dif_tq_male))

# 添加标识
mark_positions <- asd_dif_tq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_dif_tq_male ), ncol = ncol(asd_dif_tq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_dif_tq_male)
colnames(mark_matrix) <- colnames(asd_dif_tq_male)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_dif_type_male <- Heatmap(asd_dif_tsd_male[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                            name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                            width = ncol(asd_dif_tsd_male)*unit(6, "mm"), 
                            height = nrow(asd_dif_tsd_male)*unit(6, "mm"),column_title = "ASD Diff type",
                            cell_fun = text_fun)


##########################################################################################
#BD Diff type Male
bd_dif_tsd_male <- tsd[grep("^BD_Diff", tsd$CellType), 13:20]
rownames(bd_dif_tsd_male) <- gsub("BD_Diff-","", rownames(bd_dif_tsd_male))
bd_dif_tq_male <- tq[grep("^BD_Diff", tq$CellType), 13:20]
rownames(bd_dif_tq_male) <- gsub("BD_Diff-","", rownames(bd_dif_tq_male))

# 添加标识
mark_positions <- bd_dif_tq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_dif_tq_male ), ncol = ncol(bd_dif_tq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_dif_tq_male)
colnames(mark_matrix) <- colnames(bd_dif_tq_male)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
bd_dif_type_male <- Heatmap(bd_dif_tsd_male[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                           name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                           width = ncol(bd_dif_tsd_male)*unit(6, "mm"), 
                           height = nrow(bd_dif_tsd_male)*unit(6, "mm"),column_title = "BD Diff type",
                           cell_fun = text_fun)


##########################################################################################
#scz Diff type Male
scz_dif_tsd_male <- tsd[grep("^SCZ_Diff", tsd$CellType), 13:20]
rownames(scz_dif_tsd_male) <- gsub("SCZ_Diff-","", rownames(scz_dif_tsd_male))
scz_dif_tq_male <- tq[grep("^SCZ_Diff", tq$CellType), 13:20]
rownames(scz_dif_tq_male) <- gsub("SCZ_Diff-","", rownames(scz_dif_tq_male))

# 添加标识
mark_positions <- scz_dif_tq_male < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_dif_tq_male ), ncol = ncol(scz_dif_tq_male ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_dif_tq_male)
colnames(mark_matrix) <- colnames(scz_dif_tq_male)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
scz_dif_type_male <- Heatmap(scz_dif_tsd_male[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                            name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                            width = ncol(scz_dif_tsd_male)*unit(6, "mm"), 
                            height = nrow(scz_dif_tsd_male)*unit(6, "mm"),column_title = "SCZ Diff type",
                            cell_fun = text_fun)



#EWCE_All_Diff_Stages_Type_Male
asd_dif_type_male + bd_dif_type_male + scz_dif_type_male



##########################################################################################
#ASD Diff type female
asd_dif_tsd_female <- tsd[grep("^ASD_Diff", tsd$CellType), 21:28]
rownames(asd_dif_tsd_female) <- gsub("ASD_Diff-","", rownames(asd_dif_tsd_female))
asd_dif_tq_female <- tq[grep("^ASD_Diff", tq$CellType), 21:28]
rownames(asd_dif_tq_female) <- gsub("ASD_Diff-","", rownames(asd_dif_tq_female))

# 添加标识
mark_positions <- asd_dif_tq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(asd_dif_tq_female ), ncol = ncol(asd_dif_tq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(asd_dif_tq_female)
colnames(mark_matrix) <- colnames(asd_dif_tq_female)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
asd_dif_type_female <- Heatmap(asd_dif_tsd_female[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                              name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                              width = ncol(asd_dif_tsd_female)*unit(6, "mm"), 
                              height = nrow(asd_dif_tsd_female)*unit(6, "mm"),column_title = "ASD Diff type",
                              cell_fun = text_fun)


##########################################################################################
#BD Diff type
bd_dif_tsd_female <- tsd[grep("^BD_Diff", tsd$CellType), 21:28]
rownames(bd_dif_tsd_female) <- gsub("BD_Diff-","", rownames(bd_dif_tsd_female))
bd_dif_tq_female <- tq[grep("^BD_Diff", tq$CellType), 21:28]
rownames(bd_dif_tq_female) <- gsub("BD_Diff-","", rownames(bd_dif_tq_female))

# 添加标识
mark_positions <- bd_dif_tq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(bd_dif_tq_female ), ncol = ncol(bd_dif_tq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(bd_dif_tq_female)
colnames(mark_matrix) <- colnames(bd_dif_tq_female)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
bd_dif_type_female <- Heatmap(bd_dif_tsd_female[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                             name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                             width = ncol(bd_dif_tsd_female)*unit(6, "mm"), 
                             height = nrow(bd_dif_tsd_female)*unit(6, "mm"),column_title = "BD Diff type",
                             cell_fun = text_fun)


##########################################################################################
#scz Diff type
scz_dif_tsd_female <- tsd[grep("^SCZ_Diff", tsd$CellType), 21:28]
rownames(scz_dif_tsd_female) <- gsub("SCZ_Diff-","", rownames(scz_dif_tsd_female))
scz_dif_tq_female <- tq[grep("^SCZ_Diff", tq$CellType), 21:28]
rownames(scz_dif_tq_female) <- gsub("SCZ_Diff-","", rownames(scz_dif_tq_female))

# 添加标识
mark_positions <- scz_dif_tq_female < 0.05
mark_matrix <- matrix("", nrow = nrow(scz_dif_tq_female ), ncol = ncol(scz_dif_tq_female ))
mark_matrix[mark_positions] <- "*"
rownames(mark_matrix) <- rownames(scz_dif_tq_female)
colnames(mark_matrix) <- colnames(scz_dif_tq_female)
mark_matrix <- mark_matrix[row_type_name,]

# 绘制热图
scz_dif_type_female <- Heatmap(scz_dif_tsd_female[row_type_name,], cluster_rows = FALSE,cluster_columns = FALSE,
                              name = "s.d. from mean",col = colorRamp2(c(0, 2, 4, 6, 8), c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")),
                              width = ncol(scz_dif_tsd_female)*unit(6, "mm"), 
                              height = nrow(scz_dif_tsd_female)*unit(6, "mm"),column_title = "SCZ Diff type",
                              cell_fun = text_fun)



#EWCE_All_Diff_Stages_Type_female
asd_dif_type_female + bd_dif_type_female + scz_dif_type_female




##########################################################################################
#ASD SFARI sex comparsion barplot
asd_sfari_ssd_male$CellType <- rownames(asd_sfari_ssd_male)
asd_sfari_ssd_male <- as_tibble(asd_sfari_ssd_male)
asd_sfari_ssd_male_long <- asd_sfari_ssd_male %>% 
  gather(key = "group", value = "sd", trimester2ndMale:AdultMale)
asd_sfari_ssd_male_long$group <- gsub("Male", "", asd_sfari_ssd_male_long$group)

asd_sfari_sq_male$CellType <- rownames(asd_sfari_sq_male)
asd_sfari_sq_male <- as_tibble(asd_sfari_sq_male)
asd_sfari_sq_male_long <- asd_sfari_sq_male %>% 
  gather(key = "group", value = "q", trimester2ndMale:AdultMale)
asd_sfari_sq_male_long$group <- gsub("Male", "", asd_sfari_sq_male_long$group)

# merge sd and qvalue
mmerged_df <- merge(asd_sfari_ssd_male_long, asd_sfari_sq_male_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
mmerged_df <- mmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#377eb8", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
mmerged_df$sex <- "Male"


asd_sfari_ssd_female$CellType <- rownames(asd_sfari_ssd_female)
asd_sfari_ssd_female <- as_tibble(asd_sfari_ssd_female)
asd_sfari_ssd_female_long <- asd_sfari_ssd_female %>% 
  gather(key = "group", value = "sd", trimester2ndFemale:AdultFemale)
asd_sfari_ssd_female_long$group <- gsub("Female", "", asd_sfari_ssd_female_long$group)

asd_sfari_sq_female$CellType <- rownames(asd_sfari_sq_female)
asd_sfari_sq_female <- as_tibble(asd_sfari_sq_female)
asd_sfari_sq_female_long <- asd_sfari_sq_female %>% 
  gather(key = "group", value = "q", trimester2ndFemale:AdultFemale)
asd_sfari_sq_female_long$group <- gsub("Female", "", asd_sfari_sq_female_long$group)

# merge sd and qvalue
fmerged_df <- merge(asd_sfari_ssd_female_long, asd_sfari_sq_female_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
fmerged_df <- fmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#DB4D47", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
fmerged_df$sex <- "Female"

merged_df <- rbind(mmerged_df,fmerged_df)

merged_df$CellType <- factor(merged_df$CellType,levels = rev(c("EX_inter", "EX_L2_3", "EX_L4", "EX_L5", "EX_L6", "EX_L5_6_IT", "EX_SP", "EX_Progenitors", "In_INT", "In_SST", "In_RELN",
                                                           "In_SV2C", "In_VIP", "In_CCK", "In_SST_RELN", "In_CALB2", "In_PV", "In_NOS", "In_PV_MP", "In_Progenitors", "Micro",
                                                           "Oligos", "OPC", "Fibrous_astrocytes", "Protoplasmic_astrocytes", "Glial_progenitors", "Peri", "Vasc")))
merged_df$group <- factor(merged_df$group, levels = c("trimester2nd", "trimester3rd", "years0_1", "years1_2", 
                                                      "years2_4", "years4_10", "years10_20", "Adult"))

# 排除NA值对绘图的影响
merged_df$sd[is.na(merged_df$sd)] <- 0
merged_df$sd_abs[is.na(merged_df$sd_abs)] <- 0
merged_df$Significance <- ifelse(merged_df$q < 0.05 & merged_df$sd != 0, "*", "")
merged_df$Color <- ifelse(merged_df$sex == "Male", "#377eb8", "#DB4D47")

# plot the barplot
ggplot(merged_df, aes(y = CellType, x = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Significance, x = sd_abs + 0.2), position = position_dodge(width = 0.8), vjust = 0) +
  scale_fill_identity() +
  facet_wrap(~group, scales = "fixed", nrow = 1) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  )+
  labs(x = "", y = "s.d. from mean") +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))

ggsave(
  filename = "EWCE_ASD_Subtype_Sex.pdf",
  width = 10.55,
  height = 6.75, 
  units = "in",
  device = "pdf"
)




##########################################################################################
#ASD SFARI sex comparsion barplot
asd_sfari_tsd_male$CellType <- rownames(asd_sfari_tsd_male)
asd_sfari_tsd_male <- as_tibble(asd_sfari_tsd_male)
asd_sfari_tsd_male_long <- asd_sfari_tsd_male %>% 
  gather(key = "group", value = "sd", trimester2ndMale:AdultMale)
asd_sfari_tsd_male_long$group <- gsub("Male", "", asd_sfari_tsd_male_long$group)

asd_sfari_tq_male$CellType <- rownames(asd_sfari_tq_male)
asd_sfari_tq_male <- as_tibble(asd_sfari_tq_male)
asd_sfari_tq_male_long <- asd_sfari_tq_male %>% 
  gather(key = "group", value = "q", trimester2ndMale:AdultMale)
asd_sfari_tq_male_long$group <- gsub("Male", "", asd_sfari_tq_male_long$group)

# merge sd and qvalue
mmerged_df <- merge(asd_sfari_tsd_male_long, asd_sfari_tq_male_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
mmerged_df <- mmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#377eb8", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
mmerged_df$sex <- "Male"


asd_sfari_tsd_female$CellType <- rownames(asd_sfari_tsd_female)
asd_sfari_tsd_female <- as_tibble(asd_sfari_tsd_female)
asd_sfari_tsd_female_long <- asd_sfari_tsd_female %>% 
  gather(key = "group", value = "sd", trimester2ndFemale:AdultFemale)
asd_sfari_tsd_female_long$group <- gsub("Female", "", asd_sfari_tsd_female_long$group)

asd_sfari_tq_female$CellType <- rownames(asd_sfari_tq_female)
asd_sfari_tq_female <- as_tibble(asd_sfari_tq_female)
asd_sfari_tq_female_long <- asd_sfari_tq_female %>% 
  gather(key = "group", value = "q", trimester2ndFemale:AdultFemale)
asd_sfari_tq_female_long$group <- gsub("Female", "", asd_sfari_tq_female_long$group)

# merge sd and qvalue
fmerged_df <- merge(asd_sfari_tsd_female_long, asd_sfari_tq_female_long, by = c("CellType", "group"))

# add "*" and transform the sd to abs(sd)
fmerged_df <- fmerged_df %>%
  mutate(Significance = ifelse(q < 0.05, "*", ""),
         Color = ifelse(sd > 0, "#DB4D47", "white"),
         sd_abs = ifelse(sd > 0, sd, 0))
fmerged_df$sex <- "Female"

merged_df <- rbind(mmerged_df,fmerged_df)

merged_df$CellType <- factor(merged_df$CellType,levels = c("VASC","GLIALPROG","AST","OPC","OL","MG","IN","ExNeu"))
merged_df$group <- factor(merged_df$group, levels = c("trimester2nd", "trimester3rd", "years0_1", "years1_2", 
                                                      "years2_4", "years4_10", "years10_20", "Adult"))

# 排除NA值对绘图的影响
merged_df$sd[is.na(merged_df$sd)] <- 0
merged_df$sd_abs[is.na(merged_df$sd_abs)] <- 0
merged_df$Significance <- ifelse(merged_df$q < 0.05 & merged_df$sd != 0, "*", "")
merged_df$Color <- ifelse(merged_df$sex == "Male", "#377eb8", "#DB4D47")

# plot the barplot
ggplot(merged_df, aes(y = CellType, x = sd_abs, fill = Color)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Significance, x = sd_abs + 0.2), position = position_dodge(width = 0.8), vjust = 0) +
  scale_fill_identity() +
  facet_wrap(~group, scales = "fixed", nrow = 1) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold")
  )+
  labs(x = "", y = "s.d. from mean") +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))

ggsave(
  filename = "EWCE_ASD_Subtype_Sex.pdf",
  width = 10.55,
  height = 6.75, 
  units = "in",
  device = "pdf"
)