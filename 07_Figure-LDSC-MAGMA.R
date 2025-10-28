library(ggplot2)
library(reshape2)
library(tidyverse)
setwd("D:/CellType/downsample500")

####################################################################################
#Type
ldsc <- read.table("LDSC_Type.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("MAGMA_Type.txt", header = TRUE, row.names = 1, sep = "\t")

mean_log_p <- (-log10(ldsc) -log10(magma))/2

mean_log_p_melted <- melt(as.matrix(mean_log_p))

#Bonferroni correction, 10 cell types
#P = 0.05/10 = 0.005
threshold <- 0.005

ldsc_melted <- melt(as.matrix(ldsc))
magma_melted <- melt(as.matrix(magma))
combined <- merge(mean_log_p_melted, ldsc_melted, by = c("Var1", "Var2"))
combined <- merge(combined, magma_melted, by = c("Var1", "Var2"))
colnames(combined) <- c("CellType", "Stage", "Mean_Log_P", "LDSC_P", "MAGMA_P")
combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))

combined$Color <- with(combined, ifelse(Significance == "both", "#0291BA",
                                        ifelse(Significance == "ldsc", "#F4931A",
                                               ifelse(Significance == "magma", "#BD56A2", "#C8AF81"))))

all_type <- combined %>% filter(Stage == "Mean")
all_type <- separate(all_type, col = CellType, into = c("Dis", "Type"), sep = "-")
all_type$Dis <- gsub("20.*", "", all_type$Dis)
all_type$Type <- gsub("\\.", "-", all_type$Type)
all_type$Type <- factor(all_type$Type,levels=c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC"))


ggplot(all_type, aes(x = Type, y = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" 
  ) +
  facet_wrap(~ Dis, scales = "fixed", nrow = 2, strip.position = "right") +
  labs(x = "Cell Types", y = "Mean(-log10(P))") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#BD56A2", "#C8AF81")),
                             labels = c("both", "magma", "none"))) +
  scale_y_continuous(limits = c(0,6), breaks = c(0,2,4,6))

ggsave(
  filename = "LDSC-MAGMA_All_type.pdf", 
  width = 4.83,
  height = 4.49, 
  units = "in",
  device = "pdf"
)


#Bonferroni correction, 10 major cell types, 10 stages
#P = 0.05/(10*10) = 0.0005
threshold <- 0.0005
combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))
combined_sz <- combined %>% filter(str_detect(CellType, pattern = "^SCZ"), str_detect(Stage, pattern = "Mean_"))
combined_sz$CellType <- gsub("SCZ2022PGC3-", "", combined_sz$CellType)
combined_sz$Stage <- gsub("Mean_", "", combined_sz$Stage)
combined_bd <- combined %>% filter(str_detect(CellType, pattern = "^BD"), str_detect(Stage, pattern = "Mean_"))
combined_bd$CellType <- gsub("BD2025OConnel-", "", combined_bd$CellType)
combined_bd$Stage <- gsub("Mean_", "", combined_bd$Stage)

combined_sz$Stage <- factor(combined_sz$Stage, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))
combined_bd$Stage <- factor(combined_bd$Stage, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))
combined_sz$CellType <- gsub("\\.", "-", combined_sz$CellType)
combined_sz$CellType <- factor(combined_sz$CellType, levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))
combined_bd$CellType <- gsub("\\.", "-", combined_bd$CellType)
combined_bd$CellType <- factor(combined_bd$CellType, levels=rev(c("Neu-Progenitor","ExNeu","InNeu","Glia-Progenitor","Astro","OPC","Oligo","Micro", "Endo", "VLMC")))

ggplot(na.omit(combined_sz), aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" 
  ) +
  facet_wrap(~ Stage, scales = "fixed", nrow = 1) +
  labs(y = "Cell Types", x = "Mean(-log10(P))") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))
ggsave(
  filename = "LDSC-MAGMA-SCZ-stage-type.pdf", 
  width = 9.75,
  height = 3.33, 
  units = "in",
  device = "pdf"
)

ggplot(na.omit(combined_bd), aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~ Stage, scales = "fixed", nrow = 1) +
  labs(y = "Cell Types", x = "Mean(-log10(P))") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A",  "#C8AF81")),
                             labels = c("both", "ldsc", "none"))) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))
ggsave(
  filename = "LDSC-MAGMA-BD-stage-type.pdf", 
  width = 9.75,
  height = 3.33, 
  units = "in",
  device = "pdf"
)




#Subtype
ldsc <- read.table("LDSC_Subtype.txt", header = TRUE, row.names = 1, sep = "\t")
magma <- read.table("MAGMA_Subtype.txt", header = TRUE, row.names = 1, sep = "\t")

mean_log_p <- (-log10(ldsc) -log10(magma))/2

mean_log_p_melted <- melt(as.matrix(mean_log_p))

#Bonferroni correction, 28 Subtypes
#P = 0.05/28 = 0.0018
threshold <- 0.0018

ldsc_melted <- melt(as.matrix(ldsc))
magma_melted <- melt(as.matrix(magma))
combined <- merge(mean_log_p_melted, ldsc_melted, by = c("Var1", "Var2"))
combined <- merge(combined, magma_melted, by = c("Var1", "Var2"))
colnames(combined) <- c("CellType", "Stage", "Mean_Log_P", "LDSC_P", "MAGMA_P")
combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))
combined$Color <- with(combined, ifelse(Significance == "both", "#0291BA",
                                        ifelse(Significance == "ldsc", "#F4931A",
                                               ifelse(Significance == "magma", "#BD56A2", "#C8AF81"))))

all_type <- combined %>% filter(Stage == "Mean")
all_type <- separate(all_type, col = CellType, into = c("Dis", "Type"), sep = "-")
all_type$Dis <- gsub("20.*", "", all_type$Dis)
all_type$Type <- gsub("\\.", "-", all_type$Type)
all_type$Type <- factor(all_type$Type,levels=c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                               "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                               "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC"))

ggplot(all_type, aes(x = Type, y = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold"),
    #legend.position = "bottom", 
    #legend.direction = "horizontal", 
    legend.position = "none" 
  ) +
  facet_wrap(~ Dis, scales = "fixed", nrow = 2, strip.position = "right") +
  labs(x = "Cell Types", y = "Mean(-log10(P))") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_y_continuous(limits = c(0,6), breaks = c(0,2,4,6))

ggsave(
  filename = "LDSC-MAGMA_All_subtype.pdf", 
  width = 8.81,
  height = 4.49, 
  units = "in",
  device = "pdf"
)

#Bonferroni correction, 28 Subtypes, 10 stages
#P = 0.05/(28*10) = 0.00018
threshold <- 0.00018
combined$Significance <- with(combined, ifelse(LDSC_P < threshold & MAGMA_P < threshold, "both",
                                               ifelse(LDSC_P < threshold, "ldsc",
                                                      ifelse(MAGMA_P < threshold, "magma", "none"))))

combined_sz <- combined %>% filter(str_detect(CellType, pattern = "^SCZ"), str_detect(Stage, pattern = "Mean_"))
combined_sz$CellType <- gsub("SCZ2022PGC3-", "", combined_sz$CellType)
combined_sz$Stage <- gsub("Mean_", "", combined_sz$Stage)
combined_bd <- combined %>% filter(str_detect(CellType, pattern = "^BD"), str_detect(Stage, pattern = "Mean_"))
combined_bd$CellType <- gsub("BD2025OConnel-", "", combined_bd$CellType)
combined_bd$Stage <- gsub("Mean_", "", combined_bd$Stage)

combined_sz$Stage <- factor(combined_sz$Stage, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))
combined_bd$Stage <- factor(combined_bd$Stage, levels = c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"))
combined_sz$CellType <- gsub("\\.", "-", combined_sz$CellType)
combined_sz$CellType <- factor(combined_sz$CellType, levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                                  "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                                  "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))
combined_bd$CellType <- gsub("\\.", "-", combined_bd$CellType)
combined_bd$CellType <- factor(combined_bd$CellType, levels=rev(c("Ex-Progenitor", "Ex-Inter", "Ex-L23IT", "Ex-L4IT", "Ex-L5IT", "Ex-L6IT", "Ex-L6IT-Car3", "Ex-L56NP", "Ex-L5ET", "Ex-L6CT", "Ex-L6b",
                                                                  "In-Progenitor", "In-VIP", "In-SNCG", "In-LAMP5", "In-PAX6", "In-LAMP5-LHX6", "In-PVALB", "In-SST", "In-SST-CHODL", "In-Chandelier",
                                                                  "Glia-Progenitor", "Astro", "OPC", "Oligo", "Micro", "Endo", "VLMC")))

ggplot(na.omit(combined_sz), aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~ Stage, scales = "fixed", nrow = 1) +
  labs(y = "Cell Types", x = "Mean(-log10(P))") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2", "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))
ggsave(
  filename = "LDSC-MAGMA-SCZ-stage-subtype.pdf", 
  width = 9.75,
  height = 5.66, 
  units = "in",
  device = "pdf"
)

ggplot(na.omit(combined_bd), aes(y = CellType, x = Mean_Log_P, fill = Significance)) +
  geom_bar(stat = "identity",  position = position_dodge(), width = 0.75) +
  scale_fill_manual(values = c("both" = "#0291BA", "ldsc" = "#F4931A", "magma" = "#BD56A2", "none" = "#C8AF81")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#1C9E78", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" 
  ) +
  facet_wrap(~ Stage, scales = "fixed", nrow = 1) +
  labs(y = "Cell Types", x = "Mean(-log10(P))") +
  guides(fill = guide_legend(title = "Significance",
                             override.aes = list(fill = c("#0291BA", "#F4931A", "#BD56A2",  "#C8AF81")),
                             labels = c("both", "ldsc", "magma", "none"))) +
  scale_x_continuous(limits = c(0,9), breaks = c(0,3,6,9))
ggsave(
  filename = "LDSC-MAGMA-BD-stage-subtype.pdf", 
  width = 9.75,
  height = 5.66, 
  units = "in",
  device = "pdf"
)