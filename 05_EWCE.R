#2026-04-10
#wuyong0103@126.com
#Usage:
#Rscript EWCE_Stage.R

library(EWCE)
library(Seurat)
library(tidyverse)
library(data.table)

reps <- 100000
sz <- read_tsv("/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/SCHEMA_Exome.txt")
sz$bh <- p.adjust(sz$p, method = "BH")
#sz_genes <- sz$symbol[sz$bh<0.1]
sz_genes_05 <- sz$symbol[sz$p < 0.05]
sz_genes_01 <- sz$symbol[sz$p < 0.01]

bd <- read_tsv("/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/BipEx_Exome.txt")
bd$bh <- p.adjust(bd$p, method = "BH")
#bd_genes <- bd$symbol[bd$bh<0.1]
bd_genes_05 <- sz$symbol[bd$p < 0.05]
bd_genes_01 <- bd$symbol[bd$p < 0.01]

res_All_lineage_05 <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(res_All_lineage_05) <- c("CellType", "p", "fold_change", "sd_from_mean", "q", "Trait")
res_All_lineage_01 <- res_All_lineage_05
res_All_sublineage_05 <- res_All_lineage_05
res_All_sublineage_01 <- res_All_lineage_05

print("[INFO] Lineage level!")
data <- readRDS("/home/lilab/wuyong/project/scRNA/MP-R1/downsampleMin/downsampleMin_lineage.rds")
bg_genes <- Features(data)
data <- JoinLayers(data)
count <- GetAssayData(data, assay = "RNA", layer = "counts")
meta <- data@meta.data
exp <- EWCE::sct_normalize(count)
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = exp, input_species = "human", output_species = "human", level2annot = meta$lineage, no_cores=40)
annotLevels <- list(level1class=meta$lineage, level2class=meta$sub_lineage)
ctd <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = "All_lineage", savePath="/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result", no_cores=40)

load(ctd)
sz_lineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_05, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
sz_lineage_05 <- sz_lineage_05$results[order(sz_lineage_05$results$p),c(1,3:6)]
sz_lineage_05$Trait <- "SCZ"
res_All_lineage_05 <- rbind(res_All_lineage_05, sz_lineage_05)

sz_lineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_01, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
sz_lineage_01 <- sz_lineage_01$results[order(sz_lineage_01$results$p),c(1,3:6)]
sz_lineage_01$Trait <- "SCZ"
res_All_lineage_01 <- rbind(res_All_lineage_01, sz_lineage_01)

bd_lineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_05, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
bd_lineage_05 <- bd_lineage_05$results[order(bd_lineage_05$results$p),c(1,3:6)]
bd_lineage_05$Trait <- "BD"
res_All_lineage_05 <- rbind(res_All_lineage_05, bd_lineage_05)

bd_lineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_01, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
bd_lineage_01 <- bd_lineage_01$results[order(bd_lineage_01$results$p),c(1,3:6)]
bd_lineage_01$Trait <- "BD"
res_All_lineage_01 <- rbind(res_All_lineage_01, bd_lineage_01)


#sub lineage level
print("[INFO] Sublineage level!")
data_sub <- readRDS(paste0("/home/lilab/wuyong/project/scRNA/MP-R1/downsampleMin/downsampleMin_sublineage.rds"))
data_sub <- JoinLayers(data_sub)
count_sub <- GetAssayData(data_sub, assay = "RNA", layer = "counts")
meta_sub <- data_sub@meta.data
exp_sub <- EWCE::sct_normalize(count_sub)
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = exp_sub, input_species = "human", output_species = "human", level2annot = meta_sub$sub_lineage, no_cores=40)
annotLevels <- list(level1class=meta_sub$lineage, level2class=meta_sub$sub_lineage)
ctd_sub <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = "All_sublineage", savePath="/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result", no_cores=40)

load(ctd_sub)
sz_sublineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_05, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
sz_sublineage_05 <- sz_sublineage_05$results[order(sz_sublineage_05$results$p),c(1,3:6)]
sz_sublineage_05$Trait <- "SCZ"
res_All_sublineage_05 <- rbind(res_All_sublineage_05, sz_sublineage_05)

sz_sublineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_01, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
sz_sublineage_01 <- sz_sublineage_01$results[order(sz_sublineage_01$results$p),c(1,3:6)]
sz_sublineage_01$Trait <- "SCZ"
res_All_sublineage_01 <- rbind(res_All_sublineage_01, sz_sublineage_01)

bd_sublineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_05, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
bd_sublineage_05 <- bd_sublineage_05$results[order(bd_sublineage_05$results$p),c(1,3:6)]
bd_sublineage_05$Trait <- "BD"
res_All_sublineage_05 <- rbind(res_All_sublineage_05, bd_sublineage_05)

bd_sublineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_01, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
bd_sublineage_01 <- bd_sublineage_01$results[order(bd_sublineage_01$results$p),c(1,3:6)]
bd_sublineage_01$Trait <- "BD"
res_All_sublineage_01 <- rbind(res_All_sublineage_01, bd_sublineage_01)

print("[INFO] Write to output...")
fwrite(res_All_lineage_05, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/All_Lineage_05.tsv", sep = "\t")
fwrite(res_All_lineage_01, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/All_Lineage_01.tsv", sep = "\t")
fwrite(res_All_sublineage_05, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/All_Sublineage_05.tsv", sep = "\t")
fwrite(res_All_sublineage_01, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/All_Sublineage_01.tsv", sep = "\t")

#Stage level
res_lineage_05 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(res_lineage_05) <- c("CellType", "p", "fold_change", "sd_from_mean", "q", "Stage", "Trait")

res_lineage_01 <- res_lineage_05

res_sublineage_05 <- res_lineage_05
res_sublineage_01 <- res_lineage_05

for (stage in c("firsttrim", "secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder")){
    print(paste0("[INFO] Analyze ", stage, "..."))
    print("[INFO] Lineage level!")
    data <- readRDS("/home/lilab/wuyong/project/scRNA/MP-R1/downsample500/downsample500_lineage.rds")
    bg_genes <- Features(data)
    data <- subset(data, stage == stage)
    count <- GetAssayData(data, assay = "RNA", layer = "counts")
    meta <- data@meta.data
    exp <- EWCE::sct_normalize(count)
    exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = exp, input_species = "human", output_species = "human", level2annot = meta$lineage, no_cores=40)
    annotLevels <- list(level1class=meta$lineage, level2class=meta$sub_lineage)
    ctd <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = paste0(stage, "_lineage"), savePath="/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result", no_cores=40)
    
    load(ctd)
    sz_fc_lineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_05, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
    sz_fc_lineage_05 <- sz_fc_lineage_05$results[order(sz_fc_lineage_05$results$p),c(1,3:6)]
    sz_fc_lineage_05$Stage <- stage
    sz_fc_lineage_05$Trait <- "SCZ"
    res_lineage_05 <- rbind(res_lineage_05, sz_fc_lineage_05)

    sz_fc_lineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_01, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
    sz_fc_lineage_01 <- sz_fc_lineage_01$results[order(sz_fc_lineage_01$results$p),c(1,3:6)]
    sz_fc_lineage_01$Stage <- stage
    sz_fc_lineage_01$Trait <- "SCZ"
    res_lineage_01 <- rbind(res_lineage_01, sz_fc_lineage_01)

    bd_fc_lineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_05, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
    bd_fc_lineage_05 <- bd_fc_lineage_05$results[order(bd_fc_lineage_05$results$p),c(1,3:6)]
    bd_fc_lineage_05$Stage <- stage
    bd_fc_lineage_05$Trait <- "BD"
    res_lineage_05 <- rbind(res_lineage_05, bd_fc_lineage_05)

    bd_fc_lineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_01, bg = bg_genes, reps = reps, annotLevel = 1, no_cores = 40)
    bd_fc_lineage_01 <- bd_fc_lineage_01$results[order(bd_fc_lineage_01$results$p),c(1,3:6)]
    bd_fc_lineage_01$Stage <- stage
    bd_fc_lineage_01$Trait <- "BD"
    res_lineage_01 <- rbind(res_lineage_01, bd_fc_lineage_01)


    #sub lineage level
    print("[INFO] Sublineage level!")
    data_sub <- readRDS(paste0("/home/lilab/wuyong/project/scRNA/MP-R1/downsample500/downsample500_sublineage.rds"))
    data_sub <- subset(data_sub, stage == stage)
    count_sub <- GetAssayData(data_sub, assay = "RNA", layer = "counts")
    meta_sub <- data_sub@meta.data
    exp_sub <- EWCE::sct_normalize(count_sub)
    exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = exp_sub, input_species = "human", output_species = "human", level2annot = meta_sub$sub_lineage, no_cores=40)
    annotLevels <- list(level1class=meta_sub$lineage, level2class=meta_sub$sub_lineage)
    ctd_sub <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = paste0(stage, "_sublineage"), savePath="/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result", no_cores=40)

    load(ctd_sub)
    sz_fc_sublineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_05, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
    sz_fc_sublineage_05 <- sz_fc_sublineage_05$results[order(sz_fc_sublineage_05$results$p),c(1,3:6)]
    sz_fc_sublineage_05$Stage <- stage
    sz_fc_sublineage_05$Trait <- "SCZ"
    res_sublineage_05 <- rbind(res_sublineage_05, sz_fc_sublineage_05)

    sz_fc_sublineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = sz_genes_01, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
    sz_fc_sublineage_01 <- sz_fc_sublineage_01$results[order(sz_fc_sublineage_01$results$p),c(1,3:6)]
    sz_fc_sublineage_01$Stage <- stage
    sz_fc_sublineage_01$Trait <- "SCZ"
    res_sublineage_01 <- rbind(res_sublineage_01, sz_fc_sublineage_01)

    bd_fc_sublineage_05 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_05, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
    bd_fc_sublineage_05 <- bd_fc_sublineage_05$results[order(bd_fc_sublineage_05$results$p),c(1,3:6)]
    bd_fc_sublineage_05$Stage <- stage
    bd_fc_sublineage_05$Trait <- "BD"
    res_sublineage_05 <- rbind(res_sublineage_05, bd_fc_sublineage_05)

    bd_fc_sublineage_01 <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_genes_01, bg = bg_genes, reps = reps, annotLevel = 2, no_cores = 40)
    bd_fc_sublineage_01 <- bd_fc_sublineage_01$results[order(bd_fc_sublineage_01$results$p),c(1,3:6)]
    bd_fc_sublineage_01$Stage <- stage
    bd_fc_sublineage_01$Trait <- "BD"
    res_sublineage_01 <- rbind(res_sublineage_01, bd_fc_sublineage_01)
}

print("[INFO] Write to output...")
fwrite(res_lineage_05, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/Stage_Lineage_05.tsv", sep = "\t")
fwrite(res_lineage_01, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/Stage_Lineage_01.tsv", sep = "\t")
fwrite(res_sublineage_05, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/Stage_Sublineage_05.tsv", sep = "\t")
fwrite(res_sublineage_01, file = "/home/lilab/wuyong/project/scRNA/MP-R1/EWCE-result/Stage_Sublineage_01.tsv", sep = "\t")
print("[INFO] Finished!!!")

