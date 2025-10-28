library(EWCE)
library(Seurat)
library(tidyverse)
library(ggplot2)
set.seed(123)
setwd("/home/lilab/wuyong/project/scRNA/integrate/downsample500")

data_lineage <- readRDS("downsample500_lineage.rds")
meta_lineage <- read.table("downsample500_lineage.txt", sep="\t", header=T)
data_sublineage <- readRDS("downsample500_sublineage.rds")
meta_sublineage <- read.table("downsample500_sublineage.txt", sep="\t", header=T)
data_all <- readRDS("../downsampleMin/downsampleMin_lineage.rds")
meta_all <- read.table("../downsampleMin/downsampleMin_lineage.txt", sep="\t", header=T)
data_suball <- readRDS("../downsampleMin/downsampleMin_sublineage.rds")
meta_suball <- read.table("../downsampleMin/downsampleMin_sublineage.txt", sep="\t", header=T)

bg <- Features(data_lineage)

bd_ex <- read_table("../BipEx_Exome.txt")
bd_ex <- bd_ex %>% filter(p<0.05)
scz_ex <- read_table("../SCHEMA_Exome.txt")
scz_ex <- scz_ex %>% filter(p<0.05)

data_all <- JoinLayers(data_all)
mtx <- as.matrix(data_all@assays$RNA@layers$data)
colnames(mtx) <- colnames(data_all)
rownames(mtx) <- rownames(data_all)
mtx <- fix_bad_hgnc_symbols(mtx)
m <- match(colnames(mtx), meta_all$cell)
meta <- meta_all[m,]
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = mtx, sctSpecies_origin = "human", input_species = "human", output_species = "human", level2annot = meta$sub_lineage, no_cores=80)
annotLevels <- list(level1class=meta$lineage,level2class=meta$sub_lineage)
groupname <- "downsampleMin_All"
ctd_lineage <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = groupname, savePath="/home/lilab/wuyong/project/scRNA/integrate/downsample500/EWCE-result", no_cores=80)

load("./EWCE-result/ctd_downsampleMin_All.rda")
bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_ex$symbol, bg = bg, reps = 100000, annotLevel = 1, no_cores = 20)
knitr::kable(bd_results$results)
bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
bd_ex_lvl1 <- "./EWCE-result/BD_Exome_downsampleMin_lvl1_All.csv"
write.csv(bd_ewce, file = bd_ex_lvl1)
scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_ex$symbol, bg = bg, reps = 100000, annotLevel = 1, no_cores = 20)
knitr::kable(scz_results$results)
scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
scz_ex_lvl1 <- "./EWCE-result/SCZ_Exome_downsampleMin_lvl1_All.csv"
write.csv(scz_ewce, file = scz_ex_lvl1)

data_suball <- JoinLayers(data_suball)
mtx <- as.matrix(data_suball@assays$RNA@layers$data)
colnames(mtx) <- colnames(data_suball)
rownames(mtx) <- rownames(data_suball)
mtx <- fix_bad_hgnc_symbols(mtx)
m <- match(colnames(mtx), meta_suball$cell)
meta <- meta_suball[m,]
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = mtx, sctSpecies_origin = "human", input_species = "human", output_species = "human", level2annot = meta$sub_lineage, no_cores=80)
annotLevels <- list(level1class=meta$lineage,level2class=meta$sub_lineage)
groupname <- "downsampleMin_subAll"
ctd_sublineage <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = groupname, savePath="/home/lilab/wuyong/project/scRNA/integrate/downsample500/EWCE-result", no_cores=80)

load("./EWCE-result/ctd_downsampleMin_subAll.rda")
bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_ex$symbol, bg = bg, reps = 100000, annotLevel = 2, no_cores = 20)
knitr::kable(bd_results$results)
bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
bd_ex_lvl2 <- "./EWCE-result/BD_Exome_downsampleMin_lvl2_All.csv"
write.csv(bd_ewce, file = bd_ex_lvl2)
scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_ex$symbol, bg = bg, reps = 100000, annotLevel = 2, no_cores = 20)
knitr::kable(scz_results$results)
scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
scz_ex_lvl2 <- "./EWCE-result/SCZ_Exome_downsampleMin_lvl2_All.csv"
write.csv(scz_ewce, file = scz_ex_lvl2)


for (i in c("firsttrim", "secondtrim", "child", "juvenile", "youth", "kid", "midlife", "infant", "thirdtrim", "elder" )){
	data <- subset(data_lineage, stage == i)
	mtx <- as.matrix(data@assays$RNA@layers$data)
	colnames(mtx) <- colnames(data)
	rownames(mtx) <- rownames(data)
	meta <- meta_lineage[meta_lineage$stage == i, ]
	dim(mtx)
	mtx <- fix_bad_hgnc_symbols(mtx)
	m <- match(colnames(mtx), meta$cell)
	meta <- meta[m,]
	exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = mtx, sctSpecies_origin = "human", input_species = "human", output_species = "human", level2annot = meta$sub_lineage, no_cores=80)
	annotLevels <- list(level1class=meta$lineage,level2class=meta$sub_lineage)
	groupname <- paste0("downsample500_lineage_", i)
	ctd_lineage <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = groupname, savePath="/home/lilab/wuyong/project/scRNA/integrate/downsample500/EWCE-result", no_cores=80)

	data <- subset(data_sublineage, stage == i)
	mtx <- as.matrix(data@assays$RNA@layers$data)
	colnames(mtx) <- colnames(data)
	rownames(mtx) <- rownames(data)
	meta <- meta_sublineage[meta_sublineage$stage == i, ]
	dim(mtx)
	mtx <- fix_bad_hgnc_symbols(mtx)
	m <- match(colnames(mtx), meta$cell)
	meta <- meta[m,]
	exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = mtx, sctSpecies_origin = "human", input_species = "human", output_species = "human", level2annot = meta$sub_lineage, no_cores=80)
	annotLevels <- list(level1class=meta$lineage,level2class=meta$sub_lineage)
	groupname <- paste0("downsample500_sublineage_", i)
	ctd_sublineage <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = groupname, savePath="/home/lilab/wuyong/project/scRNA/integrate/downsample500/EWCE-result", no_cores=80)


	rds_lineage <- paste0("./EWCE-result/ctd_downsample500_lineage_", i, ".rda")
	load(rds_lineage)
	bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_ex$symbol, bg = bg, reps = 100000, annotLevel = 1, no_cores = 20)
	knitr::kable(bd_results$results)
	bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
	bd_ex_lvl1 <- paste0("./EWCE-result/BD_Exome_downsample500_lvl1_", i, ".csv")
	write.csv(bd_ewce, file = bd_ex_lvl1)
	scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_ex$symbol, bg = bg, reps = 100000, annotLevel = 1, no_cores = 20)
	knitr::kable(scz_results$results)
	scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
	scz_ex_lvl1 <- paste0("./EWCE-result/SCZ_Exome_downsample500_lvl1_", i, ".csv")
	write.csv(scz_ewce, file = scz_ex_lvl1)
	
	rds_sublineage <- paste0("./EWCE-result/ctd_downsample500_sublineage_", i, ".rda")
	load(rds_sublineage)
	bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_ex$symbol, bg = bg, reps = 100000, annotLevel = 2, no_cores = 20)
	knitr::kable(bd_results$results)
	bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
	bd_ex_lvl2 <- paste0("./EWCE-result/BD_Exome_downsample500_lvl2_", i, ".csv")
	write.csv(bd_ewce, file = bd_ex_lvl2)
	scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_ex$symbol, bg = bg, reps = 100000, annotLevel = 2, no_cores = 20)
	knitr::kable(scz_results$results)
	scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
	scz_ex_lvl2 <- paste0("./EWCE-result/SCZ_Exome_downsample500_lvl2_", i, ".csv")
	write.csv(scz_ewce, file = scz_ex_lvl2)
}
