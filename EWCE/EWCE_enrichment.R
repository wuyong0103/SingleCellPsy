library(EWCE)
library(ggplot2)
library(tidyverse)
deg <- read_csv("Xia2024STM.csv")
deg <- deg %>% select(ensg, hgnc, ASD.logFC, ASD.fdr, SCZ.logFC, SCZ.fdr, BD.logFC, BD.fdr, FASD.logFC, FASD.fdr, 
						FSCZ.logFC, FSCZ.fdr, FBD.logFC, FBD.fdr, MASD.logFC, MASD.fdr, MSCZ.logFC, MSCZ.fdr, MBD.logFC, MBD.fdr)
mdd <- read_csv("Daskalakis2024Science.csv")
mdd_diff <- mdd %>% filter(adj.P.Val < 0.05)

asd_diff <- deg %>% filter(ASD.fdr < 0.05)
scz_diff <- deg %>% filter(SCZ.fdr < 0.05)
bd_diff <- deg %>% filter(BD.fdr < 0.05)

asd_ex <- read_table("ASC_Exome.txt")
asd_ex <- asd_ex %>% filter(p<0.05)
asd_sfari <- read_csv("SFARI_ASD.csv")
asd_sfari <- asd_sfari %>% filter(asd_sfari$syndromic == 1 | asd_sfari$score==1)
bd_ex <- read_table("BipEx_Exome.txt")
bd_ex <- bd_ex %>% filter(p<0.05)
scz_ex <- read_table("SCHEMA_Exome.txt")
scz_ex <- scz_ex %>% filter(p<0.05)

sc_bg <- read.table("scRNA_All.txt", header = F)

for (cf in c("All", "male", "female", "trimester2nd", "trimester3rd", "years0_1", "years1_2", "years2_4", "years4_10", "years10_20", "Adult", "trimester2ndMale", "trimester2ndFemale", 
            "trimester3rdMale", "trimester3rdFemale", "years0_1Male", "years0_1Female", "years1_2Male", "years1_2Female", "years2_4Male", "years2_4Female", "years4_10Male", 
            "years4_10Female","years10_20Male", "years10_20Female", "AdultMale", "AdultFemale")){

	ctd_file <- paste("ctd_Cortex_", cf, ".rda", sep = "")
	load(ctd_file)
	
	#ASD
	asd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = asd_sfari$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(asd_results$results)
	asd_ewce <- asd_results$results[order(asd_results$results$p),3:6]
	asd_sfari_lvl2 <- paste("result/ASD_SFARI_", cf, "_lvl2.csv", sep = "")
	write.csv(asd_ewce, file = asd_sfari_lvl2)
	asd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = asd_sfari$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(asd_results$results)
	asd_ewce <- asd_results$results[order(asd_results$results$p),3:6]
	asd_sfari_lvl1 <- paste("result/ASD_SFARI_", cf, "_lvl1.csv", sep = "")
	write.csv(asd_ewce, file = asd_sfari_lvl1)

	asd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = asd_ex$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(asd_results$results)
	asd_ewce <- asd_results$results[order(asd_results$results$p),3:6]
	asd_ex_lvl2 <- paste("result/ASD_Exome_", cf, "_lvl2.csv", sep = "")
	write.csv(asd_ewce, file = asd_ex_lvl2)
	asd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = asd_ex$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(asd_results$results)
	asd_ewce <- asd_results$results[order(asd_results$results$p),3:6]
	asd_ex_lvl1 <- paste("result/ASD_Exome_", cf, "_lvl1.csv", sep = "")
	write.csv(asd_ewce, file = asd_ex_lvl1)
	
	asd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = asd_diff$hgnc, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(asd_results$results)
	asd_ewce <- asd_results$results[order(asd_results$results$p),3:6]
	asd_diff_lvl2 <- paste("result/ASD_Diff_", cf, "_lvl2.csv", sep = "")
	write.csv(asd_ewce, file = asd_diff_lvl2)
	asd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = asd_diff$hgnc, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(asd_results$results)
	asd_ewce <- asd_results$results[order(asd_results$results$p),3:6]
	asd_diff_lvl1 <- paste("result/ASD_Diff_", cf, "_lvl1.csv", sep = "")
	write.csv(asd_ewce, file = asd_diff_lvl1)
	
	#BD
	bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_ex$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(bd_results$results)
	bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
	bd_ex_lvl2 <- paste("result/BD_Exome_", cf, "_lvl2.csv", sep = "")
	write.csv(bd_ewce, file = bd_ex_lvl2)
	bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_ex$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(bd_results$results)
	bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
	bd_ex_lvl1 <- paste("result/BD_Exome_", cf, "_lvl1.csv", sep = "")
	write.csv(bd_ewce, file = bd_ex_lvl1)
	
	bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_diff$hgnc, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(bd_results$results)
	bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
	bd_diff_lvl2 <- paste("result/BD_Diff_", cf, "_lvl2.csv", sep = "")
	write.csv(bd_ewce, file = bd_diff_lvl2)
	bd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = bd_diff$hgnc, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(bd_results$results)
	bd_ewce <- bd_results$results[order(bd_results$results$p),3:6]
	bd_diff_lvl1 <- paste("result/BD_Diff_", cf, "_lvl1.csv", sep = "")
	write.csv(bd_ewce, file = bd_diff_lvl1)
	
	#SCZ
	scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_ex$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(scz_results$results)
	scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
	scz_ex_lvl2 <- paste("result/SCZ_Exome_", cf, "_lvl2.csv", sep = "")
	write.csv(scz_ewce, file = scz_ex_lvl2)
	scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_ex$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(scz_results$results)
	scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
	scz_ex_lvl1 <- paste("result/SCZ_Exome_", cf, "_lvl1.csv", sep = "")
	write.csv(scz_ewce, file = scz_ex_lvl1)
	
	scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_diff$hgnc, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(scz_results$results)
	scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
	scz_diff_lvl2 <- paste("result/SCZ_Diff_", cf, "_lvl2.csv", sep = "")
	write.csv(scz_ewce, file = scz_diff_lvl2)
	scz_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = scz_diff$hgnc, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(scz_results$results)
	scz_ewce <- scz_results$results[order(scz_results$results$p),3:6]
	scz_diff_lvl1 <- paste("result/SCZ_Diff_", cf, "_lvl1.csv", sep = "")
	write.csv(scz_ewce, file = scz_diff_lvl1)
	
	
	#MDD
	mdd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = mdd_diff$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 2, no_cores = 20)
	knitr::kable(mdd_results$results)
	mdd_ewce <- mdd_results$results[order(mdd_results$results$p),3:6]
	mdd_diff_lvl2 <- paste("result/MDD_Diff_", cf, "_lvl2.csv", sep = "")
	write.csv(mdd_ewce, file = mdd_diff_lvl2)
	mdd_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, sctSpecies = "human", genelistSpecies = "human", hits = mdd_diff$symbol, bg = sc_bg$V1, reps = 10000, annotLevel = 1, no_cores = 20)
	knitr::kable(mdd_results$results)
	mdd_ewce <- mdd_results$results[order(mdd_results$results$p),3:6]
	mdd_diff_lvl1 <- paste("result/MDD_Diff_", cf, "_lvl1.csv", sep = "")
	write.csv(mdd_ewce, file = mdd_diff_lvl1)

    rm(ctd)
    rm(ctd_file)
}
