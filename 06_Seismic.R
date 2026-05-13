#2026-04-08
#wuyong0103@126.com

#Usage:
#Rscript 05_Seismic.R /home/lilab/wuyong/project/scRNA/MP-R1/downsampleMin/downsampleMin_lineage.rds lineage None All_lineage
#args[1]: Seurat object
#args[2]: column of cell type (一定是Seurat object meta.data中包含的列名)
#args[3]：age stage, or "None"
#args[4]：prefix of output


library(Seurat)
library(SingleCellExperiment)
library(seismicGWAS)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
args <- commandArgs(trailingOnly = TRUE)

print("[INFO] Loading GWAS data in magma format...")
BD2019Stahl <- read.table("/home/lilab/wuyong/data/GWAS/BD2019Stahl/BD2019Stahl.annotated_35kbup_10_down.genes.out", header = T)
BD2021PGC3 <- read.table("/home/lilab/wuyong/data/GWAS/BD2021PGC3/BD2021PGC3.annotated_35kbup_10_down.genes.out", header = T)
BD2025OConnel <- read.table("/home/lilab/wuyong/data/GWAS/BD2025OConnel/BD2025OConnel.annotated_35kbup_10_down.genes.out", header = T)
SCZ2014PGC2 <- read.table("/home/lilab/wuyong/data/GWAS/SCZ2014PGC2/SCZ2014PGC2.annotated_35kbup_10_down.genes.out", header = T)
SCZ2019Lam <- read.table("/home/lilab/wuyong/data/GWAS/SCZ2019Lam/SCZ2019Lam.annotated_35kbup_10_down.genes.out", header = T)
SCZ2022PGC3 <- read.table("/home/lilab/wuyong/data/GWAS/SCZ2022PGC3/SCZ2022PGC3.annotated_35kbup_10_down.genes.out", header = T)
SCZ2026Bigdeli <- read.table("/home/lilab/wuyong/data/GWAS/SCZ2026Bigdeli/SCZ2026Bigdeli.annotated_35kbup_10_down.genes.out", header = T)

print("[INFO] Loading scRNA data in Seurat format...")
data <- readRDS(args[1])
data <- JoinLayers(data)
if(args[3] != "None"){
    data <- subset(data, stage == args[3])
}
data <- DietSeurat(data, assays = "RNA")

print("[INFO] Transform Seurat object to SCE object...")
sce <- as.SingleCellExperiment(data)

print("[INFO] Calculating cell type specificity...")
tmfacs_sscore <- calc_specificity(sce, ct_label_col = args[2])
tmfacs_sscore_hsa <- translate_gene_ids(tmfacs_sscore, from='hsa_symbol')

dis_list <- list(
        BD2019Stahl=BD2019Stahl, 
        BD2021PGC3=BD2021PGC3, 
        BD2025OConnel=BD2025OConnel, 
        SCZ2014PGC2=SCZ2014PGC2, 
        SCZ2019Lam=SCZ2019Lam, 
        SCZ2022PGC3=SCZ2022PGC3,
        SCZ2026Bigdeli=SCZ2026Bigdeli
        )

res_asso <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(res_asso) <- c("cell_type", "pvalue", "FDR", "Trait")

res_gene <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(res_gene) <- c("gene", "specificity", "zstat", "dfbetas", "is_influential", "symbol", "Trait", "Celltype")

res_go <- data.frame(matrix(ncol = 15, nrow = 0))
colnames(res_go) <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Enrichment", "Trait", "Celltype")

res_kegg <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(res_kegg) <- c("category", "subcategory", "ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Trait", "Celltype")

#获取单细胞表达的所有基因作为富集分析的背景基因集
bg_genes <- rownames(sce)
bg_genes <- bitr(bg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

for (dis in names(dis_list)){
    print(paste0("[INFO] Cell type enrichment of ", dis, "..."))
    asso <- get_ct_trait_associations(tmfacs_sscore_hsa, dis_list[[dis]])
    asso$Trait <- dis
    if(nrow(asso > 0)){
        res_asso <- rbind(res_asso, asso)
    }
    
    sig <- asso$cell_type[asso$FDR<0.05]

    for (ct in sig){
        #获取受影响的基因
        sig_gene <- find_inf_genes(ct, tmfacs_sscore_hsa, dis_list[[dis]])
        hsa.map <- mapIds(org.Hs.eg.db, keys = sig_gene$gene, keytype = "ENTREZID", column = "SYMBOL")
        hsa.map <- stack(hsa.map)
        colnames(hsa.map) <- c("symbol", "entrezid")
        sig_gene <- merge(sig_gene, hsa.map, by.x='gene', by.y='entrezid')
        sig_gene$Trait <- dis
        sig_gene$Celltype <- ct
        sig_gene <- sig_gene[sig_gene$is_influential == "TRUE", ]
        if(nrow(sig_gene > 0)){
            res_gene <- rbind(res_gene, sig_gene)
        }

        #富集分析
        ego_bp <- enrichGO(gene = sig_gene$gene, universe = bg_genes$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
        ego_cc <- enrichGO(gene = sig_gene$gene, universe = bg_genes$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
        ego_mf <- enrichGO(gene = sig_gene$gene, universe = bg_genes$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
        ego_kegg <- enrichKEGG(gene = sig_gene$gene, universe = bg_genes$ENTREZID, pAdjustMethod = "BH", pvalueCutoff  = 0.01)
        ego_bp <- as.data.frame(ego_bp)
        ego_cc <- as.data.frame(ego_cc)
        ego_mf <- as.data.frame(ego_mf)
        ego_kegg <- as.data.frame(ego_kegg)
        if(nrow(ego_bp > 0)){
            ego_bp$Enrichment = "BP"
            ego_bp$Trait = dis
            ego_bp$Celltype = ct
            res_go <- rbind(res_go, ego_bp)
        }
        if(nrow(ego_cc > 0)){
            ego_cc$Enrichment = "CC"
            ego_cc$Trait = dis
            ego_cc$Celltype = ct
            res_go <- rbind(res_go, ego_cc)
        }
        if(nrow(ego_mf > 0)){
            ego_mf$Enrichment = "MF"
            ego_mf$Trait = dis
            ego_mf$Celltype = ct
            res_go <- rbind(res_go, ego_mf)
        }
        if(nrow(ego_kegg > 0)){
            ego_kegg$Trait = dis
            ego_kegg$Celltype = ct
            res_kegg <- rbind(res_kegg, ego_kegg)
        }
    }
    print(paste0("[INFO] Cell type enrichment of ", dis, " finished!"))
}

print("[INFO] Writing output...")
fwrite(res_asso, file = paste0("/home/lilab/wuyong/project/scRNA/MP-R1/Seismic-result/", args[4], "_Asso.tsv"), sep = "\t")
fwrite(res_gene, file = paste0("/home/lilab/wuyong/project/scRNA/MP-R1/Seismic-result/", args[4], "_Gene.tsv"), sep = "\t")
fwrite(res_go, file = paste0("/home/lilab/wuyong/project/scRNA/MP-R1/Seismic-result/", args[4], "_GO.tsv"), sep = "\t")
fwrite(res_kegg, file = paste0("/home/lilab/wuyong/project/scRNA/MP-R1/Seismic-result/", args[4], "_KEGG.tsv"), sep = "\t")
print("[INFO] Finished all analysis!!!")
