library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(tidyverse)

prepareDataFromscRNA1 <- function (object = NULL, diffData = NULL, showAverage = TRUE, 
          cells = NULL, group.by = "ident", assays = "RNA", slot = "data", 
          scale.data = TRUE, cluster.order = NULL, keep.uniqGene = TRUE, 
          sep = "_") 
{
  markerGene <- unique(diffData$gene)
  if (showAverage == TRUE) {
    vr <- utils::compareVersion(as.character(utils::packageVersion("Seurat")), 
                                "5")
    if (vr == 1) {
      mean_gene_exp <- as.matrix(data.frame(Seurat::AverageExpression(object, 
                                                                      features = markerGene, group.by = group.by, assays = assays, 
                                                                      layer = slot)))
    }
    else {
      mean_gene_exp <- as.matrix(data.frame(Seurat::AverageExpression(object, 
                                                                      features = markerGene, group.by = group.by, assays = assays, 
                                                                      slot = slot)))
    }
    name1 <- gsub(pattern = paste0(assays, ".", sep = ""), 
                  replacement = "", colnames(mean_gene_exp))
    colnames(mean_gene_exp) <- gsub(pattern = "\\.", replacement = " ", 
                                    name1)
    if (scale.data == TRUE) {
      mean_gene_exp <- t(scale(t(mean_gene_exp)))
    }
    if (!is.null(cluster.order)) {
      mean_gene_exp <- mean_gene_exp[, cluster.order]
    }
    geneMode <- "average"
  }
  else {
    cell.order <- data.frame(cell.id = names(Seurat::Idents(object)), 
                             cell.ident = Seurat::Idents(object))
    if (is.null(cluster.order)) {
      cell.order$cell.ident <- factor(cell.order$cell.ident, 
                                      levels = levels(Seurat::Idents(object)))
    }
    else {
      cell.order$cell.ident <- factor(cell.order$cell.ident, 
                                      levels = cluster.order)
    }
    cell.order <- cell.order[order(cell.order$cell.ident), 
    ]
    getassy <- as.matrix(Seurat::GetAssayData(object = object, 
                                              slot = slot)[features = markerGene, cells = NULL, 
                                                           drop = FALSE])
    id.order <- match(cell.order$cell.id, colnames(getassy))
    getassy <- getassy[, id.order]
    colnames(getassy) <- paste(colnames(getassy), cell.order$cell.ident, 
                               sep = "|")
    mean_gene_exp <- getassy
    if (scale.data == TRUE) {
      mean_gene_exp <- t(scale(t(mean_gene_exp)))
    }
    geneMode <- "all"
  }
  merMat <- data.frame(mean_gene_exp, check.names = FALSE)
  merMat$gene <- rownames(merMat)
  cinfo.gene <- diffData[, c("cluster", "gene")]
  cn <- unique(cinfo.gene$cluster)
  wide.res <- purrr::map_df(seq_along(cn), function(x) {
    tmp <- cinfo.gene[which(cinfo.gene$cluster == cn[x]), 
    ]
    tmp2 <- dplyr::mutate(merMat[which(merMat$gene %in% tmp$gene), 
    ], cluster = as.character(x))
    return(tmp2)
  })
  if (keep.uniqGene == TRUE) {
    duplicated_genes <- duplicated(wide.res$gene)
    wide.res <- wide.res[!duplicated_genes, ]
    geneType <- paste0("unique", "|", sep)
  }
  else {
    wide.res$gene <- make.unique(wide.res$gene, sep = sep)
    geneType <- paste0("nounique", "|", sep)
  }
  df <- reshape2::melt(wide.res, id.vars = c("cluster", "gene"), 
                       variable.name = "cell_type", value.name = "norm_value")
  df$cluster_name <- paste("cluster ", df$cluster, sep = "")
  if (showAverage == FALSE) {
    df$cell_type <- vapply(strsplit(as.character(df$cell_type), 
                                    split = "\\|"), function(x) {
                                      x[2]
                                    }, character(1))
  }
  cl.info <- dplyr::arrange(dplyr::mutate(data.frame(table(wide.res$cluster)), 
                                          Var1 = as.numeric(as.character(Var1))), Var1)
  id <- unique(df$cluster_name)
  df <- purrr::map_df(seq_along(id), function(x) {
    tmp <- dplyr::filter(df, cluster_name == id[x])
    dplyr::mutate(tmp, cluster_name = paste(cluster_name, 
                                            " (", cl.info$Freq[x], ")", sep = ""))
  })
  df$cluster_name <- factor(df$cluster_name, levels = paste("cluster ", 
                                                            cl.info$Var1, " (", cl.info$Freq, ")", sep = ""))
  return(list(wide.res = wide.res, long.res = df, type = "scRNAdata", 
              geneMode = geneMode, geneType = geneType))
}

setwd("D:/CellType/downsample500/ClusterGVis")
data <- readRDS("../downsample500_sublineage.rds")

#In-PAX6
data_Pax6 <- subset(data, sub_lineage == "In-PAX6")

#differential expression
pax6.markers <- FindAllMarkers(data_Pax6, min.pct = 0.25, logfc.threshold = 0.25, group.by = "stage", only.pos = TRUE)
pax6.markers.group <- pax6.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(pax6.markers.group, "In-PAX6-stage-diff.csv")

Idents(data_Pax6) <- "stage"
st.data <- prepareDataFromscRNA1(object = data_Pax6, diffData = pax6.markers.group, showAverage = TRUE, keep.uniqGene = FALSE, group.by = "stage")

#enrichment analysis
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 1000,
                        seed = 123)
head(enrich)
enrich_top <- enrich %>% group_by(group) %>% top_n(n = -5, wt = pvalue)

enrich$group <- factor(enrich$group, labels = levels(Seurat::Idents(data_Pax6)))
write.csv(enrich, "In-PAX6-stage-enrich.csv")

#heatmap
pdf('In-Pax6.pdf',height = 12,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           annoTerm.data = enrich_top,
           line.side = "left",
           sample.order = c("secondtrim", "thirdtrim", "infant", "juvenile", "youth", "midlife", "elder"),
           cluster.order = c(5, 4, 3, 1, 2, 6, 7)
           )
dev.off()


#In-VIP
data_vip <- subset(data, sub_lineage == "In-VIP")

#differential expression
vip.markers <- FindAllMarkers(data_vip, min.pct = 0.25, logfc.threshold = 0.25, group.by = "stage", only.pos = TRUE)
vip.markers.group <- vip.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(vip.markers.group, "In-VIP-stage-diff.csv")

Idents(data_vip) <- "stage"
st.data <- prepareDataFromscRNA1(object = data_vip, diffData = vip.markers.group, showAverage = TRUE, keep.uniqGene = FALSE, group.by = "stage")

#enrichment analysis
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 1000,
                        seed = 123)
head(enrich)
enrich_top <- enrich %>% group_by(group) %>% top_n(n = -5, wt = pvalue)

enrich$group <- factor(enrich$group, labels = levels(Seurat::Idents(data_vip)))
write.csv(enrich, "In-VIP-stage-enrich.csv")

#heatmap
pdf('In-vip.pdf',height = 12,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           annoTerm.data = enrich_top,
           line.side = "left",
           sample.order = c("secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"),
           cluster.order = c(7, 6, 5, 4, 1, 2, 3, 8, 9)
)
dev.off()



#In-SST
data_sst <- subset(data, sub_lineage == "In-SST")

#differential expression
sst.markers <- FindAllMarkers(data_sst, min.pct = 0.25, logfc.threshold = 0.25, group.by = "stage", only.pos = TRUE)
sst.markers.group <- sst.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(sst.markers.group, "In-SST-stage-diff.csv")

Idents(data_sst) <- "stage"
st.data <- prepareDataFromscRNA1(object = data_sst, diffData = sst.markers.group, showAverage = TRUE, keep.uniqGene = FALSE, group.by = "stage")

#enrichment analysis
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 1000,
                        seed = 123)
head(enrich)
enrich_top <- enrich %>% group_by(group) %>% top_n(n = -5, wt = pvalue)

enrich$group <- factor(enrich$group, labels = levels(Seurat::Idents(data_sst)))
write.csv(enrich, "In-SST-stage-enrich.csv")

#heatmap
pdf('In-sst.pdf',height = 12,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           annoTerm.data = enrich_top,
           line.side = "left",
           sample.order = c("secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"),
           cluster.order = c(8, 7, 6, 4, 1, 2, 3, 5, 9)
)
dev.off()


#Ex-L56NP
data_l56np <- subset(data, sub_lineage == "Ex-L56NP")

#differential expression
l56np.markers <- FindAllMarkers(data_l56np, min.pct = 0.25, logfc.threshold = 0.25, group.by = "stage", only.pos = TRUE)
l56np.markers.group <- l56np.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(l56np.markers.group, "Ex-L56NP-stage-diff.csv")

Idents(data_l56np) <- "stage"
st.data <- prepareDataFromscRNA1(object = data_l56np, diffData = l56np.markers.group, showAverage = TRUE, keep.uniqGene = FALSE, group.by = "stage")

#enrichment analysis
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 1000,
                        seed = 123)
head(enrich)
enrich_top <- enrich %>% group_by(group) %>% top_n(n = -5, wt = pvalue)

enrich$group <- factor(enrich$group, labels = levels(Seurat::Idents(data_l56np)))
write.csv(enrich, "Ex-L56NP-stage-enrich.csv")

#heatmap
pdf('In-l56np.pdf',height = 12,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           annoTerm.data = enrich_top,
           line.side = "left",
           sample.order = c("secondtrim", "thirdtrim", "infant", "juvenile", "youth", "midlife", "elder"),
           cluster.order = c(5, 4, 3, 1, 2, 6, 7)
)
dev.off()



#Ex-L6CT
data_l6ct <- subset(data, sub_lineage == "Ex-L6CT")

#differential expression
l6ct.markers <- FindAllMarkers(data_l6ct, min.pct = 0.25, logfc.threshold = 0.25, group.by = "stage", only.pos = TRUE)
l6ct.markers.group <- l6ct.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(l6ct.markers.group, "Ex-L6CT-stage-diff.csv")

Idents(data_l6ct) <- "stage"
st.data <- prepareDataFromscRNA1(object = data_l6ct, diffData = l6ct.markers.group, showAverage = TRUE, keep.uniqGene = FALSE, group.by = "stage")

#enrichment analysis
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 1000,
                        seed = 123)
head(enrich)
enrich_top <- enrich %>% group_by(group) %>% top_n(n = -5, wt = pvalue)

enrich$group <- factor(enrich$group, labels = levels(Seurat::Idents(data_l6ct)))
write.csv(enrich, "Ex-L6CT-stage-enrich.csv")

#heatmap
pdf('In-l6ct.pdf',height = 12,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           annoTerm.data = enrich_top,
           line.side = "left",
           sample.order = c("secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder"),
           cluster.order = c(8, 7, 6, 4, 1, 2, 3, 5, 9)
)
dev.off()