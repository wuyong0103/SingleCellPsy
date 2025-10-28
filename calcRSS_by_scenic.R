library(optparse)
op_list <- list(
    make_option(c("-l", "--input_loom"), type = "character", default = NULL, action = "store", help = "The input of aucell loom file", metavar="rds"),
    make_option(c("-m", "--input_meta"), type = "character", default = NULL, action = "store", help = "The metadata of Seurat object", metavar="idents"),
    make_option(c("-s", "--stage"), type = "character", default = NULL, action = "store", help = "The colname of metadata to calculate RSS", metavar="stage"),
    make_option(c("-c", "--celltype"), type = "character", default = NULL, action = "store", help = "Cell type", metavar="label")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)

loom <- open_loom(opt$input_loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

meta <- read.table(opt$input_meta, sep="\t", header=T)
meta <- meta[meta$sub_lineage == opt$celltype, ]

cellinfo <- meta[,c(opt$stage, "nFeature_RNA", "nCount_RNA")]
colnames(cellinfo)=c('stage', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo, select = 'stage'))
selectedResolution <- "stage"

sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution])
rss=na.omit(rss)
try({
    rssPlot <- plotRSS(rss)
    rssfile <- paste0(opt$celltype, "_regulon_RSS.Rdata")
    save(regulonAUC,rssPlot,regulons,file=rssfile)
})

saveRDS(rss,paste0(opt$celltype, "_rss.rds"))
