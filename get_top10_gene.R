library(tidyverse)
Args <- commandArgs(TRUE)
exp_prefix <- Args[1] #output from Calculate_MeanExpression.py
exp_file <- paste(exp_prefix, ".csv", sep="")
exp_CT <- read_csv(exp_file)
gene_coordinates <- read.table("/home/lilab/wuyong/project/scRNA/test/NCBI37.3.gene.loc.extendedMHCexcluded", header=F, stringsAsFactors = F) %>%
    mutate(start=ifelse(V3-100000<0,0,V3-100000),end=V4+100000,V1=as.character(V1)) %>%
    select(2,start,end,1) %>%
    as.tibble() %>%
    rename(chr="V2", ENTREZ="V1") %>%
    mutate(chr=paste0("chr",chr))
exp_CT <- exp_CT %>% group_by(Gene) %>% mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean))
entrez2symbol <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL2EG) %>% rename(Gene="symbol",ENTREZ="gene_id")
exp_CT <- inner_join(exp_CT, entrez2symbol, by="Gene")
exp_CT <- inner_join(exp_CT, gene_coordinates, by="ENTREZ")
n_genes <- length(unique(exp_CT$ENTREZ))
n_genes_to_keep <- (n_genes * 0.1) %>% round()
#save(exp_CT, file = paste(exp_prefix, ".Rdata", sep=""))

magma_top10 <- function(d){
  d_spe <- d %>% group_by(Lvl5) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group_magma(.))
}

write_group_magma  = function(df) {
  df <- select(df,Lvl5,ENTREZ)
  df_name <- make.names(unique(df[1]))
  colnames(df)[2] <- df_name  
  dir.create(paste0("MAGMA/"), showWarnings = FALSE)
  select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
  write_tsv("MAGMA/top10.txt",append=T)
return(df)
}

write_group  = function(df) {
  df <- select(df,Lvl5,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}

ldsc_bedfile <- function(d){
  d_spe <- d %>% group_by(Lvl5) %>% top_n(.,n_genes_to_keep,specificity)
  d_spe %>% do(write_group(.))
}

exp_CT %>% filter(Expr_sum_mean>1) %>% magma_top10()
exp_CT %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile()
