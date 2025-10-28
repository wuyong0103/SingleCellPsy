library(Seurat)
library(Matrix)
library(tidyverse)
library(dplyr)

# read the seurat object
data <- readRDS("Integrate_RunUMAP.rds")

# read the meta file
meta <- read_tsv("Integrate-Meta.tsv")
data@meta.data$lineage <- meta$lineage
data@meta.data$sub_lineage <- meta$sub_lineage
data@meta.data$cell <- meta$cell
dim(data)

seurat_obj <- subset(data, lineage != "Unknown")
seurat_sub <- subset(data, sub_lineage != "Unknown")
rm(data)
gc()

#=========Regardless age stage, downsample to minimum cell type================
Idents(seurat_sub) <- seurat_sub$sub_lineage
tab <- table(Idents(seurat_sub))
n_target <- min(tab)
seurat_sub_Min <- subset(seurat_sub, downsample = n_target)
table(Idents(seurat_sub_Min))
write.table(seurat_sub_Min@meta.data, "./downsampleMin/downsampleMin_sublineage.txt")
saveRDS(seurat_sub_Min, "./downsampleMin/downsampleMin_sublineage.rds")

Idents(seurat_obj) <- seurat_obj$lineage
tab <- table(Idents(seurat_obj))
n_target <- min(tab)
seurat_obj_Min <- subset(seurat_obj, downsample = n_target)
table(Idents(seurat_obj_Min))
write.table(seurat_obj_Min@meta.data, "./downsampleMin/downsampleMin_lineage.txt")
saveRDS(seurat_obj_Min, "./downsampleMin/downsampleMin_lineage.rds")


#========= downsample to 500 cells ======================================
meta_filtered_sub <- seurat_sub@meta.data %>%
    group_by(stage, sub_lineage) %>%
    filter(n() >= 500)

max_attempts <- 1000
attempt <- 1
keep_cells <- generate_keep_cells_sub(meta_filtered_sub, 500)
while (attempt <= max_attempts) {
    donor_counts <- table(meta_filtered_sub$donorid[meta_filtered_sub$cell %in% keep_cells])
    single_donors <- names(donor_counts[donor_counts == 1])
    
    if (length(single_donors) == 0) {
		break
    } else {
		keep_cells <- generate_keep_cells_sub(meta_filtered_sub, 500)
		attempt <- attempt + 1
    }
}

if (attempt > max_attempts) {
    warning("达到最大尝试次数，仍有 donorid 只剩 1 个细胞。继续使用当前下采样结果（sub_lineage）。")
}

seurat_sub_500 <- subset(seurat_sub, cells = keep_cells)
seurat_sub_500 <- JoinLayers(seurat_sub_500)

saveRDS(seurat_sub_500, "./downsample500/downsample500_sublineage.rds")
with(seurat_sub_500@meta.data, table(stage, sub_lineage))
write.table(seurat_sub_500@meta.data, "./downsample500/downsample500_sublineage.txt")


#lineage downsample
meta_filtered_lineage <- seurat_obj@meta.data %>%
    group_by(stage, lineage) %>%
    filter(n() >= 500)

attempt <- 1
keep_cells <- generate_keep_cells_lineage(meta_filtered_lineage, 500)
while (attempt <= max_attempts) {
    donor_counts <- table(meta_filtered_lineage$donorid[meta_filtered_lineage$cell %in% keep_cells])
    single_donors <- names(donor_counts[donor_counts == 1])
    
    if (length(single_donors) == 0) {
        break
    } else {
        keep_cells <- generate_keep_cells_lineage(meta_filtered_lineage, 500)
        attempt <- attempt + 1
    }
}

if (attempt > max_attempts) {
      warning("达到最大尝试次数，仍有 donorid 只剩 1 个细胞。继续使用当前下采样结果（lineage）。")
}

seurat_obj_500 <- subset(seurat_obj, cells = keep_cells)
seurat_obj_500 <- JoinLayers(seurat_obj_500)

saveRDS(seurat_obj_500, "./downsample500/downsample500_lineage.rds")
with(seurat_obj_500@meta.data, table(stage, lineage))
write.table(seurat_obj_500@meta.data, "./downsample500/downsample500_lineage.txt")


# function: calculate the average expression and cell type specificity, remove cell cluster with cell number less than 20
calculate_expression_and_specificity <- function(seurat_obj, group_by, prefix, downdir) {
  # create group column (for multi-var group, merge to single column)
  if (length(group_by) > 1) {
    group_col <- apply(seurat_obj@meta.data[, group_by], 1, paste, collapse = "_")
    seurat_obj@meta.data$temp_group <- group_col
    group_var <- "temp_group"
    filter_var <- group_by[1]
  } else {
    group_var <- group_by
    filter_var <- group_var
  }

  # check the cell number in each group
  cell_counts <- table(seurat_obj@meta.data[[group_var]])
  valid_groups <- names(cell_counts[cell_counts >= 20])
  if (length(valid_groups) == 0) {
    message(paste("No groups with >= 20 cells for", paste(group_by, collapse = "_")))
    return(NULL)
  }

  # filter Seurat object
  valid_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[group_var]] %in% valid_groups]
  sub_seurat <- subset(seurat_obj, cells = valid_cells)

  # Ensure RNA assay is active and populated
  if (!"RNA" %in% names(sub_seurat@assays)) {
      message("RNA assay not found in Seurat object")
      return(NULL)
  }
  DefaultAssay(sub_seurat) <- "RNA"

  #for multi-var group
  if (length(group_by) > 1) {
    group_levels <- unique(sub_seurat@meta.data[[filter_var]])
    results <- list(avg = list(), spec = list())
    
    for (level in group_levels) {
      level_cells <- colnames(sub_seurat)[sub_seurat@meta.data[[filter_var]] == level]
      if (length(level_cells) == 0){
          message(paste("No cells found for", filter_var, "=", level))
          next
      }
      
      level_seurat <- subset(sub_seurat, cells = level_cells)
      
      level_counts <- table(level_seurat@meta.data$temp_group)
      valid_level_groups <- names(level_counts[level_counts >= 20])
      if (length(valid_level_groups) == 0){
          message(paste("No valid groups with >= 20 cells for", filter_var, "=", level))
          next
      }
      
      level_seurat <- subset(level_seurat, cells = colnames(level_seurat)[level_seurat@meta.data$temp_group %in% valid_level_groups])

      # Check if level_seurat has cells
      if (ncol(level_seurat) == 0) {
          message(paste("No cells remain after filtering for", filter_var, "=", level))
          next
      }
      
      avg_expr <- AggregateExpression(level_seurat, group.by = "temp_group", assays = "RNA")$RNA
      if (is.null(avg_expr) || nrow(avg_expr) == 0) {
          message(paste("No valid expression data for", level))
          next
      }
      genes <- rownames(avg_expr)
      avg_expr <- as_tibble(avg_expr)
      avg_expr$Gene <- genes

      #tidy the data
      exp_lvl5 <- gather(avg_expr, cell_type, sumUMI, -Gene)
      exp_lvl5 <- exp_lvl5 %>% ungroup() %>% rename(Lvl5=cell_type, Expr_sum_mean=sumUMI)

      #remove not expressed genes
      genes_to_remove <- exp_lvl5 %>% group_by(Gene) %>% summarise(sum=sum(Expr_sum_mean)) %>% filter(sum==0)
      exp_lvl5 <- filter(exp_lvl5,!Gene%in%genes_to_remove$Gene)

      #Scale to 1 million molecules
      exp_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>% mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))

      #Specificity Calculation
#      exp_lvl5 <- exp_lvl5 %>% group_by(Gene) %>% mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>% select(Gene, Lvl5, specificity) %>% spread(Lvl5, specificity)

#      write.csv(exp_lvl5, paste0(prefix, "_", level, "_specificity_matrix.csv"), row.names = FALSE)
      write.csv(exp_lvl5, paste0("./", downdir, "/", prefix, "_", level, "_SumMean.csv"), row.names = FALSE)

      results$spec[[level]] <- exp_lvl5
    }
    return(results)
  } else {
    avg_expr <- AggregateExpression(sub_seurat, group.by = group_var, assays = "RNA")$RNA
    if (is.null(avg_expr) || nrow(avg_expr) == 0) {
        message(paste("No valid expression data for", group_var))
        return(NULL)
    }
    genes <- rownames(avg_expr)
    avg_expr <- as_tibble(avg_expr)
    avg_expr$Gene <- genes

    #tidy the data
    exp_lvl5 <- gather(avg_expr, cell_type, sumUMI, -Gene)
    exp_lvl5 <- exp_lvl5 %>% ungroup() %>% rename(Lvl5=cell_type, Expr_sum_mean=sumUMI)

    #remove not expressed genes
    genes_to_remove <- exp_lvl5 %>% group_by(Gene) %>% summarise(sum=sum(Expr_sum_mean)) %>% filter(sum==0)
    exp_lvl5 <- filter(exp_lvl5,!Gene%in%genes_to_remove$Gene)

    #Scale to 1 million molecules
    exp_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>% mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))

    #Specificity Calculation
#    exp_lvl5 <- exp_lvl5 %>% group_by(Gene) %>% mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>% select(Gene, Lvl5, specificity) %>% spread(Lvl5, specificity)

#    write.csv(exp_lvl5, paste0(prefix, "_specificity_matrix.csv"), row.names = FALSE)
    write.csv(exp_lvl5, paste0("./", downdir, "/", prefix, "_SumMean.csv"), row.names = FALSE)

    return(list(spec = exp_lvl5))
  }
}

# 1. lineage level
lineage_results <- calculate_expression_and_specificity(seurat_obj_Min, "lineage", "lineage", "downsampleMin")

# 2. sub_lineage level
sublineage_results <- calculate_expression_and_specificity(seurat_sub_Min, "sub_lineage", "sub_lineage", "downsampleMin")


# 4. group by age
#seurat_obj$age <- meta$stage
age_lineage_results <- calculate_expression_and_specificity(seurat_obj_500, c("stage", "lineage"), "stage_lineage", "downsample500")
age_sublineage_results <- calculate_expression_and_specificity(seurat_sub_500, c("stage", "sub_lineage"), "stage_sub_lineage", "downsample500")