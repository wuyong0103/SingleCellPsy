library(Seurat)
library(Matrix)
library(readr)

# Read the expression matrix
expr_mat <- read_csv("matrix.csv") 


expr_df <- as.data.frame(expr_mat)
rownames(expr_df) <- expr_df[[1]] 
expr_df <- expr_df[, -1] 

#Transpose matrix
expr_df_t <- t(expr_df)

#Transform to sparse matrix
sparse_mat <- Matrix(expr_df_t, sparse = TRUE)

#Create Seurat object
Bakken2021 <- CreateSeuratObject(counts = sparse_mat, project = "Bakken2021")

#Save object
saveRDS(Bakken2021, file = "Bakken2021.rds")

meta <- read_csv("metadata.csv") %>% 
    select(sample_name, external_donor_name_label) %>% 
    rename(barcode = sample_name, donorid = external_donor_name_label) %>% 
    mutate(sex = ifelse(donorid == "H18.30.001", "F", "M"), age = ifelse(donorid =="H18.30.001", "60y", "50y"))
write.csv(meta, "meta_Bakken2021.csv")
