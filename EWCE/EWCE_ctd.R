library(EWCE)
library(Matrix)
barcodes_path <- "../barcodes.tsv.gz"
features_path <- "../features.tsv.gz"
matrix_path <- "../matrix.mtx.gz"
barcodes <- readLines(gzfile(barcodes_path))
features <- read.delim(gzfile(features_path), header = FALSE, stringsAsFactors = FALSE)
expression_matrix <- readMM(gzfile(matrix_path))
rownames(expression_matrix) <- features$V1
colnames(expression_matrix) <- barcodes
dim(expression_matrix)
meta <- read.table("../meta_subtype_EWCE.tsv", sep="\t", header=T)
rownames(meta) <- meta$cell
dim(meta)
exp <- as.matrix(expression_matrix)
exp <- fix_bad_hgnc_symbols(exp)

m <- match(colnames(exp), meta$cell)
meta <- meta[m,]
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(exp = exp, input_species = "human", output_species = "human", level2annot = meta$sub_lineage, no_cores=80)
annotLevels <- list(level1class=meta$lineage,level2class=meta$sub_lineage)
ctd <- EWCE::generate_celltype_data(exp = exp_CortexOnly_DROPPED, annotLevels = annotLevels, groupName = "Cortex_All", savePath=getwd(), no_cores=80)

meta_male <- meta[meta$sex=="Male",]
m <- match(meta_male$cell, colnames(exp))
exp_male <- exp[,m]
exp_male_dropped <- EWCE::drop_uninformative_genes(exp = exp_male, input_species = "human", output_species = "human", level2annot = meta_male$sub_lineage, no_cores=30)
annotLevels_male <- list(level1class=meta_male$lineage, level2class=meta_male$sub_lineage)
ctd_male <- EWCE::generate_celltype_data(exp = exp_male_dropped, annotLevels = annotLevels_male, groupName = "Cortex_male", savePath=getwd())

meta_female <- meta[meta$sex=="Female",]
m <- match(meta_female$cell, colnames(exp))
exp_female <- exp[,m]
exp_female_dropped <- EWCE::drop_uninformative_genes(exp = exp_female, input_species = "human", output_species = "human", level2annot = meta_female$sub_lineage, no_cores=30)
annotLevels_female <- list(level1class=meta_female$lineage, level2class=meta_female$sub_lineage)
ctd_female <- EWCE::generate_celltype_data(exp = exp_female_dropped, annotLevels = annotLevels_female, groupName = "Cortex_female", savePath=getwd())

for (age in c("trimester2nd", "trimester3rd", "years0_1", "years1_2", "years2_4", "years4_10", "years10_20", "Adult")) {
	meta_age <- meta[meta$age_range==age,]
	m <- match(meta_age$cell, colnames(exp))
	exp_age <- exp[,m]
	exp_age_dropped <- EWCE::drop_uninformative_genes(exp = exp_age, input_species = "human", output_species = "human", level2annot = meta_age$sub_lineage, no_cores=30)
	annotLevels_age <- list(level1class=meta_age$lineage, level2class=meta_age$sub_lineage)
	gn <- paste("Cortex_", age, sep = "")
	ctd_age <- EWCE::generate_celltype_data(exp = exp_age_dropped, annotLevels = annotLevels_age, groupName = gn, savePath=getwd())
	for (sex in c("Male", "Female")){
		meta_sex_age <- meta[meta$sex==sex & meta$age_range==age, ]
		m <- match(meta_sex_age$cell, colnames(exp))
		exp_sex_age <- exp[,m]
		exp_sex_age_dropped <- EWCE::drop_uninformative_genes(exp = exp_sex_age, input_species = "human", output_species = "human", level2annot = meta_sex_age$sub_lineage, no_cores=30)
		annotLevels_sex_age <- list(level1class=meta_sex_age$lineage, level2class=meta_sex_age$sub_lineage)
		gn <- paste("Cortex_", age, sex, sep = "")
		ctd_sex_age <- EWCE::generate_celltype_data(exp = exp_sex_age_dropped, annotLevels = annotLevels_sex_age, groupName = gn, savePath=getwd())
	}
}
