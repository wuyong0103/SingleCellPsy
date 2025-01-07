import sys
import gc
import pandas as pd
import scipy.io
from scipy.sparse import coo_matrix

barcodes_file = sys.argv[1] + "barcodes.tsv.gz"
features_file = sys.argv[1] + "features.tsv.gz"
matrix_file = sys.argv[1] + "matrix.mtx.gz"
meta_file = sys.argv[1] + "meta_subtype_right.tsv"

#read the matrix expression file
barcodes = pd.read_csv(barcodes_file, header=None, names = ['Cell_ID'])
features = pd.read_csv(features_file, header=None, names = ["Gene"])
matrix = scipy.io.mmread(matrix_file)
exp = matrix.toarray()
exp = pd.DataFrame(exp, index=features['Gene'], columns=barcodes['Cell_ID'])

#read the meta file
meta = pd.read_csv(meta_file, sep='\t', index_col=0)

#merge the meta and gene expression file
merged_data = pd.concat([meta, exp.T], axis=1, join='inner')

# calculate the mean expression of every gene in each cell type
type_counts = merged_data['lineage'].value_counts()
valid_types = type_counts[type_counts >= 20].index
filtered_df_lineage = merged_data[merged_data['lineage'].isin(valid_types)]
mean_expression_by_cell_type = filtered_df_lineage.iloc[:, [4] + list(range(6, filtered_df_lineage.shape[1]))].groupby('lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean.csv')

# calculate the mean expression of every gene in each <strong>sub-cell type</strong>
type_counts = merged_data['sub_lineage'].value_counts()
valid_types = type_counts[type_counts >= 20].index
filtered_df_sub = merged_data[merged_data['sub_lineage'].isin(valid_types)]
mean_expression_by_cell_type = filtered_df_sub.iloc[:,5:].groupby('sub_lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_SubCellType.csv')

# calculate the mean expression of every gene in each cell type in Male and Female
type_age_counts = merged_data.groupby(['lineage', 'sex']).size()
valid_combinations = type_age_counts[type_age_counts >= 20].index
filtered_df_lineagesex = merged_data[merged_data.set_index(['lineage', 'sex']).index.isin(valid_combinations)]
mean_expression_by_cell_type = filtered_df_lineagesex.loc[filtered_df_lineagesex['sex']=='Male'].iloc[:, [4] + list(range(6, filtered_df_lineagesex.shape[1]))].groupby('lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_Male.csv')
mean_expression_by_cell_type = filtered_df_lineagesex.loc[filtered_df_lineagesex['sex']=='Female'].iloc[:, [4] + list(range(6, filtered_df_lineagesex.shape[1]))].groupby('lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_Female.csv')

# calculate the mean expression of every gene in each sub-cell type in Male and Female
type_age_counts = merged_data.groupby(['sub_lineage', 'sex']).size()
valid_combinations = type_age_counts[type_age_counts >= 20].index
filtered_df_subsex = merged_data[merged_data.set_index(['sub_lineage', 'sex']).index.isin(valid_combinations)]
mean_expression_by_cell_type = filtered_df_subsex.loc[filtered_df_subsex['sex']=='Male'].iloc[:,5:].groupby('sub_lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_SubCellType_Male.csv')
mean_expression_by_cell_type = filtered_df_subsex.loc[filtered_df_subsex['sex']=='Female'].iloc[:,5:].groupby('sub_lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_SubCellType_Female.csv')

del filtered_df_lineage
del filtered_df_sub
del filtered_df_lineagesex
del filtered_df_subsex
gc.collect()

type_age_counts = merged_data.groupby(['lineage', 'age_range']).size()
valid_combinations = type_age_counts[type_age_counts >= 20].index
filtered_df_lineageage = merged_data[merged_data.set_index(['lineage', 'age_range']).index.isin(valid_combinations)]

type_age_counts = merged_data.groupby(['sub_lineage', 'age_range']).size()
valid_combinations = type_age_counts[type_age_counts >= 20].index
filtered_df_subage = merged_data[merged_data.set_index(['sub_lineage', 'age_range']).index.isin(valid_combinations)]

type_sex_age_counts = merged_data.groupby(['lineage', 'sex', 'age_range']).size()
valid_combinations = type_sex_age_counts[type_sex_age_counts >= 20].index
filtered_df_lineagesexage = merged_data[merged_data.set_index(['lineage', 'sex', 'age_range']).index.isin(valid_combinations)]

type_sex_age_counts = merged_data.groupby(['sub_lineage', 'sex', 'age_range']).size()
valid_combinations = type_sex_age_counts[type_sex_age_counts >= 20].index
filtered_df_subsexage = merged_data[merged_data.set_index(['sub_lineage', 'sex', 'age_range']).index.isin(valid_combinations)]

for age in ["trimester2nd", "trimester3rd", "years0_1", "years1_2", "years2_4", "years4_10", "years10_20", "Adult"]:
    mean_expression_by_cell_type = filtered_df_lineageage.loc[filtered_df_lineageage['age_range']==age].iloc[:, [4] + list(range(6, filtered_df_lineageage.shape[1]))].groupby('lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_' + age + '.csv')
    mean_expression_by_cell_type = filtered_df_subage.loc[filtered_df_subage['age_range']==age].iloc[:,5:].groupby('sub_lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_SubCellType_' + age + '.csv')

del filtered_df_lineageage
del filtered_df_subage
gc.collect()

for age in ["trimester2nd", "trimester3rd", "years0_1", "years1_2", "years2_4", "years4_10", "years10_20", "Adult"]:    
    mean_expression_by_cell_type = filtered_df_lineagesexage.loc[(filtered_df_lineagesexage['age_range']==age) & (filtered_df_lineagesexage['sex']=='Male')].iloc[:, [4] + list(range(6, filtered_df_lineagesexage.shape[1]))].groupby('lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_Male_' + age + '.csv')
    mean_expression_by_cell_type = filtered_df_lineagesexage.loc[(filtered_df_lineagesexage['age_range']==age) & (filtered_df_lineagesexage['sex']=='Female')].iloc[:, [4] + list(range(6, filtered_df_lineagesexage.shape[1]))].groupby('lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_Female_' + age + '.csv')
    
    mean_expression_by_cell_type = filtered_df_subsexage.loc[(filtered_df_subsexage['age_range']==age) & (filtered_df_subsexage['sex']=='Male')].iloc[:,5:].groupby('sub_lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_SubCellType_Male_' + age + '.csv')
    mean_expression_by_cell_type = filtered_df_subsexage.loc[(filtered_df_subsexage['age_range']==age) & (filtered_df_subsexage['sex']=='Female')].iloc[:,5:].groupby('sub_lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_SubCellType_Female_' + age + '.csv')
