import sys
import pandas as pd
import scipy.io
from scipy.sparse import coo_matrix

barcodes_file = sys.argv[1] + "barcodes.tsv.gz"
features_file = sys.argv[1] + "features.tsv.gz"
matrix_file = sys.argv[1] + "matrix.mtx.gz"
meta_file = sys.argv[1] + "meta_subtype.tsv"

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
mean_expression_by_cell_type = merged_data.groupby('lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean.csv')
# calculate the mean expression of every gene in each <strong>sub-cell type</strong>
mean_expression_by_cell_type = merged_data.groupby('sub_lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_SubCellType.csv')

# calculate the mean expression of every gene in each cell type in Male and Female
mean_expression_by_cell_type = merged_data.loc[merged_data['sex']=='Male'].groupby('lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_Male.csv')
mean_expression_by_cell_type = merged_data.loc[merged_data['sex']=='Female'].groupby('lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_Female.csv')

# calculate the mean expression of every gene in each sub-cell type in Male and Female
mean_expression_by_cell_type = merged_data.loc[merged_data['sex']=='Male'].groupby('sub_lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_SubCellType_Male.csv')
mean_expression_by_cell_type = merged_data.loc[merged_data['sex']=='Female'].groupby('sub_lineage').mean()
mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
mean_expression.T.to_csv('Mean_SubCellType_Female.csv')

for age in ["trimester2nd", "trimester3rd", "years0_1", "years1_2", "years2_4", "years4_10", "years10_20", "Adult"]:
    mean_expression_by_cell_type = merged_data.loc[merged_data['age_range']==age].groupby('lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_' + age + '.csv')
    mean_expression_by_cell_type = merged_data.loc[merged_data['age_range']==age].groupby('sub_lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_SubCellType_' + age + '.csv')
    
    mean_expression_by_cell_type = merged_data.loc[(merged_data['age_range']==age) & (merged_data['sex']=='Male')].groupby('lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_Male_' + age + '.csv')
    mean_expression_by_cell_type = merged_data.loc[(merged_data['age_range']==age) & (merged_data['sex']=='Female')].groupby('lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_Female_' + age + '.csv')
    
    mean_expression_by_cell_type = merged_data.loc[(merged_data['age_range']==age) & (merged_data['sex']=='Male')].groupby('sub_lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_SubCellType_Male_' + age + '.csv')
    mean_expression_by_cell_type = merged_data.loc[(merged_data['age_range']==age) & (merged_data['sex']=='Female')].groupby('sub_lineage').mean()
    mean_expression = mean_expression_by_cell_type.loc[:,'FAM138A':'MT-CYB']
    mean_expression.T.to_csv('Mean_SubCellType_Female_' + age + '.csv')
