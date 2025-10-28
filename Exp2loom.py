import os, sys
import argparse
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
parser = argparse.ArgumentParser(description="Transform the expression file to loom file")
parser.add_argument('--celltype', required=True, help="Path to expression matrix")
args = parser.parse_args()

in_mtx = args.celltype + "_mtx.csv"
x=sc.read_csv(in_mtx);
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
out_loom = args.celltype + ".loom"
lp.create(out_loom, x.X.transpose(), row_attrs, col_attrs);
