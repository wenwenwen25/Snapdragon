import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("WV_sce_exp.csv");#R中导出的表达矩阵
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("WV_sce.loom",x.X.transpose(),row_attrs,col_attrs)
