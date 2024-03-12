import sys
import pandas as pd
import numpy as np


input_df = sys.argv[1]
output_df = sys.argv[2]
genes_without_benchmark = sys.argv[3]

df = pd.read_csv(input_df, sep='\t')
genes_df = pd.read_csv(genes_without_benchmark, sep='\t').iloc[:, 0].to_list()
df.loc[df.Gene.isin(genes_df), ['Precision', 'Sensitivity', 'F-measure']] = "NaN"
df.to_csv(output_df, sep='\t', index=False, na_rep='NaN')
