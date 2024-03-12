#!/usr/bin/env python
import sys
import pandas as pd

if len(sys.argv) < 4:
    print(f"The command should be line:\n{sys.argv[0]} inputDf InputGenesDf OuputDfName")
    exit(1)

input_df = sys.argv[1]
input_uncovered_genes_df = sys.argv[2]
output_df = sys.argv[3]

input_df = pd.read_csv(input_df, sep = '\t')
genes_df = pd.read_csv(input_uncovered_genes_df, sep = '\t',  names=['Gene', 'Unocvered variants'])

new_df = pd.merge(input_df, genes_df, on='Gene', how='left').fillna(0, downcast='infer')

new_df.to_csv(output_df, sep='\t', index=False, na_rep='NaN')
