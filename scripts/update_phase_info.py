import sys
import pandas as pd
import numpy as np


input_df = sys.argv[1]
phase_df = sys.argv[2]
output_df = sys.argv[3]


df = pd.read_csv(input_df, sep='\t')
phase_df = pd.read_csv(phase_df, sep='\t')
df = pd.merge(df, phase_df, on='Gene',  how="left")
df.to_csv(output_df, sep='\t', index=False, na_rep='NaN')
