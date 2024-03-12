#!/usr/bin/env python

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path



stat_file = sys.argv[1]
output_dir = sys.argv[2]
genes_without_benchmark = sys.argv[3]


def main():
    # Average coverage
    sample_name = ".".join(stat_file.split("/")[-1].split(".")[:-1])
    df = pd.read_csv(stat_file, sep='\t')
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(ncols=1, figsize=(16,9),)
    _ = sns.barplot(
            data=df,
            x='Gene',
            y='Average_cov',
            order=df.sort_values('Average_cov').Gene,
            color='blue'
        )
    axs.tick_params(axis='x', rotation=90)
    plt.xlabel("Gene", size=20, fontdict={'fontweight': 'bold'})
    plt.ylabel("Log Average coverage", size=20,fontdict={'fontweight': 'bold',})
    plt.tick_params(axis='x', labelsize=5, color='black', labelcolor='black')
    plt.tick_params(axis='y', labelsize=15, labelcolor='Blue')
    plt.yscale('log')
    plt.axhline(20, c='black', linestyle='--')
    plt.text(0, 21, "20x coverage", fontdict={'fontweight': 'bold', 'size': 12, 'color': '#13008B'})
    plt.title(sample_name, size=20, fontdict={'fontweight': 'bold', 'color': '#215F8F'})
    plt.tight_layout()
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    plt.savefig(os.path.join(output_dir, sample_name + "_average_gene_coverage.pdf"),
               bbox_inches = 'tight',
               dpi=300,
               )

    # F1-score
    genes_df = pd.read_csv(genes_without_benchmark, sep='\t').iloc[:, 0].to_list()
    df_273 = df.loc[~df.Gene.isin(genes_df)]
    sns.set_context(context='paper')
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(ncols=1, figsize=(16,9),)
    _ = sns.barplot(
            data=df_273,
            x='Gene',
            y='F-measure',
            order=df_273.sort_values('F-measure').Gene,
            color='blue'
        )
    axs.tick_params(axis='x', rotation=90)
    plt.xlabel("Gene", size=20, fontdict={'fontweight': 'bold'})
    plt.ylabel("F1-score", size=20,fontdict={'fontweight': 'bold',})
    plt.tick_params(axis='x', labelsize=5, color='black', labelcolor='black')
    plt.tick_params(axis='y', labelsize=15, labelcolor='Blue')
    plt.axhline(.5, c='black', linestyle='--')
    plt.text(0, 0.55, "50% F1-score", fontdict={'fontweight': 'bold', 'size': 12, 'color': '#13008B'})
    plt.title(sample_name, size=20, fontdict={'fontweight': 'bold', 'color': '#215F8F'})
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, sample_name + "_f1_score.pdf"),
               bbox_inches = 'tight',
               dpi=300,
               )

    # Less than 8x coverage
    df['lt8_per'] = df.lt_8_cover/df.Len
    sns.set_context(context='paper')
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(ncols=1, figsize=(16,9),)
    _ = sns.barplot(
            data=df,
            x='Gene',
            y=1 - df.lt8_per,
            order=df.sort_values('lt8_per', ascending=False).Gene,
            color='blue',
        alpha=.7
        )
    axs.tick_params(axis='x', rotation=90)
    plt.xlabel("Gene", size=20, fontdict={'fontweight': 'bold'})
    plt.ylabel("Percentage of gene with >8x base coverage", size=20,fontdict={'fontweight': 'bold',})
    plt.tick_params(axis='x', labelsize=5, color='black', labelcolor='black')
    plt.tick_params(axis='y', labelsize=15, labelcolor='Blue')
    plt.axhline(0.9, c='red', linestyle='--')
    plt.text(0, 0.92, "90% of Gene bases", fontdict={'fontweight': 'bold', 'size': 12, 'color': '#13008B'})
    plt.axhline(0.5, c='black', linestyle='--')
    plt.text(0, 0.52, "50% of Gene bases", fontdict={'fontweight': 'bold', 'size': 12, 'color': '#13008B'})
    plt.title(sample_name+' >8x bases%', size=20, fontdict={'fontweight': 'bold', 'color': '#215F8F'})
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, sample_name + "_lt8_genes_all.pdf"),
               bbox_inches = 'tight',
               dpi=300,
               )

if __name__ == "__main__":
    main()
