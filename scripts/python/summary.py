#import rna
#from rna import GeneExpression as ge
import rna
import numpy as np
import pandas as pd
import pickle

import importlib
importlib.reload(rna)

# plot
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette("muted")
plots_path = "plots_analysis/"

# Load data
print("Start loading data")
# Option 1 - raw data files
input_data = "/Users/owner/Box Sync/UW/research/scRna/data/cellranger_10x/train"
all_experiments = rna.GeneExpression.loadFrom(input_data)
# rna.pickleDumpLargeFile(all_experiments,
#                        "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/raw_data.pkl")

# Option 2 - raw data pickle - does not work
# all_experiments = pickle.load(open(
#   "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/raw_data.pkl", 'rb'))

# Option 3 - removed genes pickle
# all_experiments = pickle.load(open(
#    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/removed_genes.pkl", 'rb'))
print("End loading data")

for experiment in all_experiments:
    experiment.info()

# Distribution of 0 cells and genes - plot with a red dotted vertical line
# thresholds = np.linspace(0, 1, 100)
# n_0_cells = []
# n_0_genes = []
# for experiment in all_experiments:
#     print(experiment.condition + "." + experiment.sample)
#     for threshold in thresholds:
#         cells = experiment.getToRemove(what_remove="cells",
#                                        percentage_of_0=threshold,
#                                        sign=">")
#         n_0_cells.append(len(cells))

#         genes = experiment.getToRemove(what_remove="genes",
#                                        percentage_of_0=threshold,
#                                        sign=">")
#         n_0_genes.append(len(genes))

#     # Detect threshold to remove data
#     ind_cell = np.where(n_0_cells - np.max(n_0_cells) < 0)[0][0]
#     threshold_cell = thresholds[ind_cell]

#     ind_gene = np.where(n_0_genes - np.max(n_0_genes) < 0)[0][0]
#     threshold_gene = thresholds[ind_gene]

#     # plot
#     plt.subplot(2, 1, 1)
#     plt.plot(thresholds, n_0_cells, color="black")
#     plt.axvline(x=threshold_cell, color="red", linestyle='--')
#     plt.title("cells")
#     plt.ylabel("count")
#     plt.tick_params(
#         axis='x',          # changes apply to the x-axis
#         which='both',      # both major and minor ticks are affected
#         bottom='off',     # ticks along the bottom edge are off
#         labelbottom='off'
#     )

#     plt.subplot(2, 1, 2)
#     plt.plot(thresholds, n_0_genes, color="black")
#     plt.axvline(x=threshold_gene, color="red", linestyle='--')
#     plt.title("genes")
#     plt.ylabel("count")
#     plt.xlabel("percentage of 0")

#     plt.suptitle(experiment.condition + "." + experiment.sample)
#     plt.savefig("0_of_" + experiment.condition +
#                 "." + experiment.sample + ".pdf")
#     n_0 = np.array([n_0_cells, n_0_genes])
#     np.save(experiment.condition + "." + experiment.sample, n_0)

# Depth per cell - all experiments
what_compare = "cells"
seq_depth = []
for i in range(0, len(all_experiments)):
    experiment = all_experiments[i]
    ind_of_0 = (experiment.expression == 0).sum(axis=0)\
        / experiment.expression.shape[0]
    sns.distplot(ind_of_0, kde=False, label=experiment.name)
    seq_depth.append(ind_of_0)

plt.xlabel("Proportion of genes with 0 expression")
plt.ylabel("Number of cells")
plt.legend()
plt.savefig(plots_path + "depth_cell.pdf")
plt.close()

seq_depth = pd.DataFrame(seq_depth).transpose()
seq_depth.describe()
for i in seq_depth.columns:
    print("%d. %.5f" % (i, np.nanmedian(seq_depth[i])))

# Depth per cell - mutant
for i in range(0, 2):
    experiment = all_experiments[i]
    ind_of_0 = (experiment.expression == 0).sum(axis=0)\
        / experiment.expression.shape[0]
    sns.distplot(ind_of_0, kde=False, label=experiment.sample)

plt.xlabel("Proportion of genes with 0 expression")
plt.ylabel("Number of cells")
plt.legend()
plt.savefig(plots_path + "depth_cell_mutant.pdf")
plt.close()

# Depth per cell - wild type
for i in range(2, 4):
    experiment = all_experiments[i]
    ind_of_0 = (experiment.expression == 0).sum(axis=0)\
        / experiment.expression.shape[0]
    sns.distplot(ind_of_0, kde=False, label=experiment.sample)

plt.xlabel("Proportion of genes with 0 expression")
plt.ylabel("Number of cells")
plt.legend()
plt.savefig(plots_path + "depth_cell_wt.pdf")
plt.close()

# Distribution of total expression per cell
seq_depth = []
for i in range(0, 5):
    experiment = all_experiments[i]
    sum_genes = experiment.expression.sum(axis=0)
    sum_total = sum_genes.sum(axis=0)
    sum_genes /= sum_total
    sns.distplot(sum_genes, kde=False, label=experiment.name)
    seq_depth.append(sum_genes)

plt.xlabel("Gene expression proportion")
plt.ylabel("Number of cells")
plt.legend()
plt.savefig(plots_path + "gene_expression.pdf")
plt.close()

seq_depth = pd.DataFrame(seq_depth).transpose()
seq_depth.describe()
for i in seq_depth.columns:
    print("%d. %.5f" % (i, np.nanmedian(seq_depth[i])))


# Sequencing depth per cell - violin and boxplot
sum_genes = {}
for i in range(0, 5):
    experiment = all_experiments[i]
    sum_genes[experiment.name] = experiment.expression.sum(axis=0)
sum_genes = pd.DataFrame(sum_genes)
sns.violinplot(data=sum_genes)
plt.ylabel("Number of reads per cell")
plt.legend()
plt.savefig(plots_path + "sequence_depth_violin.pdf")
plt.close()

sns.boxplot(data=sum_genes)
plt.ylabel("Number of reads per cell")
plt.legend()
plt.savefig(plots_path + "sequencing_depth_boxplot.pdf")
plt.close()

# Heatmap of proportion of genes that have min read in # of cells
n_cells = list(range(1, 21))
min_reads = list(range(1, 21))
for i in range(0, 5):
    experiment = all_experiments[i]
    genes_distr = []
    for j in min_reads:
        good_genes = (experiment.expression >= j).sum(axis=1)
        genes_distr.append(good_genes.value_counts().sort_index()[n_cells]
                           / experiment.expression.shape[0])
    genes_distr = pd.concat(genes_distr, axis=1)
    genes_distr.columns = min_reads
    sns.heatmap(genes_distr,
                cmap='seismic', vmin=0, vmax=genes_distr.max().max(),
                xticklabels=5, yticklabels=5)
    # plt.yticks(n_cells)
    plt.title(experiment.name)
    plt.xlabel("minimum total count per gene")
    plt.ylabel("number of cells")
    plt.savefig(plots_path + "atleast_counts_" + experiment.name + ".pdf")
    plt.close()
    print(i)

# Heatmap of proportion of genes that have min read in min # of cells
n_cells = list(range(1, 21))
min_reads = list(range(1, 21))
for i in range(0, len(all_experiments)):
    experiment = all_experiments[i]
    perc_bad_genes = []
    for j in min_reads:
        n_bad_genes_per_count = []
        good_genes = (experiment.expression >= j).sum(axis=1)
        good_genes = good_genes.value_counts().sort_index()
        perc_bad_genes_tmp = [
            100 - good_genes[k:].sum() / experiment.expression.shape[0] * 100
            for k in n_cells]
        perc_bad_genes.append(perc_bad_genes_tmp)

    rem_genes_perc = pd.DataFrame(perc_bad_genes).T
    rem_genes_perc.columns = min_reads
    rem_genes_perc.index = n_cells
    sns.heatmap(rem_genes_perc,
                cmap='gnuplot2', vmin=0, vmax=100,
                xticklabels=5, yticklabels=5,
                cbar_kws={'label': 'Percentage of genes to remove'})
    # plt.yticks(n_cells)
    plt.title(experiment.name)
    plt.xlabel("Minimum count per cell")
    plt.ylabel("Minimum number of cells")
    plt.savefig(plots_path + "atleast_counts_cells_" +
                experiment.name + ".pdf")
    plt.close()
    print(i)


# Heatmap of proportion of genes that have min read in # of cells
n_cells = list(range(1, 21))
min_reads = list(range(1, 21))
for i in range(0, len(all_experiments)):
    experiment = all_experiments[i]
    genes_distr = []
    for j in min_reads:
        good_genes = (experiment.expression >= j).sum(axis=1)
        genes_distr.append(good_genes.value_counts().sort_index()[n_cells]
                           / experiment.expression.shape[0])
    genes_distr = pd.concat(genes_distr, axis=1)
    genes_distr.columns = min_reads
    sns.heatmap(genes_distr,
                cmap='seismic', vmin=0, vmax=genes_distr.max().max(),
                xticklabels=5, yticklabels=5)
    # plt.yticks(n_cells)
    plt.title(experiment.name)
    plt.xlabel("minimum total count per gene")
    plt.ylabel("number of cells")
    plt.savefig(plots_path + "atleast_counts_" + experiment.name + ".pdf")
    plt.close()
    print(i)


# Violin of percentiles - primary done for filtered genes
info = {}
for experiment in all_experiments:
    data_tmp = experiment.expression.transpose()
    info_tmp = data_tmp.describe()
    info[experiment.name] = info_tmp.iloc[4:7, :]

for stat_name in ["25%", "50%", "75%"]:
    data_plot = [info[name].loc[stat_name, :] for name in info.keys()]
    data_plot = pd.DataFrame(data_plot).transpose()
    data_plot.columns = info.keys()
    sns.violinplot(data=data_plot)
    plt.ylabel(stat_name + " quantile of gene expression")
    plt.ylim([-0.15, 2.15])
    plt.savefig(plots_path + "quantile_" + stat_name + ".pdf")
    plt.close()
