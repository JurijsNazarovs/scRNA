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
plots_path = "plots/"

print("Start loading data")
input_data = "/Users/owner/Box Sync/UW/research/scRna/data/cellranger_10x/train"
all_experiments = rna.GeneExpression.loadFrom(input_data)
# all_experiments = pickle.load(open(
#    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/removed_genes.pkl", 'rb'))
print("End loading data")

# Distribution of 0 cells and genes
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

# Other statistics
# what_compare = "cells"
# distributions = {}
# fig, ax_hist = plt.subplots(nrows=3, ncols=2,
#                             sharex=True, sharey=True)
# ax_hist = ax_hist.flatten()
# for i in range(0, len(all_experiments)):
#     experiment = all_experiments[i]
#     if what_compare.lower() == "genes":
#         axis = 1
#     else:
#         axis = 0

#     ind_of_0 = (experiment.expression == 0).sum(
#         axis=axis) / experiment.expression.shape[axis]
#     distributions[experiment.name] = (ind_of_0)

#     # Histogram
#     ax = ax_hist[i]
#     ax.hist(ind_of_0, bins=20)
#     ax.set_title(experiment.name)
#     pass

# ax_hist[-1].axis('off')
# fig.text(0.01, 0.5, 'Number of ' + what_compare,
#          va='center', rotation='vertical')
# plt.suptitle("Percentage of 0 in " + what_compare)
# plt.xlim([0, 1])
# # plt.show()
# plt.savefig("comparison_of_" + what_compare + ".pdf")
# plt.close()

# Number of genes which are 0 everywhere
# Number of cells which are 0 everywhere
# Histogram of cells based on number of 0 per cell => get distribution
# histogram of genes based on number of 0 per gene => get distribtuion
# Compare these distributions among each experiment
# Create subplots with overlapped distributions
# Create boxplots and violins

# distributions = pd.DataFrame(distributions)
# distributions.describe()

# Depth per cell
what_compare = "cells"

# for i in range(0, len(all_experiments)):
for i in range(0, 5):
    experiment = all_experiments[i]
    ind_of_0 = (experiment.expression == 0).sum(axis=0)\
        / experiment.expression.shape[0]
    sns.distplot(ind_of_0, kde=False, label=experiment.name)

plt.title("Depth per cell")
plt.xlabel("ratio of zero expressed genes")
plt.ylabel("number of cells")
plt.legend()
#plt.xlim((0.8, 1))
# plt.show()
plt.savefig(plots_path + "depth_cell.pdf")
plt.close()


for i in range(0, 2):
    experiment = all_experiments[i]
    ind_of_0 = (experiment.expression == 0).sum(axis=0)\
        / experiment.expression.shape[0]
    sns.distplot(ind_of_0, kde=False, label=experiment.sample)

plt.title("Depth per cell - " + experiment.condition)
plt.xlabel("ratio of zero expressed genes")
plt.ylabel("number of cells")
plt.legend()
#plt.xlim((0.8, 1))
# plt.show()
plt.savefig(plots_path + "depth_cell_wt.pdf")
plt.close()

for i in range(2, 4):
    experiment = all_experiments[i]
    ind_of_0 = (experiment.expression == 0).sum(axis=0)\
        / experiment.expression.shape[0]
    sns.distplot(ind_of_0, kde=False, label=experiment.sample)

plt.title("Depth per cell - " + experiment.condition)
plt.xlabel("ratio of zero expressed genes")
plt.ylabel("number of cells")
plt.legend()
#plt.xlim((0.8, 1))
# plt.show()
plt.savefig(plots_path + "depth_cell_mutant.pdf")
plt.close()

# Active genes - hist
for i in range(0, 5):
    experiment = all_experiments[i]
    sum_genes = experiment.expression.sum(axis=0)
    sum_total = sum_genes.sum(axis=0)
    sum_genes /= sum_total
    sns.distplot(sum_genes, kde=False, label=experiment.name)

plt.title("Percentage of expressed genes wrt to total expression")
plt.xlabel("gene expression")
plt.ylabel("number of cells")
plt.legend()
# plt.show()
plt.savefig(plots_path + "gene_expression.pdf")
plt.close()


# Sequencing depth per cell - violin and boxplot
sum_genes = {}
for i in range(0, 5):
    experiment = all_experiments[i]
    sum_genes[experiment.name] = experiment.expression.sum(axis=0)
sum_genes = pd.DataFrame(sum_genes)
sns.violinplot(data=sum_genes)
plt.title("Sequencing depth")
plt.ylabel("total depth per cell")
plt.legend()
# plt.show()
plt.savefig(plots_path + "sequence_depth_violin.pdf")
plt.close()

sns.boxplot(data=sum_genes)
plt.title("Sequencing depth")
plt.ylabel("total depth per cell")
plt.legend()
# plt.show()
plt.savefig(plots_path + "sequencing_depth_boxplot.pdf")
plt.close()


# # Heatmap of percentage of genes that have min read in # of cells
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
