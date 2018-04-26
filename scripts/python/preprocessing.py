import rna
import numpy as np
import pickle
import os

import importlib
importlib.reload(rna)

path_input = "/Users/owner/Box Sync/UW/research/scRna/data/cellranger_10x/train/"
path_result = "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/"
if not os.path.exists(path_result):
    os.makedirs(path_result)

print("Start loading data")
all_experiments = rna.GeneExpression.loadFrom(path_input)
for experiment in all_experiments:
    print(experiment.name)
print("End loading data")

# Remove genes with min.reads = 2 in at least 2 cells
print("Remove genes which do not have at least 2 genes in at least 2 cells")
genes_remove = []
for experiment in all_experiments:
    genes_remove_tmp = experiment.getToRemoveMinCount(what_remove="genes",
                                                      min_reads=2,
                                                      n_cells=2)
    print(len(genes_remove_tmp) / experiment.expression.shape[0] * 100)
    genes_remove.append(genes_remove_tmp)

# Get intersection of removing genes
genes_remove = rna.intersect(genes_remove)
for experiment in all_experiments:
    experiment.remove(genes_remove, what_remove="genes")
    print(experiment.expression.shape[0])
pickle.dump(all_experiments, open(
    path_result + "removed_genes.pkl", 'wb'))

for experiment in all_experiments:
    experiment.expression.to_csv(path_result + experiment.name + "_rg.csv",
                                 index=True,
                                 header=True)

# Median normalization
rna.normalize(all_experiments)
# Log normalization
for experiment in all_experiments:
    experiment.expression = np.log10(1 + experiment.expression)
pickle.dump(all_experiments, open(
    path_result + "removed_genes_med_log.pkl", 'wb'))

for experiment in all_experiments:
    experiment.expression.to_csv(path_result + experiment.name + "_rg_med_log.csv",
                                 index=True,
                                 header=True)

# Dimension reduction for combined data
combined_experiments = rna.GeneExpression.combine(all_experiments)

print("pca")
combined_experiments.reduceDimension(p=6000, method="tsvd")
combined_experiments.expression = combined_experiments.expression.T
pickle.dump(combined_experiments, open(path_result + "pca_combined.pkl", 'wb'))

print("t-sne")
combined_experiments.reduceDimension(p=2, method="tsne")
pickle.dump(combined_experiments, open(
    path_result + "pca+tsne_combined.pkl", 'wb'))

# Splitting reduced dimension experiments
reduced_experiments = []
prev_step = 0
for experiment in all_experiments:
    new_step = prev_step + experiment.expression.shape[1]
    tmp = rna.GeneExpression.copyFrom(experiment)
    tmp.expression = combined_experiments.expression.iloc[prev_step:new_step, :]
    prev_step = new_step
    reduced_experiments.append(tmp)

# Check dimensions
for i in range(0, len(reduced_experiments)):
    print(reduced_experiments[i].expression.shape[0],
          all_experiments[i].expression.shape[1])

pickle.dump(reduced_experiments, open(
    path_result + "norm_log_pca_tsne.pkl", 'wb'))
