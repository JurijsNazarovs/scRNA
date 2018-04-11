import rna
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from itertools import chain

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import importlib
importlib.reload(rna)

import seaborn as sns
plots_path = "/Users/owner/Box Sync/UW/research/scRna/plots_analysis/"
if not os.path.exists(plots_path):
    os.makedirs(plots_path)
path_result = "/Users/owner/Box Sync/UW/research/scRna/data_analysis/"
if not os.path.exists(path_result):
    os.makedirs(path_result)

all_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/removed_genes.pkl", 'rb'))
reduced_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/norm_log_pca_tsne.pkl", 'rb'))


# ------ Below is a work with selected DE genes ------
# Select genes with distribution as Gata2
import scipy.stats as scst

all_genes = all_experiments[0].expression.index
res_p_val = []
for gene in all_genes:
    tmp_p_val = []
    for pair in [[0, 1], [2, 3], [0, 2], [0, 3], [1, 2], [1, 3]]:
        #[A, B], [F, I], [A,F], [A, I], [B, F], [B, I]
        tmp = (scst.ranksums(all_experiments[pair[0]].expression.loc[gene, :],
                             all_experiments[pair[1]].expression.loc[gene, :]))
        tmp_p_val.append(tmp[1])

    res_p_val.append(tmp_p_val)

res_p_val = pd.DataFrame(res_p_val)
res_p_val.index = all_genes

rep_cut = min(res_p_val.loc["Gata2", 0:2])
trt_cut = max(res_p_val.loc["Gata2", 2:7])
final_genes = all_genes[(res_p_val.loc[:, 0:1].min(axis=1) > rep_cut) &
                        (res_p_val.loc[:, 2:6].max(axis=1) <= trt_cut)]
print("Number of selected genes:" + str(len(final_genes)))

edit_experiments = []
for experiment in all_experiments:
    experiment_tmp = rna.GeneExpression.copyFrom(experiment)
    experiment_tmp.expression = experiment_tmp.expression.loc[final_genes, :]
    edit_experiments.append(experiment_tmp)

pickle.dump(edit_experiments,  open(path_result + "norm_log_DE.pkl", "wb"))


# Clustering: Select k based on silhoette
for experiment in edit_experiments:
    info = []
    print("----------------")
    print(experiment.name)
    ks = []
    for k in range(1, min(20, len(final_genes))):
        clusterer = KMeans(n_clusters=k, random_state=10)
        labels = clusterer.fit_predict(experiment.expression.T)
        if len(set(labels)) > 1:
            ks.append(k)
            silhouette_avg = silhouette_score(
                experiment.expression.T, labels)
            print("for k = %d: %.5f" % (k, silhouette_avg))
            info.append(silhouette_avg)
    print("----------------")
    plt.plot(ks, info, label=experiment.name)

plt.legend()
plt.ylabel("Average silhouette")
plt.xlabel("k")
plt.savefig(plots_path + "k_means_avg_sil.pdf")
plt.close()


importlib.reload(rna)
# selected_k = [3, 3, 6, 6]
selected_k = [4, 4, 6, 4]
labels_all = {}
for i in range(0, len(selected_k)):
    experiment = edit_experiments[i]
    clusterer = KMeans(n_clusters=selected_k[i], random_state=10)
    labels = clusterer.fit_predict(experiment.expression.T)
    labels_all[experiment.name] = labels

    rna.plotClusterExpression(experiment, labels)
    plt.savefig(plots_path + "cluster_" + experiment.name + ".png")
    plt.close()

# Write in csv for R
edit_experiments[0].expression.to_csv("single_expression.csv",
                                      index=False,
                                      header=True)


# Other heatmaps variation

# 1)  row by cells (where all 4 sample  are concatenated) - allow clustering of
# the genes and cells but add a column label with 4 distinct colors so that we can
# see how samples are getting mixed up (see below for an example)
combined_experiments = rna.GeneExpression.combine(edit_experiments)
col_labels = [i.name for i in edit_experiments]
col_labels = list(chain.from_iterable(
    [i.name] * i.expression.shape[1] for i in edit_experiments))
rna.plotClusterExpression(combined_experiments, col_labels, row_cluster=True)
plt.savefig(plots_path + "combined_rows-cluster.png")
plt.close()
# 2)  row by samples (where for each gene we only have 4*4 values that depict its
# (25% quantile, median, mean, 75th quantile) expressionfor each sample) - allow
# clustering of the genes, but not the columns.

info = {}
for experiment in edit_experiments:
    data_tmp = experiment.expression.transpose()
    info_tmp = data_tmp.describe()
    info[experiment.name] = info_tmp.iloc[[1, 4, 5, 6], :].T

info = pd.concat(info, axis=1)
info.columns = [info.columns.levels[1][i] for i in info.columns.labels[1]]
col_labels = list(chain.from_iterable(
    [[i.name] * 4 for i in edit_experiments]))
combined_experiments_tmp = rna.GeneExpression(info)
rna.plotClusterExpression(combined_experiments_tmp, col_labels, row_cluster=True,
                          xticks=info.columns)
plt.savefig(plots_path + "combined_rows-cluster_4stat.png")
plt.close()


# 3)  row by cell for sample A and then take that row ordering and plot the other
# samples with that ordering by allowing clustering of their cells (this will be 4
# different plots where each has the same rows).
clustered_genes = sns.clustermap(
    edit_experiments[0].expression, row_cluster=True, col_cluster=F).data2d.index

# Heatmap
for experiment in edit_experiments:
    sns.clustermap(experiment.expression.loc[clustered_genes],
                   col_cluster=True, row_cluster=False,
                   xticklabels=False, yticklabels=False)
    plt.savefig(plots_path + "fixed_rows_" + experiment.name + ".png")
    plt.close()

# Another heatmap
n_colors = plt.cm.gist_heat.N
colors1 = np.flip(plt.cm.gist_heat(np.linspace(0, 1, n_colors)), axis=0)
range1 = np.linspace(0, 0.3, n_colors)
colors2 = "b"
range2 = 0.6
colors3 = "g"
range3 = 1

colorsT = [i for i in zip(range1, colors)]
colorsT.append((range2, colors2))
colorsT.append((range3, colors3))
cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colorsT)


class TwoInnerPointsNormalize(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, low=None, up=None, clip=False):
        self.low = low
        self.up = up
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.low, self.up, self.vmax], [0, 0.3, 0.8, 1]
        return np.ma.masked_array(np.interp(value, x, y))


norm = TwoInnerPointsNormalize(vmin=0, vmax=4, low=0.6, up=1.2)
for experiment in edit_experiments:
    sns.clustermap(experiment.expression.loc[clustered_genes],
                   col_cluster=True, row_cluster=False,
                   xticklabels=False, yticklabels=False,
                   cmap=cmap, norm=norm)
    plt.savefig(plots_path + "fixed_rows_" +
                experiment.name + "_another_map.png")
    plt.close()
