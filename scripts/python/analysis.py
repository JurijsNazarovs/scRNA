import rna
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import importlib
importlib.reload(rna)

import seaborn as sns
plots_path = "plots_analysis/"

# Removed genes
all_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/removed_genes.pkl", 'rb'))
gene_of_interest = ["Gata2", "Exosc3", "Samd14", "Kit", "Hdc", "Grb10"]

# Distributions
for gene in gene_of_interest:
    rna.GeneExpression.plotDistribution(
        all_experiments, gene, "violin", orient="h")
    plt.savefig(plots_path + "expression_violin_" + gene + ".pdf")
    plt.close()
# rna.GeneExpression.plotDistribution(
#     all_experiments, gene_of_interest, "hist")

# Correlation
# rna.GeneExpression.plotCorrelation(
#     all_experiments, gene_of_interest, "heatmap")
# rna.GeneExpression.plotCorrelation(
#     all_experiments, gene_of_interest, "scatter")


# # Reduced dimension data
reduced_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/tsne_2.pkl", 'rb'))
# gene_of_interest = "Gata2"
gene_of_interest = ["Gata2"]

# Genes Pattern
# for i in range(0, len(reduced_experiments)):
#     rna.plotGenesPatterns(
#         all_experiments[i], reduced_experiments[i], gene_of_interest)
#     plt.savefig(plots_path + "pattern." + str(i) + ".pdf")
#     plt.close()


# Genes Comparison
# rna.plotGenesComparison(
#    all_experiments[0], reduced_experiments[0], gene_of_interest)


# tmp = rna.GeneExpression.copyFrom(all_experiments[1])
# tmp.expression = tmp.expression.loc[gene_of_interest]
# rna.plotCorrelationMatrix(tmp, gene_of_interest)


# Clustering
labels = []
for experiment in reduced_experiments:
    k = 3
    label = rna.cluster(experiment, k=k)
    labels.append(label)

    rna.GeneExpression.plot2DLabels(experiment, label)
    plt.savefig(plots_path + "k_" + str(k) + "." +
                experiment.condition + "_" + experiment.sample + ".pdf")
    plt.close()

# rna.plotCorrelationClusters(all_experiments[0], labels[0],
#                             all_experiments[1], labels[1])


# matr, _ = rna.confusionMatrixUsingGenes(all_experiments[0], labels[0],
#                                         all_experiments[1], labels[1])


# Working with cells. Decrese dimension and clustering cells.
# PCA
from sklearn.decomposition import TruncatedSVD

# Plot explained variance ratio
p = 10
svd = TruncatedSVD(n_components=p, n_iter=10)
for experiment in all_experiments:
    svd.fit(experiment.expression.T)
    plt.plot(np.array(range(0, p)) + 1,
             svd.explained_variance_ratio_, label=experiment.name)

plt.legend()
plt.ylabel("Explained variance ratio")
plt.xlabel("Component")
plt.savefig(plots_path + "pca_marginal_variance.pdf")
plt.close()

# plot sum of explained variance ratio
for experiment in all_experiments:
    svd.fit(experiment.expression.T)
    explained_var = [sum(svd.explained_variance_ratio_[0: i + 1])
                     for i in range(0, len(svd.explained_variance_ratio_))]
    plt.plot(np.array(range(0, p)) + 1,
             explained_var, label=experiment.name)

plt.legend()
plt.ylabel("Cumulative explained variance ratio")
plt.xlabel("Component")
plt.savefig(plots_path + "pca_cum_variance.pdf")
plt.close()


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

edit_experiments = []
for experiment in all_experiments:
    experiment_tmp = rna.GeneExpression.copyFrom(experiment)
    experiment_tmp.expression = experiment_tmp.expression.loc[final_genes, :]
    edit_experiments.append(experiment_tmp)


# Select k based on silhoette
for experiment in edit_experiments:
    info = []
    print("----------------")
    print(experiment.name)
    ks = []
    for k in range(1, min(11, len(final_genes))):
        clusterer = KMeans(n_clusters=k, random_state=10)
        labels = clusterer.fit_predict(experiment.expression.T)
        if len(set(labels)) > 1:
            ks.append(k)
            silhouette_avg = silhouette_score(
                experiment.expression.T, labels)
            print("for k = %d: %.5f" % (k, silhouette_avg))
            info.append(silhouette_avg)
    print("----------------")
    plt.plot(ks,
             info, label=experiment.name)

plt.legend()
plt.ylabel("Average silhouette")
plt.xlabel("k")
plt.savefig(plots_path + "k_means_avg_sil.pdf")
plt.close()


importlib.reload(rna)
#selected_k = [3, 3, 6, 6]
selected_k = [7, 7, 4, 4]
for i in range(0, len(selected_k)):
    experiment = edit_experiments[i]
    clusterer = KMeans(n_clusters=selected_k[i], random_state=10)
    labels = clusterer.fit_predict(experiment.expression.T)

    rna.plotClusterExpression(experiment, labels)
    plt.savefig(plots_path + "cluster_" + experiment.name + ".png")
    plt.close()
