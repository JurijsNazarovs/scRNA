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
plots_path = "/Users/owner/Box Sync/UW/research/scRna/plots_new/"
if not os.path.exists(plots_path):
    os.makedirs(plots_path)

# Input data
gene_of_interest = ["Gata2", "Socs2"]

all_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/removed_genes.pkl", 'rb'))
reduced_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/norm_log_pca_tsne.pkl", 'rb'))


# Distributions
for gene in gene_of_interest:
    rna.GeneExpression.plotDistribution(
        all_experiments, gene, "violin", orient="h")
    plt.savefig(plots_path + "expression_violin_" + gene + ".pdf")
    plt.close()

# Genes Pattern
for i in range(0, len(reduced_experiments)):
    for j in range(1, len(gene_of_interest)):
        print(i, j)
        rna.plotGenesPatterns(
            all_experiments[i], reduced_experiments[i],
            [gene_of_interest[0], gene_of_interest[j]],
            is_normalized=True)
        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        plt.savefig(plots_path + "pattern." + all_experiments[i].name +
                    gene_of_interest[0] + "_" + gene_of_interest[j] + ".pdf")
        plt.close()

        # Genes Comparison
        rna.plotGenesComparison(all_experiments[i], reduced_experiments[i],
                                [gene_of_interest[0], gene_of_interest[j]])

        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        plt.savefig(plots_path + "comparison." + all_experiments[i].name +
                    gene_of_interest[0] + "_" + gene_of_interest[j] + ".pdf")
        plt.close()

# Gene comparison with threshold
for i in range(0, len(reduced_experiments)):
    for j in range(1, len(gene_of_interest)):
        rna.plotGenesComparisonWithThreshold(all_experiments[i],
                                             reduced_experiments[i],
                                             [gene_of_interest[0],
                                              gene_of_interest[j]])

        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        plt.savefig(plots_path + "comparison_threshold." + all_experiments[i].name +
                    gene_of_interest[0] + "_" + gene_of_interest[j] + ".pdf")
        plt.close()


# Select k clusters for k-mean based on silhoette
for experiment in reduced_experiments:
    info = []
    print("----------------")
    print(experiment.name)
    ks = []
    for k in range(1, min(30, experiment.expression.shape[0])):
        clusterer = KMeans(n_clusters=k, random_state=10)
        labels = clusterer.fit_predict(experiment.expression)
        if len(set(labels)) > 1:
            ks.append(k)
            silhouette_avg = silhouette_score(
                experiment.expression, labels)
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


# Drop data to R to do a pheatmap
clusterer = KMeans(n_clusters=5, random_state=10)
labels = clusterer.fit_predict(reduced_experiments[0].expression)
rna.GeneExpression.plot2DLabels(reduced_experiments[0], labels)
plt.savefig(plots_path + "k = 5.pdf")
plt.close()


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
