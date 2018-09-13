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
if len(gene_of_interest) < 2:
    print("we have to have at least 2 genes")
    exit(1)


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
