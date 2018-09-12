import rna
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import importlib
importlib.reload(rna)
import sys
import os


def make_plot(gene_of_interest,
              all_experiments,
              reduced_experiments,
              is_normalized=False,
              plots_path="./plots_hot/"):
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)

    # cmap
    cmap = plt.get_cmap("hot")  # here just to define cmap.N
    all_counts = []
    for i in range(0, len(all_experiments)):
        original_expression = all_experiments[i]
        cells = original_expression.expression.loc[gene_of_interest] > 0
        counts = np.array(
            original_expression.expression.loc[gene_of_interest, cells].tolist())
        if is_normalized:
            counts = np.round(np.power(10, counts)) - 1
        all_counts.extend(counts)

    norm = mpl.colors.BoundaryNorm(
        np.arange(0.5, max(all_counts) + 1, 1), cmap.N)

    # Genes Pattern
    #fig = plt.figure(figsize=(10, 8))
    fig, axes = plt.subplots(figsize=(10, 8),
                             nrows=int(np.ceil(len(reduced_experiments) / 2)),
                             ncols=2)

    for i, ax in enumerate(axes.flat):
        print(i)
        plt.sca(ax)
        rna.plotGenesPatterns(
            all_experiments[i], reduced_experiments[i],
            gene_of_interest,
            is_normalized=is_normalized,
            cmap_norm=[norm],
            is_legend=False)

        # plt.axis('off')
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off
        plt.yticks([])

        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        # plt.show()

    # Legend
    max_n_on_cbar = 10
    #plt.subplots(1, 1)
    cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.02])
    if (max(all_counts) - 1 < max_n_on_cbar):
        plt.colorbar(ticks=np.arange(1, max(all_counts) + 1, 1),
                     orientation="horizontal", cax=cbar_ax)
    else:
        plt.colorbar(ticks=np.arange(
            1, max(all_counts) + 1, max(all_counts) // max_n_on_cbar),
            orientation="horizontal",
            cax=cbar_ax)

    # plt.show()
    plt.savefig(plots_path + "pattern_" + gene_of_interest + ".pdf")
    plt.close()


# Input data
if __name__ == "__main__":
    gene_of_interest = sys.argv[1]
    all_experiments = pickle.load(open(sys.argv[2], 'rb'))
    reduced_experiments = pickle.load(open(sys.argv[3], 'rb'))
    plots_path = sys.argv[4]
    if len(sys.argv) > 5:
        is_normalized = sys.argv[5]
    else:
        is_normalized = False

    # gene_of_interest = "Gata2"
    # all_experiments = pickle.load(open(
    #     "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/removed_genes.pkl", 'rb'))
    # reduced_experiments = pickle.load(open(
    #     "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/norm_log_pca_tsne.pkl", 'rb'))

    make_plot(gene_of_interest,
              all_experiments,
              reduced_experiments,
              plots_path=plots_path,
              is_normalized=is_normalized)
