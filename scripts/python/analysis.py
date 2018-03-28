import rna
import numpy as np
import pickle
import matplotlib.pyplot as plt

import importlib
importlib.reload(rna)

# Removed genes
all_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/removed_genes.pkl", 'rb'))
gene_of_interest = "Gata2"

# Distributions
rna.GeneExpression.plotDistribution(
    all_experiments, gene_of_interest, "violin")
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
#gene_of_interest = "Gata2"
gene_of_interest = ["Gata2"]

# Pattern
for i in range(0, len(reduced_experiments)):
    rna.plotGenesPatterns(
        all_experiments[i], reduced_experiments[i], gene_of_interest)
    plt.savefig("pattern." + str(i) + ".pdf")
    plt.close()


# Genes Comparison
# rna.plotGenesComparison(
#    all_experiments[0], reduced_experiments[0], gene_of_interest)


#tmp = rna.GeneExpression.copyFrom(all_experiments[1])
#tmp.expression = tmp.expression.loc[gene_of_interest]
#rna.plotCorrelationMatrix(tmp, gene_of_interest)


# Clustering
labels = []
for experiment in reduced_experiments:
    k = 3
    label = rna.cluster(experiment, k=k)
    labels.append(label)

    rna.GeneExpression.plot2DLabels(experiment, label)
    plt.savefig("k_" + str(k) + "." +
                experiment.condition + "_" + experiment.sample + ".pdf")
    plt.close()

rna.plotCorrelationClusters(all_experiments[0], labels[0],
                            all_experiments[1], labels[1])


# matr, _ = rna.confusionMatrixUsingGenes(all_experiments[0], labels[0],
#                                         all_experiments[1], labels[1])
