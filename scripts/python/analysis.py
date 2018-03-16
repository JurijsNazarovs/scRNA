import rna
import numpy as np
import pickle
import matplotlib.pyplot as plt

import importlib
importlib.reload(rna)

# Removed genes
all_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/proceededData/removed_genes.pkl", 'rb'))
gene_of_interest = ""

# Distributions
# rna.GeneExpression.plotDistribution(
#     all_experiments, gene_of_interest, "violin")
# rna.GeneExpression.plotDistribution(
#     all_experiments, gene_of_interest, "hist")

# Correlation
# rna.GeneExpression.plotCorrelation(
#     all_experiments, gene_of_interest, "heatmap")
# rna.GeneExpression.plotCorrelation(
#     all_experiments, gene_of_interest, "scatter")


# # Reduced dimension data
reduced_experiments = pickle.load(open(
    "/Users/owner/Box Sync/UW/research/scRna/proceededData/tsne_2.pkl", 'rb'))
#gene_of_interest = "Gata2"
gene_of_interest = ["Grb10", "Gata2"]

# Pattern
for i in range(0, len(reduced_experiments)):
    rna.plotGenesPatterns(
        all_experiments[i], reduced_experiments[i], gene_of_interest)
    plt.savefig("pattern." + str(i) + ".pdf")
    plt.close()

# Clustering
# labels = []
# for experiment in reduced_experiments:
#     k = 3
#     label = (rna.cluster(experiment, k=k))
#     labels.append(label)
#     rna.plot2DLabels(experiment, label)
#     plt.savefig("k_" + str(k) + "." +
#                 experiment.condition + "_" + experiment.sample + ".pdf")
#     plt.close()

# for i in range(0, len(labels) - 1):
#     for j in range(i + 1, len(labels)):
#         conf_matr, gene_set = rna.confusionMatrixUsingGenes(all_experiments[i],
#                                                             labels[i],
#                                                             all_experiments[j],
#                                                             labels[j])
#         # import pdb
#         # pdb.set_trace()

#         name = all_experiments[i].condition + \
#             "." + all_experiments[i].sample +\
#             "_" + all_experiments[j].condition + \
#             "." + all_experiments[j].sample

#         # f = open("conf_matr." + name + ".txt", 'w')
#         # f.write(conf_matr)
#         # f.close()

#         # plot all clusters genes
#         k = 3
#         # plt.subplot(k, k)
#         iii = 1
#         fig = plt.figure()
#         for ii in range(0, k):
#             for jj in range(0, k):
#                 # import pdb
#                 # pdb.set_trace()

#                 plt.subplot(k, k, iii)
#                 gene_set_tmp = gene_set.get((ii, jj))
#                 tmp = rna.GeneExpression.copyFrom(all_experiments[i])
#                 tmp.expression = tmp.expression.loc[gene_set_tmp[0]]
#                 rna.plotCorrelationMatrix(tmp, gene_set_tmp[1])

#                 iii += 1

#         plt.savefig("cluster_" + name + ".1.pdf")
#         plt.close()
