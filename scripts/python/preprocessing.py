import rna
import numpy as np
import pickle

import importlib
importlib.reload(rna)

print("Start loading data")
input_data = "/Users/owner/Box Sync/UW/research/scRna/data/cellranger_10x/train"
all_experiments = rna.GeneExpression.loadFrom(input_data)
print("End loading data")


# Remove genes
thresholds_genes = [0.98989899, 0.98989899, 0.98989899, 0.98989899]
for i in range(0, len(all_experiments)):
    experiment = all_experiments[i]
    genes_remove = experiment.getToRemove(what_remove="genes",
                                          percentage_of_0=thresholds_genes[i],
                                          sign=">")
    genes_remove = rna.removeElementFromList(genes_remove, "Gata2")
    experiment.remove(genes_remove, what_remove="genes")

# Save data with removed genes
pickle.dump(all_experiments, open(
    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/removed_genes.pkl", 'wb'))


for experiment in all_experiments:
    experiment.reduceDimension(p=2, method="tsne")

pickle.dump(all_experiments, open(
    "/Users/owner/Box Sync/UW/research/scRna/proceeded_data/tsne_2.pkl", 'wb'))


for experiment in all_experiments:
    print(experiment.expression.shape)
