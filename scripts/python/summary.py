import rna
import numpy as np

# plot
import matplotlib.pyplot as plt

print("Start loading data")
input_data = "/Users/owner/Box Sync/UW/research/scRna/data/cellranger_10x/train"
all_experiments = rna.GeneExpression.loadFrom(input_data)
print("End loading data")

# Distribution of 0 cells and genes
thresholds = np.linspace(0, 1, 100)
n_0_cells = []
n_0_genes = []
for experiment in all_experiments:
    print(experiment.condition + "." + experiment.sample)
    for threshold in thresholds:
        cells = experiment.getToRemove(what_remove="cells",
                                       percentage_of_0=threshold,
                                       sign=">")
        n_0_cells.append(len(cells))

        genes = experiment.getToRemove(what_remove="genes",
                                       percentage_of_0=threshold,
                                       sign=">")
        n_0_genes.append(len(genes))

    # Detect threshold to remove data
    ind_cell = np.where(n_0_cells - np.max(n_0_cells) < 0)[0][0]
    threshold_cell = thresholds[ind_cell]

    ind_gene = np.where(n_0_genes - np.max(n_0_genes) < 0)[0][0]
    threshold_gene = thresholds[ind_gene]

    # plot
    plt.subplot(2, 1, 1)
    plt.plot(thresholds, n_0_cells, color="black")
    plt.axvline(x=threshold_cell, color="red", linestyle='--')
    plt.title("cells")
    plt.ylabel("count")
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',     # ticks along the bottom edge are off
        labelbottom='off'
    )

    plt.subplot(2, 1, 2)
    plt.plot(thresholds, n_0_genes, color="black")
    plt.axvline(x=threshold_gene, color="red", linestyle='--')
    plt.title("genes")
    plt.ylabel("count")
    plt.xlabel("percentage of 0")

    plt.suptitle(experiment.condition + "." + experiment.sample)
    plt.savefig("0_of_" + experiment.condition +
                "." + experiment.sample + ".pdf")
    n_0 = np.array([n_0_cells, n_0_genes])
    np.save(experiment.condition + "." + experiment.sample, n_0)
