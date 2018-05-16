import pandas as pd
import numpy as np
import scipy as sc
import scipy.io
import os
import csv
import glob
import copy
import pickle
import sys

# stat
from scipy.stats.stats import pearsonr
from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans

# plot
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# import importlib
# importlib.reload(rna)


class GeneExpression(object):
    def __init__(self,
                 expression_matrix,
                 gene_names=None,
                 cell_names=None,
                 condition="condition",
                 sample="sample"):
        if isinstance(expression_matrix, pd.DataFrame):
            self.expression = expression_matrix
        else:
            self.expression = pd.DataFrame(expression_matrix)

        if not gene_names is None:
            self.expression.index = gene_names
        if not cell_names is None:
            self.expression.columns = cell_names

        self.condition = condition
        self.sample = sample
        self.name = self.condition + "." + self.sample

    # Other constructors
    @classmethod
    def combine(cls, list_gene_expression):
        expression = [i.expression for i in list_gene_expression]
        condition = [i.condition for i in list_gene_expression]
        sample = [i.sample for i in list_gene_expression]

        expression = pd.concat(expression, axis=1)
        condition = "_".join(set(condition))
        sample = "_".join(set(sample))

        new_expression = cls(expression_matrix=expression,
                             condition=condition,
                             sample=sample)
        return(new_expression)

    @classmethod
    def copyFrom(cls, gene_expression):
        new_expression = copy.deepcopy(gene_expression)
        return(new_expression)

    # Class methods
    @classmethod
    def loadFrom(cls,
                 dir_path,
                 expression_postfix="matrix.mtx",
                 genes_list_postfix="genes.tsv",
                 cells_list_postfix="barcodes.tsv"):
        # dir_path contains directories for every sample
        # every directory inside is cosidered as sample
        # every directory might have a file type.wt, where wt goes to
        # condition of gene expression object
        gene_expressions = []  # resulted list
        for sample in next(os.walk(dir_path))[1]:  # directory
            # Detect pathes
            sample_dir_path = dir_path + "/" + sample
            genes_path = sample_dir_path + "/" + genes_list_postfix
            cells_path = sample_dir_path + "/" + cells_list_postfix
            expression_path = sample_dir_path + "/" + expression_postfix

            # Read files to create a GeneExpression object
            gene_names = [row[1]
                          for row in csv.reader(open(genes_path), delimiter="\t")]
            cell_names = [row[0]
                          for row in csv.reader(open(cells_path), delimiter="\t")]
            expression_mat = scipy.io.mmread(expression_path).toarray()
            condition = glob.glob(sample_dir_path + "/type.*")
            condition = condition[0].split(sample_dir_path + "/type.")[1] if\
                len(condition) > 0 else ""
            # condition = condition.split(
            #    sample_dir_path + "/type.")[1]

            # Create a GeneExpression object
            gene_expression = GeneExpression(expression_matrix=expression_mat,
                                             gene_names=gene_names,
                                             cell_names=cell_names,
                                             condition=condition,
                                             sample=sample)
            # Filter
            gene_expression.collapse()
            gene_expressions.append(gene_expression)

        return(gene_expressions)

    # Instance Methods
    def info(self):
        print("Sample: " + self.sample)
        print("Condition: " + self.condition)
        print("Unique name: " + self.name)
        print("Number fo genes: %d" % (self.expression.shape[0]))
        print("Number of cells: %d" % (self.expression.shape[1]))

    def remove(self, list_remove, what_remove="genes"):
        # Removes cells/genes from list_to_remove
        if what_remove.lower() == "cells":
            final_list_remove = intersect(
                [list_remove, self.expression.columns])
            axis = 1
        else:
            final_list_remove = intersect(
                [list_remove, self.expression.index])
            axis = 0

        if len(final_list_remove) != 0:
            print(self.condition + "." + self.sample + ": removing " +
                  str(len(final_list_remove) / self.expression.shape[axis] * 100) +
                  "% of " + what_remove)
            self.expression.drop(final_list_remove, axis=axis, inplace=True)

    def getToRemove(self, what_remove="cells", percentage_of_0=0.5, sign=">"):
        # Returns list of cells/genes which are expressed in *sign* then
        # *percentage_of_0* of genes/cells
        what_remove = what_remove.lower()
        if what_remove.lower() == "cells":
            axis = 0
            shape_ind = 0
        elif what_remove.lower() == "genes":
            axis = 1
            shape_ind = 1
        else:
            print("what_remove  = ", what_remove, " is not supported")
            return(None)

        # Detect indicies to remove
        ind_remove = (self.expression == 0).sum(axis=axis)

        if sign == "<":
            ind_remove = ind_remove < round(
                percentage_of_0 * self.expression.shape[shape_ind])
        elif sign == ">":
            ind_remove = ind_remove > round(
                percentage_of_0 * self.expression.shape[shape_ind])
        elif sign == "=":
            ind_remove = ind_remove == round(
                percentage_of_0 * self.expression.shape[shape_ind])
        else:
            print("sign ", sign, " is not supoprted")
            return(None)

        return(list(ind_remove[ind_remove].index))

    def getToRemoveMinCount(self, what_remove="genes", min_reads=2, n_cells=2):
        # Returns list of cells/genes which are expressed in *sign* then
        # *percentage_of_0* of genes/cells
        what_remove = what_remove.lower()
        if what_remove.lower() == "cells":
            axis = 0
            shape_ind = 0
        elif what_remove.lower() == "genes":
            axis = 1
            shape_ind = 1
        else:
            print("what_remove  = ", what_remove, " is not supported")
            return(None)

        # Detect indicies to remove
        ind_remove = (self.expression >= min_reads).sum(axis=axis)
        ind_remove = ~(ind_remove >= n_cells)

        return(list(ind_remove[ind_remove].index))

    def collapse(self, method="sum"):
        all_genes = list(self.expression.index)
        duplicated_genes = set(
            [i for i in all_genes if all_genes.count(i) > 1])
        collapsed_expression = []
        for gene in duplicated_genes:
            if method == "sum":
                collapsed_expression.append(self.expression.loc[gene].sum())
            else:
                print("method is not supported")
                return(Null)

        collapsed_expression = pd.DataFrame(collapsed_expression)
        collapsed_expression.index = duplicated_genes

        new_expression = self.expression.drop(duplicated_genes)
        new_expression.append(collapsed_expression)
        self.expression = new_expression

    def reduceDimension(self, p=2, method="tsne", what_reduce="rows"):
        if what_reduce == "rows":
            matrix = self.expression.T
        else:
            matrix = self.expression

        if method == "tsvd":
            svd = TruncatedSVD(n_components=p, n_iter=7)
            matrix_transform = svd.fit_transform(matrix)
        if method == "tsne":
            tsne = TSNE(n_components=p, perplexity=30, n_iter=300, verbose=1)
            matrix_transform = tsne.fit_transform(matrix)

        new_expression = pd.DataFrame(matrix_transform, index=matrix.index)
        self.expression = new_expression

    def correlation(self, gene):
        all_genes = list(self.expression.index)

        if not gene in self.expression.index:
            correlation = [np.nan for i in range(0, len(all_genes))]
        else:
            correlation = []
            for gene_tmp in all_genes:
                gene_expression = self.expression.loc[gene].values
                gene_tmp_expression = self.expression.loc[gene_tmp].values

                p_corr = pearsonr(gene_expression, gene_tmp_expression)[0]
                correlation.append(p_corr)

        correlation = pd.DataFrame(correlation, columns=[gene])
        correlation.index = all_genes
        return(correlation)

    # Plots
    @classmethod
    def plotDistribution(cls,
                         list_gene_expression,
                         gene,
                         plot_type=["violin"],
                         orient="v"):

        if isinstance(plot_type, str):
            plot_type = [plot_type]

        if isinstance(list_gene_expression, GeneExpression):  # does not work
            list_gene_expression = [list_gene_expression]

        uniq_id = [i.name for i in list_gene_expression]
        gene_expression = [i.expression.loc[gene]
                           for i in list_gene_expression]

        for plt_type in plot_type:
            if plt_type == "violin":
                ax = sns.violinplot(data=gene_expression, orient=orient)
                if orient == "v":
                    ax.set(xticklabels=uniq_id,
                           ylabel="count")
                else:
                    ax.set(yticklabels=uniq_id,
                           xlabel="count")
            elif plt_type == "boxplot":
                ax = sns.boxplot(data=gene_expression, orient=orient)
                if orient == "v":
                    ax.set(xticklabels=uniq_id,
                           ylabel="count")
                else:
                    ax.set(yticklabels=uniq_id,
                           xlabel="count")
            elif plt_type == "hist":
                for i in range(0, len(uniq_id)):
                    ge_tmp = gene_expression[i].values
                    step = np.diff(np.unique(ge_tmp)).min()
                    left_of_first_bin = ge_tmp.min() - float(step) / 2
                    right_of_last_bin = ge_tmp.max() + float(step) / 2
                    bins = np.arange(left_of_first_bin,
                                     right_of_last_bin, step)

                    plt.hist(ge_tmp,
                             bins=bins,
                             alpha=0.5,
                             label=uniq_id[i])
                    plt.xlabel("counts")
                    plt.ylabel("n of cells")

                plt.legend(loc='upper right')
            elif plt_type == "kdeplot":
                for i in range(0, len(uniq_id)):
                    sns.kdeplot(gene_expression[i].values,
                                label=uniq_id[i], shade=True)

                plt.legend(loc='upper right')
            elif plt_type == "distplot":
                for i in range(0, len(uniq_id)):
                    sns.distplot(gene_expression[i].values,
                                 label=uniq_id[i])

                plt.legend(loc='upper right')

            plt.title(gene)
            plt.show(block=False)

    @classmethod
    def plotCorrelation(cls,
                        list_gene_expression,
                        gene,
                        plot_type=["heatmap"]):
        # plots correlation of one gene with other for different expressions

        if isinstance(plot_type, str):
            plot_type = [plot_type]

        if isinstance(list_gene_expression, GeneExpression):  # does not work
            list_gene_expression = [list_gene_expression]

        uniq_id = [i.name for i in list_gene_expression]
        correlation_set = []
        for gene_expression in list_gene_expression:
            correlation_set.append(gene_expression.correlation(gene))

        correlation_set = pd.concat(correlation_set, axis=1)
        correlation_set.columns = uniq_id

        for plt_type in plot_type:
            if plt_type == "heatmap":
                with sns.axes_style('dark'):
                    sns.heatmap(correlation_set, yticklabels=False,
                                cmap='seismic', vmin=-1, vmax=1)
            elif plt_type == "scatter":
                for col in correlation_set.columns:
                    y = correlation_set.loc[:, col].values
                    x = np.arange(0, len(y))
                    sns.regplot(x, y, fit_reg=False, dropna=True,
                                label=col)
                    plt.legend()
        plt.title("Correlation of " + gene)
        plt.xticks()
        plt.show(block=False)

    @classmethod
    def plot2DLabels(cls, expression, labels):
        if expression.expression.shape[1] != 2 or \
                expression.expression.shape[0] == 0:
            print("expression should have shape (x !=0, 2)")
            return None

        cmap = plt.cm.jet
        norm = mpl.colors.BoundaryNorm(
            np.arange(labels.min() - 0.5, labels.max() + 1, 1), cmap.N)

        legend_color = []
        for label in set(labels):
            ind = labels == label
            plt.scatter(expression.expression.loc[ind, 0],
                        expression.expression.loc[ind, 1],
                        label=str(label + 1),
                        c=labels[ind],
                        cmap=cmap,
                        norm=norm,
                        edgecolors="black")
            legend_color.append(cmap(norm(label)))

        # Set legend color
        plt.legend(loc='upper right')
        ax = plt.gca()
        leg = ax.get_legend()
        for i in range(0, len(legend_color)):
            leg.legendHandles[i].set_color(legend_color[i])

        plt.title(expression.name + ":  k = " + str(labels.max() + 1))
        plt.show(block=False)


# Different methods not related to a class
def intersect(list_of_lists):
    if not isinstance(list_of_lists, list):
        print("Provided argumnet is not a list")
        return(None)
    new_list = set(list_of_lists[0])

    for i in range(1, len(list_of_lists)):
        new_list = new_list & set(list_of_lists[i])
    return(list(new_list))


def removeElementFromList(list_, element):
    if element in list_:
        list_ = [i for i in list_ if i != element]
    return list_


def plotGenesPatterns(original_expression, reduced_expression, genes,
                      is_normalized=True, cmap_norm=None, is_legend=True):

    # put it in plots
    cmaps = ["hot", 'Blues', 'Reds', 'Greens', 'Purples']
    if reduced_expression.expression.shape[1] != 2 or \
            reduced_expression.expression.shape[0] == 0:
        print("reduced_expression should have shape (x !=0, 2)")
        return None
    if isinstance(genes, str):
        genes = [genes]
    if len(genes) > len(cmaps):
        print("Just " + str(len(cmaps)) + " genes can be plotted")
        genes = genes[0:len(cmaps)]
        print("Following genes will be plotted:")
        print(genes)

    # Main plot
    plt.scatter(reduced_expression.expression.loc[:, 0],
                reduced_expression.expression.loc[:, 1],
                label="",
                facecolors='none',
                edgecolor='none')  # change for grey
    plt.title(original_expression.name)

    size_step = 100
    size = len(genes) * size_step + 10
    legend_color = []
    for i in range(0, len(genes)):
        if not genes[i] in original_expression.expression.index:
            print("No gene " + genes[i] + " is presented")
            continue
        cells = original_expression.expression.loc[genes[i]] > 0

        counts = np.array(
            original_expression.expression.loc[genes[i], cells].tolist())
        if is_normalized:
            counts = np.round(np.power(10, counts)) - 1

        cmap = plt.get_cmap(cmaps[i])
        colors_of_cmap = cmap(np.linspace(
            0.2, 1, cmap.N))  # del white part
        colors_of_cmap = np.flip(colors_of_cmap, axis=0)

        # Create a new colormap from those colors
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            None, colors_of_cmap)

        if cmap_norm == None:
            norm = mpl.colors.BoundaryNorm(
                np.arange(0.5, counts.max() + 1, 1), cmap.N)
        else:
            norm = cmap_norm[i]

        plt.scatter(reduced_expression.expression.loc[cells, 0],
                    reduced_expression.expression.loc[cells, 1],
                    c=counts,
                    norm=norm,
                    label="%s: %.2f" % (genes[i], counts.max()),
                    edgecolor="black",
                    cmap=cmap,
                    alpha=0.9,
                    s=size)
        size -= size_step

        if is_legend:
            legend_color.append(cmap(norm(cmap.N)))
            max_n_on_cbar = 10
            if len(genes) == 1:
                if (counts.max() - 1 < max_n_on_cbar):
                    plt.colorbar(ticks=np.arange(1, counts.max() + 1, 1),
                                 orientation="horizontal")
                else:
                    plt.colorbar(ticks=np.arange(
                        1, counts.max() + 1, counts.max() // max_n_on_cbar),
                        orientation="horizontal")

    # Set legend color
    if is_legend:
        plt.legend(loc='upper right')
        ax = plt.gca()
        leg = ax.get_legend()
        for i in range(0, len(legend_color)):
            leg.legendHandles[i].set_color(legend_color[i])

    # plt.show(block=False)


def plotGenesComparison(original_expression, reduced_expression, genes):
    # Plot genes and their overlap on decreased dimension
    colors = ["blue", "red", "yellow"]
    if reduced_expression.expression.shape[1] != 2 or \
            reduced_expression.expression.shape[0] == 0:
        print("reduced_expression should have shape (x !=0, 2)")
        return None
    if isinstance(genes, str):
        genes = [genes]
    if len(genes) != 2:
        print("Just 2 genes can be plotted")
        return None

    # Main plot
    plt.scatter(reduced_expression.expression.loc[:, 0],
                reduced_expression.expression.loc[:, 1],
                label="",
                facecolors='none',
                edgecolor='grey')
    plt.title(original_expression.name)

    # Every gene separately
    legend_color = []
    cells_set = []
    for i in range(0, len(genes)):
        if not genes[i] in original_expression.expression.index:
            print("No gene " + genes[i] + " is presented")
            return None

        cells = original_expression.expression.loc[genes[i]] > 0
        cells_set.append(cells.index[cells].tolist())

        plt.scatter(reduced_expression.expression.loc[cells, 0],
                    reduced_expression.expression.loc[cells, 1],
                    c=colors[i],
                    label=genes[i],
                    edgecolor="black",
                    alpha=0.9)

    # Genes overlap
    cells_set = intersect(cells_set)
    plt.scatter(reduced_expression.expression.loc[cells_set, 0],
                reduced_expression.expression.loc[cells_set, 1],
                c=colors[len(colors) - 1],
                label=" & ".join(genes),
                edgecolor="black",
                alpha=0.9)

    # Set legend color
    plt.legend(loc='upper right')
    plt.show(block=False)


def plotGenesComparisonWithThreshold(original_expression, reduced_expression, genes):
        # Plot genes and their overlap on decreased dimension
    colors = ["blue", "red", "yellow"]
    if reduced_expression.expression.shape[1] != 2 or \
            reduced_expression.expression.shape[0] == 0:
        print("reduced_expression should have shape (x !=0, 2)")
        return None
    if isinstance(genes, str):
        genes = [genes]
    if len(genes) != 2:
        print("Just 2 genes can be plotted")
        return None

    # Main plot
    plt.scatter(reduced_expression.expression.loc[:, 0],
                reduced_expression.expression.loc[:, 1],
                label="",
                facecolors='none',
                edgecolor='grey')
    plt.title(original_expression.name)

    # Every gene separately
    legend_color = []
    thresholds = []
    for i in range(0, len(genes)):
        if not genes[i] in original_expression.expression.index:
            print("No gene " + genes[i] + " is presented")
            return None
        thresholds.append(original_expression.expression.loc[genes[i]].mean())

    # Blue
    cells = (original_expression.expression.loc[genes[0]] >= thresholds[0]) &\
            (original_expression.expression.loc[genes[1]] < thresholds[1])

    plt.scatter(reduced_expression.expression.loc[cells, 0],
                reduced_expression.expression.loc[cells, 1],
                c=colors[0],
                label=("%s >= %.2f & %s < %.2f" %
                       (genes[0], thresholds[0], genes[1], thresholds[1])),
                edgecolor="black",
                alpha=0.9)

    # Red
    cells = (original_expression.expression.loc[genes[0]] <= thresholds[0]) &\
            (original_expression.expression.loc[genes[1]] > thresholds[1])

    plt.scatter(reduced_expression.expression.loc[cells, 0],
                reduced_expression.expression.loc[cells, 1],
                c=colors[1],
                label=("%s <= %.2f & %s > %.2f" %
                       (genes[0], thresholds[0], genes[1], thresholds[1])),
                edgecolor="black",
                alpha=0.9)

    # Yellow
    cells = (original_expression.expression.loc[genes[0]] >= thresholds[0]) &\
            (original_expression.expression.loc[genes[1]] > thresholds[1])

    plt.scatter(reduced_expression.expression.loc[cells, 0],
                reduced_expression.expression.loc[cells, 1],
                c=colors[2],
                label=("%s >= %.2f & %s > %.2f" %
                       (genes[0], thresholds[0], genes[1], thresholds[1])),
                edgecolor="black",
                alpha=0.9)

    # Set legend color
    plt.legend(loc='lower left')
    plt.show(block=False)


def normalize(experiments, method="median"):
    if method == "median":
        col_sums = [np.sum(exp.expression, axis=0) for exp in experiments]
        col_sums = [s for sublist in col_sums for s in sublist]

        median_sum = np.median(col_sums)
        for experiment in experiments:
            col_sum = np.sum(experiment.expression, axis=0)
            experiment.expression = experiment.expression * median_sum / col_sum


def cluster(experiment, k=2, method="kmeans"):
    if method == "kmeans":
        kmeans_model = KMeans(
            n_clusters=k, random_state=1).fit(experiment.expression)
        labels = kmeans_model.labels_

        return(labels)


def confusionMatrix(rep1, labels1, rep2, labels2):
    confusion_matrix = []

    for i in set(labels1):
        elements1 = rep1.expression.index[labels1 == i]
        row_tmp = []
        for j in set(labels2):
            elements2 = rep2.expression.index[labels2 == j]
            row_tmp.append(len(intersect([elements1, elements2])))

        confusion_matrix.append(row_tmp)
    confusion_matrix = pd.DataFrame(confusion_matrix,
                                    columns=set(labels2),
                                    index=set(labels1))
    return(confusion_matrix)


def confusionMatrixUsingGenes(rep1, labels1, rep2, labels2, threshold=0):
    # chose genes in cluster expression of which > threshold in every cell
    # inside the cluster
    confusion_matrix = []
    cluster_genes = {}

    for i in set(labels1):
        genes_id = (
            rep1.expression.loc[:, labels1 == i] > threshold).all(axis=1)
        genes1 = rep1.expression.index[genes_id]
        row_tmp = []
        for j in set(labels2):
            genes_id = (
                rep2.expression.loc[:, labels2 == j] > threshold).all(axis=1)
            genes2 = rep2.expression.index[genes_id]

            common_genes = intersect([genes1, genes2])
            cluster_genes[i, j] = [genes1.tolist(), genes2.tolist()]

            row_tmp.append(len(common_genes))
        confusion_matrix.append(row_tmp)
    confusion_matrix = pd.DataFrame(confusion_matrix,
                                    columns=set(labels2),
                                    index=set(labels1))
    return([confusion_matrix, cluster_genes])


def plotCorrelationMatrix(self, genes):
    if not isinstance(genes, list):
        genes = list(genes)
    correlation_set = []
    for gene in genes:
        correlation_set.append(self.correlation(gene))

    correlation_set = pd.concat(correlation_set, axis=1)
    correlation_set.columns = genes

    yticklabels, xticklabels = False, False
    if len(correlation_set.columns) < 30:
        yticklabels = True
    if len(correlation_set.index) < 30:
        xticklabels = True

    with sns.axes_style('dark'):
        sns.heatmap(correlation_set,
                    yticklabels=yticklabels, xticklabels=xticklabels,
                    cmap='seismic', vmin=-1, vmax=1)

    plt.title(self.name)
    plt.show(block=False)


def plotClusterExpression(self, cluster_labels,
                          row_cluster=False, col_clusters=False,
                          yticks=False, xticks=False):
    # Plots heatmap with labels for columns sorted together
    expression = pd.DataFrame.copy(self.expression)
    expression.columns = cluster_labels
    expression.sort_index(axis=1, inplace=True)

    lut = dict(zip(set(expression.columns),
                   sns.mpl_palette("hsv", len(set(expression.columns)))))
    col_colors = pd.DataFrame(expression.columns)[0].map(lut)

    expression.columns.name = "Cells clusters"
    expression.index.name = "Genes"
    sns_plot = sns.clustermap(expression, col_colors=col_colors.values,
                              col_cluster=col_clusters, row_cluster=row_cluster,
                              cmap="gnuplot2",
                              yticklabels=yticks, xticklabels=xticks)
    #sns_plot.ax_heatmap.set_xlabel = "Cells cluster"
    #sns_plot.ax_heatmap.set_ylabel = "Genes"

    # Add legend for clusters
    for label in sorted(set(cluster_labels)):
        sns_plot.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                       label=label, linewidth=0)
        sns_plot.ax_col_dendrogram.legend(loc="center", ncol=6)

    # Move color bar
    sns_plot.cax.set_position([.15, .2, .03, .45])


def plotCorrelationClusters(experiment_main, label_main,
                            experiment_relative, label_relative):
    _, gene_set = confusionMatrixUsingGenes(experiment_main,
                                            label_main,
                                            experiment_relative,
                                            label_relative)
    # main - what we use to do correlation
    k_main = max(label_main) + 1
    k_relative = max(label_relative) + 1

    fig, ax = plt.subplots(nrows=k_main, ncols=k_relative,
                           sharex=True, sharey=True)
    cbar_ax = fig.add_axes([.91, .3, .025, .4])
    # cbar_ax = fig.add_axes([.1, .01, .8, .04])
    n_plots = 1
    for i in range(0, k_main):
        for j in range(0, k_relative):
            gene_set_tmp = gene_set.get((i, j))

            reduced_experiment = GeneExpression.copyFrom(experiment_main)
            reduced_experiment.expression =\
                reduced_experiment.expression.loc[gene_set_tmp[0]]

            # Create correlation map. Not using plotCorrelationMatrix
            # because of a lot of differences in plot
            correlation_set = []
            for gene in gene_set_tmp[1]:
                correlation_set.append(reduced_experiment.correlation(gene))

            correlation_set = pd.concat(correlation_set, axis=1)
            correlation_set.columns = gene_set_tmp[1]

            # Plot
            yticklabels, xticklabels = False, False

            if len(correlation_set.columns) < 20:
                yticklabels = True
            if len(correlation_set.index) < 20:
                xticklabels = True

            # plt.subplot(k_main, k_relative, n_plots)
            with sns.axes_style('dark'):
                sns.heatmap(correlation_set, ax=ax[i][j],
                            yticklabels=yticklabels, xticklabels=xticklabels,
                            cmap='seismic', vmin=-1, vmax=1,
                            cbar=True,
                            cbar_ax=cbar_ax)
                # cbar_kws={"orientation": "horizontal"}

            n_plots += 1

    plt.suptitle(experiment_main.name)
    fig.text(0.04, 0.5, "genes of " +
             experiment_main.name, rotation='vertical')
    fig.text(0.5, 0.04, "genes of " + experiment_relative.name, ha='center')

    # mappable = ax[0][k_relative - 1].get_children()[0]
    # plt.colorbar(mappable, ax, orientation='vertical')
    plt.show(block=False)


def pickleDumpLargeFile(obj, filepath):
    # This is a defensive way to write pickle.write, allowing for very large
    # files on all platforms
    max_bytes = 2**31 - 1
    bytes_out = pickle.dumps(obj)
    n_bytes = sys.getsizeof(bytes_out)
    with open(filepath, 'wb') as f_out:
        for idx in range(0, n_bytes, max_bytes):
            f_out.write(bytes_out[idx:idx + max_bytes])
