import pandas as pd
import numpy as np
import scipy as sc
import scipy.io
import os
import csv
import glob
import copy
import pickle

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
                 condition=None,
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
            condition = glob.glob(sample_dir_path + "/type.*")[0]
            condition = condition.split(
                sample_dir_path + "/type.")[1]

            # Create a GeneExpression object
            gene_expression = GeneExpression(expression_matrix=expression_mat,
                                             gene_names=gene_names,
                                             cell_names=cell_names,
                                             condition=condition,
                                             sample=sample)
            # Filter
            geneExpression.collapse()
            gene_expressions.append(gene_expression)

        return(geneExpressions)

    # Instance Methods
    def remove(self, list_remove, what_remove="genes",):
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
            self.expression.drop(final_list_remove, axis=axis, inplace=True)
            print(self.condition + "." + self.sample + ": removed " +
                  len(final_list_remove) + " " + what_remove)

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
        else:
            pring("sign ", sign, " is not supoprted")
            return(None)

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
        # import pdb
        # pdb.set_trace()
        all_genes = list(self.expression.index)

        if not gene in self.expression.index:
            correlation = [np.nan for i in range(0, len(all_genes))]
            correlation = pd.DataFrame(correlation, columns=[gene])
            correlation.index = all_genes
            return(correlation)

        # del rest_genes[rest_genes.index(gene)]
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
                         plot_type=["violin"]):

        if isinstance(plot_type, str):
            plot_type = [plot_type]

        if isinstance(list_gene_expression, GeneExpression):  # does not work
            list_gene_expression = [list_gene_expression]

        uniq_id = [i.condition + "." + i.sample for i in list_gene_expression]
        gene_expression = [i.expression.loc[gene]
                           for i in list_gene_expression]

        for plt_type in plot_type:
            if plt_type == "violin":
                ax = sns.violinplot(data=gene_expression)
                ax.set(xticklabels=uniq_id,
                       ylabel="count")

            if plt_type == "boxplot":
                ax = sns.boxplot(data=gene_expression)
                ax.set(xticklabels=uniq_id,
                       ylabel="count")

            if plt_type == "hist":
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

            if plt_type == "kdeplot":
                for i in range(0, len(uniq_id)):
                    sns.kdeplot(gene_expression[i].values,
                                label=uniq_id[i], shade=True)

                plt.legend(loc='upper right')

            if plt_type == "distplot":
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

        if isinstance(plot_type, str):
            plot_type = [plot_type]

        if isinstance(list_gene_expression, GeneExpression):  # does not work
            list_gene_expression = [list_gene_expression]

        uniq_id = [i.condition + "." + i.sample for i in list_gene_expression]
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
            if plt_type == "scatter":
                for col in correlation_set.columns:
                    y = correlation_set.loc[:, col].values
                    x = np.arange(0, len(y))
                    sns.regplot(x, y, fit_reg=False, dropna=True,
                                label=col)
                    plt.legend()
        plt.title("Correlation of " + gene)
        plt.xticks()
        plt.show(block=False)


# Different methods not related to a class
def intersect(list_of_lists):
    if not isinstance(list_of_lists, list):
        print("Provided argumnet is not a list")
        return(None)
    new_list = set(list_of_lists[0])
    for i in range(1, len(list_of_lists)):
        new_list = list(new_list & set(list_of_lists[i]))
    return(new_list)


def removeElementFromList(list_, element):
    if element in list_:
        list_ = [i for i in list_ if i != element]
    return list_


def plotGenesPatterns(original_expression, reduced_expression, genes):
    cmaps = ['Blues', 'Reds', 'Greens', 'Purples']
    if reduced_expression.expression.shape[1] != 2 or \
            reduced_expression.expression.shape[0] == 0:
        print("reduced_expression should have shape (x !=0, 2)")
        return None
    if isinstance(genes, str):
        genes = [genes]
    if len(genes) > len(cmaps):
        print("Just " + len(cmaps) + " genes can be plotted")
        genes = genes[0:len(cmaps)]

    # Main plot
    plt.scatter(reduced_expression.expression.loc[:, 0],
                reduced_expression.expression.loc[:, 1],
                label="",
                facecolors='none',
                edgecolor='grey')
    plt.title(original_expression.condition + "." + original_expression.sample)

    size_step = 100
    size = len(genes) * size_step + 10
    legend_color = []
    for i in range(0, len(genes)):
        cells = original_expression.expression.loc[genes[i]] > 0
        if sum(cells) == 0:
            print("No gene " + genes[i] + " is presented")
            continue

        colors = np.array(
            original_expression.expression.loc[genes[i], cells].tolist())

        cmap = plt.get_cmap(cmaps[i])
        colors_of_cmap = cmap(np.linspace(0.2, 1, cmap.N))  # del white part

        # Create a new colormap from those colors
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            None, colors_of_cmap)

        norm = mpl.colors.BoundaryNorm(
            np.arange(0.5, colors.max() + 1, 1), cmap.N)

        plt.scatter(reduced_expression.expression.loc[cells, 0],
                    reduced_expression.expression.loc[cells, 1],
                    c=colors,
                    norm=norm,
                    label=(genes[i] + ": " + str(colors.max())),
                    edgecolor="black",
                    cmap=cmap,
                    alpha=0.9,
                    s=size)
        size -= size_step
        legend_color.append(cmap(norm(cmap.N)))
        if len(genes) == 1:
            if (colors.max() - 1 < 10):
                plt.colorbar(ticks=np.arange(1, colors.max() + 1, 1),
                             orientation="horizontal")
            else:
                plt.colorbar(ticks=np.arange(
                    1, colors.max() + 1, colors.max() // 10),
                    orientation="horizontal")

    # Set legend color
    plt.legend(loc='upper right')
    ax = plt.gca()
    leg = ax.get_legend()
    for i in range(0, len(legend_color)):
        leg.legendHandles[i].set_color(legend_color[i])

    plt.show(block=False)


def plot2DLabels(expression, labels):
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

    # plt.colorbar(ticks=np.arange(labels.min(), labels.max() + 1, 1))
    plt.title(expression.condition + "." +
              expression.sample + ":  k = " + str(labels.max() + 1))
    plt.show(block=False)


def cluster(self, k=2, method="kmeans"):
    if method == "kmeans":
        kmeans_model = KMeans(
            n_clusters=k, random_state=1).fit(self.expression)
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
    # chose genes in cluster expression of which > threshold
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

    with sns.axes_style('dark'):
        sns.heatmap(correlation_set, yticklabels=False,
                    cmap='seismic', vmin=-1, vmax=1)
    plt.title(self.condition + "." + self.sample)
    plt.xticks()
    # plt.show(block=False)


# hui = GeneExpression.loadFrom()
# # GeneExpression.plotDistribution([hui[0]], "Gata2", "hist")
# # GeneExpression.plotCorrelation(hui, "Gata2")


# rg2 = GeneExpression.copyFrom(hui[0])
# genes = rg2.getToRemove(percentage_of_0=0.20, what_remove="genes")
# genes = removeElementFromList(genes, "Gata2")
# rg2.remove(genes, what_remove="genes")
# cells = rg2.getToRemove(percentage_of_0=0.20, what_remove="cells")
# rg2.remove(cells, what_remove="cells")
# for express in hui:
#     express.remove(genes, what_remove="genes")

# GeneExpression.plotCorrelation(hui, "Gata2", "heatmap")
# GeneExpression.plotDistribution(hui, "Gata2", "violin")

# su2 = GeneExpression.copyFrom(rg2)
# rg2.reduceDimension(p=2, method="tsne")

# #plotGenesPatterns(su, rg, ["Rpl13", "Gata2", "Mrpl30"])
# #plotGenesPatterns(su, rg, ["Gata2"])
# plotGenesPatterns(su2, rg2, ["Rpl13", "Gata2"])

# #PIK = "pickle.dat"
# #pickle.dump(rg, open(PIK, "wb"))
# #pizda = pickle.load(open(PIK, "rb"))


# labels2 = cluster(rg2, p=3)
# plot2DLabels(rg2, labels2)

# cm = confusionMatrix(rg, labels, rg2, labels2)


# cm = confusionMatrixUsingGenes(su2, labels2, su2, labels2, threshold=3)
# print(cm)


# plotCorrelationMatrix(su2, su2.expression.index[1:10])
# su2.expression.index