library(pheatmap)
library(ggplot2)
library(cluster)
library(factoextra)
library(Rtsne)
library(viridis)

## Input parameters
args = commandArgs(trailingOnly=TRUE)
experiment_path <- args[1]
experiment_name <- unlist(strsplit(experiment_path, "/"))
experiment_name <- experiment_name[length(experiment_name)]

out_path <- args[2]
dir.create(out_path)
out_file <- paste0(out_path, "/", experiment_name)

load(paste0(experiment_path))
#data_magic <- read.csv(paste0(experiment_name, ".csv"), stringsAsFactors = F,
#                        header = T, row.names = 1)
data_magic <- t(data_magic)


# Silhouette
pdf(paste0(out_file, '_silhouette.pdf'))
fviz_nbclust(data_magic, kmeans, method = "silhouette")
dev.off()

# k-means
n_clusters <- 10
kmm <- kmeans(data_magic, n_clusters)$cluster

# Compare k-means of raw data
# data_hui <- t(read.csv(paste0("wt.F_rg", ".csv"), stringsAsFactors = F,
#                        header = T, row.names = 1))
# kmm_raw <- kmeans(data_hui, n_clusters)$cluster
# print(mean(kmm == kmm_raw))


# t-sne
data_tsne <- Rtsne(data_magic, 2)$Y
data_plot <- data.frame(x = data_tsne[, 1], y = data_tsne[, 2], lab = factor(kmm))
pdf(paste0(out_file, "_tsne_kmeans_", n_clusters, ".pdf"))
ggplot(data_plot, aes(x=x, y=y, color = lab)) +
  geom_point()
dev.off()


if (1 == 0){
    ## Create breaks for pheatmap
    quantile_breaks <- function(xs, n = 10) {
        breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
        breaks[!duplicated(breaks)]
    }
    magic_breaks <- quantile_breaks(t(as.matrix(data_MAGIC)), n = 20)


    png(paste0(experiment_name, '_heatmap.png'))
    pheatmap(t(as.matrix(data_MAGIC)), 
             show_colnames = FALSE, 
             show_rownames = FALSE,
             scale = "row", #kmeans_k = 5,
             cluster_rows = FALSE,
             color             = inferno(length(magic_breaks) - 1),
             breaks            = magic_breaks,
             border_color      = NA,
             drop_levels       = TRUE,
             fontsize          = 14)
    dev.off()
}
