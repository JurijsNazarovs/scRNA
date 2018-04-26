library(Rmagic)
library(pheatmap)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(viridis)

## Input parameters
args = commandArgs(trailingOnly=TRUE)
experiment_path = args[1]
experiment_name = args[2]
print(experiment_path)
print(experiment_name)
exit()

## Data preparation
experiment <- read.csv(experiment_path, row.names = 1)
dat <- experiment

data_MAGIC <- run_magic(dat, rescale_percent=0.99)
rownames(data_MAGIC) <- rownames(dat)

## Create breaks for pheatmap
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
magic_breaks <- quantile_breaks(t(as.matrix(data_MAGIC)), n = 20)
dat_breaks <- quantile_breaks(t(dat), n = 20)


## pheatmap
png(paste0(experiment_name, '_raw.jpg'))
pheatmap(t(dat),  show_colnames = FALSE, show_rownames = FALSE, scale = "row", #kmeans_k = 5,
         color             = viridis(length(dat_breaks) - 1),
         breaks            = dat_breaks,
         border_color      = NA,
         drop_levels       = TRUE,
         fontsize          = 14)
dev.off()


png(paste0(experiment_name, '_magic.jpg'))
pheatmap(t(as.matrix(data_MAGIC)),  show_colnames = FALSE, show_rownames = FALSE, scale = "row", #kmeans_k = 5,
         color             = inferno(length(magic_breaks) - 1),
         breaks            = magic_breaks,
         border_color      = NA,
         drop_levels       = TRUE,
         fontsize          = 14)
dev.off()
