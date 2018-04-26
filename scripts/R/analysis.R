library(Rmagic)
library(pheatmap)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)

data_path <- "/Users/owner/Box Sync/UW/research/scRna/data_analysis/"
experiment <- read.csv(paste0(data_path, "mutant.A.csv"), row.names = 1)

dat <- experiment

data_MAGIC <- run_magic(dat, rescale_percent=0.99)
rownames(data_MAGIC) <- rownames(dat)

b.dat <- data.frame(g1 = dat[ which(rownames(dat) == "Gata2"), ], g2 = dat[ which(rownames(dat) == "Pde6d"), ])
a.dat <- data.frame(g1 =  as.matrix(data_MAGIC)[ which(rownames(data_MAGIC) == "Gata2"), ], g2 = as.matrix(data_MAGIC)[ which(rownames(data_MAGIC) == "Pde6d"), ])

ggplot(b.dat) + geom_point(aes(g1, g2)) + scale_colour_gradient(low = 'purple', high='yellow')
ggplot(a.dat) + geom_point(aes(g1, g2)) + scale_colour_gradient(low = 'purple', high='yellow')

# kmeans on the columns (cells)

library(viridis)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

magic_breaks <- quantile_breaks(t(as.matrix(data_MAGIC)), n = 20)
dat_breaks <- quantile_breaks(t(dat), n = 20)


jpeg('/Users/keles/Research_Collab/Bresnick-scRNA-seq/28March2018/sample_i_magic_row_scaled.jpg')
pheatmap(t(dat),  show_colnames = FALSE, show_rownames = FALSE, scale = "row", #kmeans_k = 5,
         color             = viridis(length(dat_breaks) - 1),
         breaks            = dat_breaks,
         border_color      = NA,
         drop_levels       = TRUE,
         fontsize          = 14)
dev.off()


jpeg('/Users/keles/Research_Collab/Bresnick-scRNA-seq/28March2018/sample_i_magic_row_scaled.jpg')
pheatmap(t(as.matrix(data_MAGIC)),  show_colnames = FALSE, show_rownames = FALSE, scale = "row", #kmeans_k = 5,
         color             = inferno(length(magic_breaks) - 1),
         breaks            = magic_breaks,
         border_color      = NA,
         drop_levels       = TRUE,
         fontsize          = 14)
dev.off()