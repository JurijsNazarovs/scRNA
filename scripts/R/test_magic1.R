## ==================================
## Read in cellRanger processed data
## ==================================
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(cellrangerRkit)
library(SparseDC)

## Function to read in data: Neeed pipestance_path <-"/store01/scRna/results/10x/F/"
readRaw <- function(pipestance_path){
    gbm <- load_cellranger_matrix(pipestance_path)
    return(gbm)
}

## Function to filter genes based on min.read, min.cell and also low variability
filter_gene <- function(raw.gbm, min.read = 3, min.cell = 3, var.quantile = 0.05){
    cat("--------------------------------------", "\n")
    cat("Dimension before filter ", dim(raw.gbm), "\n")
    filter1 <- rowSums(exprs(raw.gbm)>min.read)>min.cell
    table(filter1)
    sub.gbm <- raw.gbm[filter1, ]
    cat("Dimension after gene min.read-min.cell filter ", dim(sub.gbm), "\n")
    
    as.matrix(exprs(sub.gbm)) %>% log1p %>% rowVars -> vars
    names(vars) <- rownames(exprs(sub.gbm))
    #vars <- vars[vars > quantile(vars, var.quantile)]
    sub.gbm.f <- sub.gbm[vars > quantile(vars, var.quantile), ] #discard bottom 5%

    cat("Dimension after var filter ", dim(sub.gbm.f), "\n")
    return(sub.gbm.f)
}

## Function to filter cells based on sequencing depth
filter_cell <- function(raw.gbm, min.cell.depth = 4000){
    cat("--------------------------------------", "\n")
    cat("Dimension before filter ", dim(raw.gbm), "\n")
    filter1 <- colSums(exprs(raw.gbm))>min.cell.depth
    table(filter1)
    sub.gbm <- raw.gbm[, filter1]
    cat("Dimension before filter ", dim(sub.gbm), "\n")
    return(sub.gbm)
}




# Raw counts
a.gbm <- readRaw("/store01/scRna/results/10x/A/")
b.gbm <- readRaw("/store01/scRna/results/10x/B/")
i.gbm <- readRaw("/store01/scRna/results/10x/I/")
f.gbm <- readRaw("/store01/scRna/results/10x/F/")

# Gene filtering: returns filtered genes by unfiltered cells
a.gbm.g.temp <- filter_gene(a.gbm, min.read = 3, min.cell = 3, var.quantile = 0.05)
b.gbm.g.temp <- filter_gene(b.gbm, min.read = 3, min.cell = 3, var.quantile = 0.05)
i.gbm.g.temp <- filter_gene(i.gbm, min.read = 3, min.cell = 3, var.quantile = 0.05)
f.gbm.g.temp <- filter_gene(f.gbm, min.read = 3, min.cell = 3, var.quantile = 0.05)

gSet <- unique(c(a.gbm.g.temp@featureData$id, b.gbm.g.temp@featureData$id, i.gbm.g.temp@featureData$id, f.gbm.g.temp@featureData$id))

# Cell filtering: returns unfiltered genes by filtered cells
a.gbm.c.temp <- filter_cell(a.gbm,  min.cell.depth = 4000)
b.gbm.c.temp <- filter_cell(b.gbm,  min.cell.depth = 4000)
i.gbm.c.temp <- filter_cell(i.gbm,  min.cell.depth = 4000)
f.gbm.c.temp <- filter_cell(f.gbm,  min.cell.depth = 4000)

# Get filtered genes from the cell filtering output
a.gbm.f<- a.gbm.c.temp[gSet, ]
b.gbm.f <- b.gbm.c.temp[gSet, ]
i.gbm.f <- i.gbm.c.temp[gSet, ]
f.gbm.f <- f.gbm.c.temp[gSet, ]


cat("Dimensions after all filtering ", dim(a.gbm.f), dim(b.gbm.f), dim(i.gbm.f), dim(f.gbm.f), "\n")

# save a version before normalization
save(a.gbm.f, b.gbm.f, i.gbm.f, f.gbm.f, file = "/store01/scRna/SK/a-b-i-f-filter-3-3-var0.05-cellD-4000-beforeNormalization.RData")


## ==================================
## Experiment with MAGIC

library(Rmagic)
load("/store01/scRna/SK/a-b-i-f-filter-3-3-var0.05-cellD-4000-beforeNormalization.RData")
dat <- as.matrix(exprs(i.gbm.f))
rownames(dat) <- i.gbm.f@featureData@data[, 2]

data_MAGIC <- run_magic(dat, rescale_percent=0.99)
rownames(data_MAGIC) <- i.gbm.f@featureData@data[, 2]

save(dat, data_MAGIC, file = "/store01/scRna/SK/magic_i.RData")

setwd("/Users/keles/Research_Collab/Bresnick-scRNA-seq/28March2018/")
load("magic_i.RData")
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



## ==================================

