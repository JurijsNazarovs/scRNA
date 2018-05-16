library("ggplot2")
setwd ("/Users/owner/tmp/rg_pathways")

all_files <- list.files(path = ".", pattern ="*.csv")
data_plot <- data.frame(val = double(), name = factor())
for (file in all_files){
  data_tmp <- read.csv(file,  header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  data_tmp <- data_tmp[, -1]
  val <- colSums(data_tmp == 0)/dim(data_tmp)[1]
  data_plot <- rbind(data_plot, 
                     data.frame(val = val, 
                                name = as.factor(rep(strsplit(file, "_")[[1]][1], length(val)))))
}

pdf("rg_pathways.pdf")
ggplot(data_plot, aes(x = val, fill = name)) + 
  geom_histogram(alpha = 0.5, bins = 50, position = "identity") +
  labs( x="Proprotion of genes with 0 expresison", y="Number of cells")
dev.off()      