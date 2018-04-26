library(cellrangerRkit)
#pipestance_path <- "/Users/owner/Box Sync/UW/research/scRna/data/cellranger_10x/test/A"
pipestance_path <- "/Users/owner/Box Sync/UW/research/scRna/scripts/R"
# download_sample(
#   sample_name = "pbmc3k", sample_dir = pipestance_path,
#   host = "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/"
# )
gbm <- load_cellranger_matrix(pipestance_path)
analysis_results <- load_cellranger_analysis_results(pipestance_path)


tsne_proj <- analysis_results$tsne
visualize_umi_counts(gbm, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(3, 4), marker_size = 0.05)


# Filtering
use_genes <- get_nonzero_genes(gbm) # completely 0 genes in every cell
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes, ]) # - median/sum
gbm_log <- log_gene_bc_matrix(gbm_bcnorm, base = 10) # log(1 + x)
print(dim(gbm_log))

n_clu <- 2:10
km_res <- analysis_results$clustering # load pre-computed kmeans results
clu_res <- sapply(n_clu, function(x) km_res[[paste("kmeans", x, "clusters", sep = "_")]]$Cluster)
colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans", x, sep = "."))
visualize_clusters(clu_res, tsne_proj[c("TSNE.1", "TSNE.2")])

example_K <- 3
# number of clusters (use "Set3" for brewer.pal below if example_K > 8) 
example_col <- rev(brewer.pal(example_K,"Set2")) # customize plotting colors
cluster_result <- analysis_results$clustering[[paste("kmeans", example_K, "clusters", sep = "_")]]
visualize_clusters(cluster_result$Cluster, tsne_proj[c("TSNE.1", "TSNE.2")], colour = example_col)
# sort the cells by the cluster labels
cells_to_plot <- order_cell_by_clusters(gbm, cluster_result$Cluster)
# order the genes from most up-regulated to most down-regulated in each cluster
prioritized_genes <- prioritize_top_genes(gbm, cluster_result$Cluster, "sseq", min_mean = 0.5)

output_folder <- "/path_to_your_local_folder/pbmc_data_public/pbmc3k/gene_sets"
write_cluster_specific_genes(prioritized_genes, output_folder, n_genes = 10)

# create values and axis annotations for pheatmap
gbm_pheatmap(log_gene_bc_matrix(gbm), prioritized_genes, cells_to_plot, 
             n_genes = 3, colour = example_col, limits = c(-1, 2))
