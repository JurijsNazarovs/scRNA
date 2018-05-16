#!/bin/bash
inp_dir="/Users/owner/Box Sync/UW/research/scRna/selected_genes"
out_dir="/Users/owner/Box Sync/UW/research/scRna/plots_tmp/"
mkdir -p ${out_dir}

for file in "${inp_dir}"/*; do
  echo "> $file"
  all_exper="/Users/owner/Box Sync/UW/research/scRna/data_proceeded/removed_genes.pkl"
  reduced_exper="/Users/owner/Box Sync/UW/research/scRna/data_proceeded/norm_log_pca_tsne.pkl"
  while read gene; do
    echo "$gene"
    python3 bubble_plots.py "$gene" "$all_exper" "$reduced_exper" "${out_dir}/$(basename "$file")/"
  done < "$file"
done
