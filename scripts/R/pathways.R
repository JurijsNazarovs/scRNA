## Script to crate pathways accroding to a gskb library
## Based on pathways_name pathways
## Input: filtered experiment
## Output: experiment with corresponding to genes pathways in .csv format
##         with summ of values for same pathways
library("gskb")

loadData <- function(data_name){
    ## Function to load from pathways_name
    ## otherwise I cant use pathways_name as a variable
    ## for pathways
    e <- new.env()
    name <- data(list = data_name, envir = e)[1]
    e[[name]]
}

getPathways <- function(pathways, gene){
    res_pathways <- c()
    for (i in 1:length(pathways)){
        is_found <- tolower(gene) %in% pathways[[i]]
        #is_found <- tolower(gene) %in% lapply(pathways[[i]], tolower)
        if (is_found == TRUE){
            res_pathways <- c(res_pathways, names(pathways)[i])
        }
    }
    return(res_pathways)
}

## Input parameters
pathways_name <- "mm_pathway"

## experiment_path <- "/u/n/a/nazarovs/private/scRNA/data/mutant.A_rg.csv"
## out_file <- "/u/n/a/nazarovs/private/scRNA/data/mutant.A_rg_pathways.csv"

args = commandArgs(trailingOnly=TRUE)
experiment_path <- args[1]
experiment_name <- unlist(strsplit(experiment_path, "/"))
experiment_name <- experiment_name[length(experiment_name)]
out_path <- args[2]
out_file <- paste0(out_path, "/", experiment_name)

print(paste0("> Started for ", experiment_path)

## Read data
pathways <- loadData(pathways_name) # possible pathways:
for (i in 1:length(pathways)){
     pathways[i] <- lapply(pathways[i], tolower)
}

data <- read.csv(experiment_path, header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
genes <- data[, 1]
data <- as.matrix(data[, -1])

## Create a new table with pathways
all_pathways <- c()
table_w_pathways <- data.frame(matrix(nrow=0, ncol=dim(data)[2]))
for (i in 1:dim(data)[1]){
    pathways_tmp <- getPathways(pathways, genes[i])
    if (length(pathways_tmp) > 0){
        all_pathways <- c(all_pathways, pathways_tmp)
        for (j in 1:length(pathways_tmp)){
            table_w_pathways[dim(table_w_pathways)[1] + 1, ] <- data[i, ]
        }
    }
}

table_w_pathways$pathways <- all_pathways

## Sum wrt pathways (same pathways)
table_w_pathways <- aggregate(. ~ pathways, data = table_w_pathways, FUN = sum)
write.table(table_w_pathways, file = out_file, quote = FALSE, row.names=FALSE, sep=",")
write.table(table_w_pathways$pathways, file = paste0(out_file, "_listOfpathways"), quote=FALSE,
            row.names=FALSE, col.names=FALSE, sep=",")
save(table_w_pathways, file = paste0(out_file, ".RData"))
