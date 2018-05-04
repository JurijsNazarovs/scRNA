library("gskb")

loadData <- function(data_name){
    e <- new.env()
    name <- data(list = data_name, envir = e)[1]
    e[[name]]
}

getPathways <- function(pathways, gene){
    res_pathways <- c()
    for (i in 1:length(pathways)){
        # is_found <- tolower(gene) %in% pathways[[i]]
        is_found <- tolower(gene) %in% lapply(pathways[[i]], tolower)
        if (is_found == TRUE){
            res_pathways <- c(res_pathways, names(pathways)[i])
        }
    }
    return(res_pathways)
}

## Input parameters
inp_file <- "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/mutant.A_rg.csv"
out_file <- "/Users/owner/Box Sync/UW/research/scRna/data_proceeded/mutant.A_rg_pathways.csv"
pathways_name <- "mm_pathway"

## Read data
pathways <- loadData(pathways_name) # possible pathways:
for (i in 1:length(pathways)){
     pathways[i] <- lapply(pathways[i], tolower)
}
#data <- read.csv(inp_file, header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
genes <- data[, 1]
data <- data[, -1]

## Create new table with pathways
all_pathways <- c()
table_w_pathways <- data.frame(matrix(nrow=0, ncol=dim(data)[2]))
for (i in 1:dim(data)[1]){
    pathways_tmp <- getPathways(pathways, genes[i])
    if (length(pathways_tmp) > 0){
        all_pathways <- c(all_pathways, pathways_tmp)
        for (j in 1:length(pathways_tmp)){
            table_w_pathways[dim(table_w_pathways)[1] + j, ] <- data[i, ]
        }
    }
}

names(table_w_pathways) <- all_pathways

## Summ wrt pathways (same pathways)
hui <- aggregate(table_w_pathnames, names(table_w_pathways), sum) 
