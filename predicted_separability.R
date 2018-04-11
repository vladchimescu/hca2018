#!/usr/bin/env Rscript

# This script takes a loom object, builds a model to predict separability, and
# outputs the separability predictions for 1k-50k cells (by 1k)

################################ Arguments #####################################
# args[1]    Path to loom file
# args[2]    Name of cluster 1 
# args[3]    Fraction of cluster 1 in the full dataset
# args[4]    Name of cluster 2 
# args[5]    Fraction of cluster 2 in the full dataset
# args[6]    Path for output tsv file
################################################################################

suppressMessages(library(Seurat))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(rhdf5))
suppressMessages(library(RANN))
suppressMessages(library(tidyr))
suppressMessages(library(randomForest))

load("final_model.RData") 


args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
cluster_name1 <- args[2]
freq1 <- args[3]
cluster_name2 <- args[4]
freq2 <- args[5]
outputfile <- args[6]

freq1 = as.double(freq1)
freq2 = as.double(freq2)


all_data <- NULL # will be filled up later

################################ Fxn Definitions ###############################
PredictSeparability <- function(num_cells) {
  features <- all_data
  features$ncell <- num_cells
  features$log2_ncell <- log2(num_cells)
  predict(final_model, features)
}

### Process Data
ProcessData <- function(object) {
  object <- FindVariableGenes(object, do.plot = FALSE, display.progress = FALSE)
  object@var.genes <- rownames(head(object@hvg.info, 1000))
  object <- ScaleData(object, genes.use = object@var.genes,
                      display.progress = FALSE, check.for.norm = FALSE)
  object <- RunPCA(object, pcs.compute = 30, do.print = FALSE)
  return(GetCellEmbeddings(object, reduction.type = "pca", dims.use = 1:30))
}

# calculaets true positive rate of cluster assignment
# clustering using k-means
get_tp <- function(cells_1, cells_2, data_for_knn) {
  if (length(cells_1) < 3) {
    #warning("Fewer than 3 cells in the first group (cells.1)")
    return(NA)
  }
  if (length(cells_2) < 3) {
    #warning("Fewer than 3 cells in the second group (cells.2)")
    return(NA)
  }
  # run kmeans with k=2 and 50 random starts
  km.results = kmeans(data_for_knn[c(cells_1, cells_2), ],  centers = 2, nstart = 50)
  # identify cells in each cluster (i=1,2)
  km.cl1 = names(km.results$cluster[km.results$cluster == 1])
  km.cl2 = names(km.results$cluster[km.results$cluster == 2])
  
  # compute true-positive rate
  tp = max(mean(c(sum(km.cl1 %in% cells_1)/ length(cells_1),
                  sum(km.cl2 %in% cells_2)/ length(cells_2))),
           mean(c(sum(km.cl1 %in% cells_2)/ length(cells_2),
                  sum(km.cl2 %in% cells_1)/ length(cells_1))))
  tp
}

# function for distance between kmeans centroids of two clusters
get_kmdist <- function(cells_1, cells_2, data_for_knn) {
  if (length(cells_1) < 3) {
    #warning("Fewer than 3 cells in the first group (cells.1)")
    return(NA)
  }
  if (length(cells_2) < 3) {
    #warning("Fewer than 3 cells in the second group (cells.2)")
    return(NA)
  }
  # run kmeans with k=2 and 50 random starts
  km.results = kmeans(data_for_knn[c(cells_1, cells_2), ],  centers = 2, nstart = 50)
  # write centroid distance
  vdiff = km.results$centers[1,] - km.results$centers[2,]
  vecnorm(vdiff)
  
}

# Euclidean vector norm
vecnorm <- function(x) sqrt(sum(x^2))
################################################################################

data <- h5read(filename, name="/")

# we can read this using loomR, but since that might not be installed, we'll read using HDF5
data_matrix <- t(data$matrix)
dimnames(data_matrix) <- list(data$row_attrs$gene_names, data$col_attrs$cell_names)
cluster_ids <- data$col_attrs$cluster; 
names(cluster_ids) <- colnames(data_matrix)
rm(data)
data_matrix <- Matrix(data_matrix,sparse=T)

object <- CreateSeuratObject(raw.data = data_matrix)
object <- NormalizeData(object,display.progress = F,scale.factor = 1e4,normalization.method = "LogNormalize")


### generate features
cluster_separabilities <- tibble(ncell = numeric(), 
                                 clusters = character(), 
                                 separability = double(),
                                 ncell_cluster1 = double(),
                                 ncell_cluster2 = double(),
                                 mean_nUMI = double(),
                                 tp = double(),
                                 kdist = double())

cluster_separabilities_local <- tibble(ncell = numeric(), 
                                         clusters = character(), 
                                         separability = double(),
                                         ncell_cluster1 = double(),
                                         ncell_cluster2 = double(),
                                         mean_nUMI = double(),
                                         tp = double(),
                                         kdist = double())
  
subsample = object

# mean library size
mean_nUMI = mean(subsample@meta.data[subsample@cell.names, 'nUMI'])
  
# Run ProcessData
data_for_knn = ProcessData(subsample)
  
subsample_cells <- subsample@cell.names

cells_1 <- intersect(subsample_cells,names(cluster_ids[which(cluster_ids==cluster_name1)]))
cells_2 <- intersect(subsample_cells,names(cluster_ids[which(cluster_ids==cluster_name2)]))
cluster_sep <- 0
tp = get_tp(data_for_knn = data_for_knn, cells_1 = cells_1, cells_2 = cells_2)

# k-means
kdist = get_kmdist(data_for_knn = data_for_knn, cells_1 = cells_1, cells_2 = cells_2)
cluster_comparison <- paste(cluster_name1,cluster_name2,sep=":")  
row_n <- nrow(cluster_separabilities_local)+1
cluster_separabilities_local[row_n,"ncell"] <- 0; # fill
cluster_separabilities_local[row_n,"separability"] <- 0  # do not compute for testing
cluster_separabilities_local[row_n,"clusters"] <- paste(cluster_name1,cluster_name2,sep=":")
# compute cluster balance
cluster_separabilities_local[row_n,"ncell_cluster1"] <- length(cells_1)
cluster_separabilities_local[row_n,"ncell_cluster2"] <- length(cells_2)
# compute mean coverage
cluster_separabilities_local[row_n,"mean_nUMI"] <- mean_nUMI
depth1s = subsample@meta.data[cells_1, 'nUMI']
depth2s = subsample@meta.data[cells_2, 'nUMI']
cluster_separabilities_local[row_n,"mean_depth_cluster1"] <- mean(depth1s)
cluster_separabilities_local[row_n,"median_depth_cluster1"] <- median(depth1s)
cluster_separabilities_local[row_n,"var_depth_cluster1"] <- var(depth1s)
cluster_separabilities_local[row_n,"mean_depth_cluster2"] <- mean(depth2s)
cluster_separabilities_local[row_n,"median_depth_cluster2"] <- median(depth2s)
cluster_separabilities_local[row_n,"var_depth_cluster2"] <- var(depth2s)

# k-means features
cluster_separabilities_local[row_n, "tp"] <- tp
cluster_separabilities_local[row_n, "kdist"] <- kdist
      
cluster_separabilities =  cluster_separabilities_local

## DE genes and fold changes
object <- SetIdent(object, names(cluster_ids),cluster_ids)
avg_cluster <- AverageExpression(object,return.seurat = T,show.progress = F)

NumDEGenes <- function(avg_cluster, cluster_comparison) {
  cluster_1 <- ExtractField(cluster_comparison,1,":")
  cluster_2 <- ExtractField(cluster_comparison,2,":")
  return(length(which(abs(avg_cluster@data[,cluster_1]-avg_cluster@data[,cluster_2])>log(2))))
}

foldChange <- function(avg_cluster, cluster_comparison) {
  cluster_1 <- ExtractField(cluster_comparison,1,":")
  cluster_2 <- ExtractField(cluster_comparison,2,":")
  fc = abs(avg_cluster@data[,cluster_1] - avg_cluster@data[,cluster_2])
  mean_fold_change = mean(fc)
  max_fold_change = max(fc)
  median_fold_change = median(fc)
  var_fold_change = var(fc)
  mad_fold_change = mad(fc)
  return(c(mean_fold_change=mean_fold_change,
           max_fold_change = median_fold_change,
           median_fold_change = median_fold_change,
           var_fold_change = var_fold_change,
           mad_fold_change = mad_fold_change))
}

all_data <- filter(cluster_separabilities, !is.na(separability))

num_de_genes <- sapply(1:nrow(all_data),function(x) NumDEGenes(avg_cluster, cluster_comparison = as.character(all_data$clusters[x])))
all_data$num_de_genes <- (num_de_genes)
all_data$frequency_cluster1 <- freq1
all_data$frequency_cluster2 <- freq2

fold_changes <- sapply(1:nrow(all_data),function(x) foldChange(avg_cluster, cluster_comparison = as.character(all_data$clusters[x])))
all_data$mean_fold_change <- fold_changes['mean_fold_change',]
all_data$max_fold_change <- fold_changes['max_fold_change',]
all_data$median_fold_change <- fold_changes['median_fold_change',]
all_data$var_fold_change <- fold_changes['var_fold_change',]
all_data$mad_fold_change <- fold_changes['mad_fold_change',]

### compute cluster balance
all_data$ncell_big_cluster <- pmax(all_data$ncell_cluster1, all_data$ncell_cluster2)
all_data$ncell_small_cluster <- pmin(all_data$ncell_cluster1, all_data$ncell_cluster2)
all_data$cluster_balance <- all_data$ncell_small_cluster / all_data$ncell_big_cluster

prep.dataset <- function(dataset){
    dataset$total.frequency <- dataset$frequency_cluster1 + dataset$frequency_cluster2
    dataset$frequency_big_cluster <- pmax(dataset$frequency_cluster1, dataset$frequency_cluster2)
    dataset$min.depth <- pmin(dataset$mean_depth_cluster1, dataset$mean_depth_cluster2)
    dataset$total.depth <- dataset$mean_depth_cluster1 + dataset$mean_depth_cluster2
    dataset$total.cells <- dataset$ncell_cluster1 + dataset$ncell_cluster2
    dataset$freq.ratios <- ifelse(dataset$mean_depth_cluster1 > dataset$mean_depth_cluster2, dataset$mean_depth_cluster1/dataset$mean_depth_cluster2, dataset$mean_depth_cluster2/dataset$mean_depth_cluster1)
    dataset$log2_mean_nUMI <- log2(dataset$mean_nUMI)
    dataset
}

all_data <- prep.dataset(all_data)

ncell_use <- seq(1000, 100000, 1000)

predictions <- lapply(ncell_use, PredictSeparability)
predictions <- do.call(rbind, predictions)

output <- data.frame(ncell_use, predictions)
colnames(output) <- c("Number of Cells", "Predicted Separability")
write.table(x = output, file = outputfile, sep = "\t", quote = FALSE, row.names = FALSE)
