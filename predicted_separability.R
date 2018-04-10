#!/usr/bin/env Rscript

# This script takes a loom object, builds a model to predict separability, and
# outputs the separability predictions for 1k-100k cells (by 1k)

################################ Arguments #####################################
# args[1]    Path to loom file
# args[2]    Fraction of cluster 1 in the full dataset
# args[3]    Fraction of cluster 2 in the full dataset
# args[4]    Path for output tsv file
################################################################################

suppressMessages(library(Seurat))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
library(rhdf5)
library(RANN)
library(tidyr)

all_data <- NULL # will be filled up later
################################ Fxn Definitions ###############################
PredictSeparability <- function(num_cells) {
  # features <- tibble(ncell = numeric(), 
  #   clusters = character(), 
  #   separability = double(),
  #   ncell_cluster1 = double(),
  #   ncell_cluster2 = double(),
  #   mean_nUMI = double(),
  #   tp = double(),
  #   kdist = double())  
  # 
  # subsample <- SubsetData(object,cells.use = sample(object@cell.names,size = num_cells))
  
  #print(cluster_ids);
  #print(num_cells)
  #num_cells
  
  features = all_data
  features$ncell = num_cells
  
}
################################################################################


args <- commandArgs(trailingOnly = TRUE)
filename = args[1]
filename = "/home/jovyan/data/tasks/how_many_cells/human_pancreas_activated_quiescent_stellate.loom"

data <- h5read(filename, name="/")

# we can read this using loomR, but since that might not be installed, we'll read using HDF5
data_matrix <- t(data$matrix)
dimnames(data_matrix) <- list(data$row_attrs$gene_names, data$col_attrs$cell_names)
cluster_ids <- data$col_attrs$cluster; names(cluster_ids) <- colnames(data_matrix)
rm(data)
data_matrix <- Matrix(data_matrix,sparse=T)

#view a few lines of the UMI matrix
head(data_matrix[20:25,1:3])


object <- CreateSeuratObject(raw.data = data_matrix)
object <- NormalizeData(object,display.progress = F,scale.factor = 1e4,normalization.method = "LogNormalize")

### process data
ProcessData <- function(object) {
  object <- FindVariableGenes(object, do.plot = FALSE, display.progress = FALSE)
  object@var.genes <- rownames(head(object@hvg.info, 1000))
  object <- ScaleData(object, genes.use = object@var.genes,
                      display.progress = FALSE, check.for.norm = FALSE)
  object <- RunPCA(object, pcs.compute = 30, do.print = FALSE)
  return(GetCellEmbeddings(object, reduction.type = "pca", dims.use = 1:30))
}

### functions
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


### generate features
cluster_separabilities <- tibble(ncell = numeric(), 
                                 clusters = character(), 
                                 separability = double(),
                                 ncell_cluster1 = double(),
                                 ncell_cluster2 = double(),
                                 mean_nUMI = double(),
                                 tp = double(),
                                 kdist = double())
subsample_cells_store <- list()

  cluster_separabilities_local <- tibble(ncell = numeric(), 
                                         clusters = character(), 
                                         separability = double(),
                                         ncell_cluster1 = double(),
                                         ncell_cluster2 = double(),
                                         mean_nUMI = double(),
                                         tp = double(),
                                         kdist = double())
  
  #Downsample Seurat object
  # downsample_cells <- ncell_use[n_cell]
  #subsample <- SubsetData(object,cells.use = sample(object@cell.names,size = downsample_cells))
  subsample = object
  # mean library size
  mean_nUMI = mean(subsample@meta.data[subsample@cell.names, 'nUMI'])
  
  #Run ProcessData
  data_for_knn = ProcessData(subsample)
  
  #Now, compute all pairs of separabilities for all clusters
  clusters_compute <- unique(cluster_ids)
  subsample_cells <- subsample@cell.names
  # subsample_cells_store[[as.character(n_cell)]] <- subsample_cells
  
  for(i in 1:length(clusters_compute)) {
    for(j in 1:i) {
      if (j<i) {
        cells_1 <- intersect(subsample_cells,names(cluster_ids[which(cluster_ids==clusters_compute[i])]))
        cells_2 <- intersect(subsample_cells,names(cluster_ids[which(cluster_ids==clusters_compute[j])]))
        cluster_sep <- 0
        tp = get_tp(data_for_knn = data_for_knn, cells_1 = cells_1, cells_2 = cells_2)
        # k-means
        kdist = get_kmdist(data_for_knn = data_for_knn, cells_1 = cells_1, cells_2 = cells_2)
        cluster_comparison <- paste(clusters_compute[i],clusters_compute[j],sep=":")  
        row_n <- nrow(cluster_separabilities_local)+1
        cluster_separabilities_local[row_n,"ncell"] <- 0; # fill
        cluster_separabilities_local[row_n,"separability"] <- 0  # do not compute for testing
        cluster_separabilities_local[row_n,"clusters"] <- paste(clusters_compute[i],clusters_compute[j],sep=":")
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
      }
    }
  }
  
#cluster_separabilities <- mclapply(ncell_use, fn1, mc.cores=n_cores)
#cluster_separabilities <- do.call(rbind, cluster_separabilities)

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
all_data$frequency_cluster1 <- sapply(1:nrow(all_data),function(x) length(which(cluster_ids==ExtractField(as.character(all_data$clusters[x]),1,":"))) / length(cluster_ids))
all_data$frequency_cluster2 <- sapply(1:nrow(all_data),function(x) length(which(cluster_ids==ExtractField(as.character(all_data$clusters[x]),2,":"))) / length(cluster_ids))

fold_changes <- sapply(1:nrow(all_data),function(x) foldChange(avg_cluster, cluster_comparison = as.character(all_data$clusters[x])))
all_data$mean_fold_change <- fold_changes['mean_fold_change',]
all_data$max_fold_change <- fold_changes['max_fold_change',]
all_data$median_fold_change <- fold_changes['median_fold_change',]
all_data$var_fold_change <- fold_changes['var_fold_change',]
all_data$mad_fold_change <- fold_changes['mad_fold_change',]

print(head(all_data))

### compute cluster balance
all_data$ncell_big_cluster <- pmax(all_data$ncell_cluster1, all_data$ncell_cluster2)
all_data$ncell_small_cluster <- pmin(all_data$ncell_cluster1, all_data$ncell_cluster2)
all_data$cluster_balance <- all_data$ncell_small_cluster / all_data$ncell_big_cluster
all_data[1:5,]


ncell_use <- seq(1000, 100000, 1000)

predictions <- mclapply(ncell_use, PredictSeparability, mc.cores=4)
