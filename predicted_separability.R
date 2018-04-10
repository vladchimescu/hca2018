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

################################ Fxn Definitions ###############################
PredictSeparability <- function(num_cells) {
  features <- tibble(ncell = numeric(), 
    clusters = character(), 
    separability = double(),
    ncell_cluster1 = double(),
    ncell_cluster2 = double(),
    mean_nUMI = double(),
    tp = double(),
    kdist = double())  
  
  subsample <- SubsetData(object,cells.use = sample(object@cell.names,size = num_cells))
  
  #print(cluster_ids);
  #print(num_cells)
  #num_cells
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

ncell_use <- seq(1000, 100000, 1000)

predictions <- mclapply(ncell_use, PredictSeparability, mc.cores=4)
