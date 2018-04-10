## our version of separability metric
ComputeSeparability <- function(input.data, cells.1, cells.2) {
  if (length(cells.1) < 3) {
    #warning("Fewer than 3 cells in the first group (cells.1)")
    return(NA)
  }
  if (length(cells.2) < 3) {
    #warning("Fewer than 3 cells in the second group (cells.2)")
    return(NA)
  }
  k <- min(c(length(cells.1), length(cells.2)))
  tnn <- nn2(data = input.data[c(cells.1, cells.2), ], k = k + 1)
  idx <- tnn$nn.idx[, -1]
  rownames(idx) <- c(cells.1, cells.2)
  correct_neighbors_c1 <- sapply(cells.1, function(x)
    {
      length(which(rownames(idx)[idx[x,]] %in% cells.1))
    }
  )
  correct_neighbors_c2 <- sapply(cells.2, function(x)
    {
      length(which(rownames(idx)[idx[x,]] %in% cells.2))
    }
  )
  return(mean(c(correct_neighbors_c1, correct_neighbors_c2)) / k)
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

