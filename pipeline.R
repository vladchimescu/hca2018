library(Seurat)
library(rhdf5)
library(RANN)
library(tidyr)
library(dplyr)
library(argparser)
source('scratch/Group7/functions.R')

## to run the script in terminal:
# Rscript pipeline.R -loom '/home/jovyan/data/tasks/how_many_cells/pbmc.loom' \
# -ncell '250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500'
# -o /home/jovyan/scratch/Group7/results/user

args <- arg_parser("program");
args <- add_argument(args, '-loom',
                     help='loom file',
                     #default='/home/jovyan/data/tasks/how_many_cells/human_pancreas.loom')
                     #default='/home/jovyan/data/tasks/how_many_cells/bipolar.loom')
                     default='/home/jovyan/data/tasks/how_many_cells/pbmc.loom')
# args <- add_argument(args, '-separability',
#                      help='type of separability function: rahul / tp',
#                      default='rahul')
args <- add_argument(args, '-ncell',
                     help='number of cells to downsample -- comma separated',
                     #default='250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500')
                     #default="250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000")
                     default="250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,25000,26000,27000,28000,29000,30000")
args <- add_argument(args, '-o',
                     help='output prefix.',
                     # default='/home/jovyan/scratch/Group7/ashis/human_pancreas/human_pancreas')
                     #default='/home/jovyan/scratch/Group7/ashis/bipolar/bipolar')
                     #default='/home/jovyan/scratch/Group7/ashis/pbmc/pbmc')
                     default='/home/jovyan/scratch/Group7/ashis/output')


argv = parse_args(args)
loom_file <- argv$loom
#separability_func <- argv$separability
ncell_input = argv$ncell
out_prefix <- argv$o

### parse ncell input
parts <- strsplit(ncell_input, ',')[[1]]
valid_parts <- as.vector(sapply(parts, function(p) nchar(p)>0))
ncell_use <- as.integer(parts[valid_parts])

### prepare plot file
plt_fn = paste0(out_prefix, '_plots.pdf')
pdf(plt_fn)

#loom_file <- "/home/jovyan/data/tasks/how_many_cells/human_pancreas.loom"

# we can read this using loomR, but since that might not be installed, we'll read using HDF5
data <- h5read(loom_file,name = "/")
data_matrix <- t(data$matrix)
dimnames(data_matrix) <- list(data$row_attrs$gene_names, data$col_attrs$cell_names)
cluster_ids <- data$col_attrs$cluster; names(cluster_ids) <- colnames(data_matrix)
rm(data)
data_matrix <- Matrix(data_matrix,sparse=T)

#view a few lines of the UMI matrix
head(data_matrix[20:25,1:3])
print(dim(data_matrix))
data.frame(table(cluster_ids))

### make seurat and process data
object <- CreateSeuratObject(raw.data = data_matrix)
object <- NormalizeData(object,display.progress = F,scale.factor = 1e4,normalization.method = "LogNormalize")

ProcessData <- function(object) {
  object <- FindVariableGenes(object, do.plot = FALSE, display.progress = FALSE)
  object@var.genes <- rownames(head(object@hvg.info, 1000))
  object <- ScaleData(object, genes.use = object@var.genes,
                      display.progress = FALSE, check.for.norm = FALSE)
  object <- RunPCA(object, pcs.compute = 30, do.print = FALSE)
  return(GetCellEmbeddings(object, reduction.type = "pca", dims.use = 1:30))
}

ComputeSeparability <- function(input.data, cells.1, cells.2, k = 20) {
  if (length(cells.1) < 3) {
    #warning("Fewer than 3 cells in the first group (cells.1)")
    return(NA)
  }
  if (length(cells.2) < 3) {
    #warning("Fewer than 3 cells in the second group (cells.2)")
    return(NA)
  }
  k <- min(c(k, length(cells.1), length(cells.2)))
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

results=list()
cluster_separabilities <- tibble(ncell = numeric(), 
                                 clusters = character(), 
                                 separability = double(),
                                 ncell_cluster1 = double(),
                                 ncell_cluster2 = double(),
                                 mean_nUMI = double(),
                                 tp = double(),
                                 kdist = double())
subsample_cells_store <- list()

for(n_cell in 1:length(ncell_use)) {
  #Downsample Seurat object
  downsample_cells <- ncell_use[n_cell]
  subsample <- SubsetData(object,cells.use = sample(object@cell.names,size = downsample_cells))
  
  # mean library size
  mean_nUMI = mean(subsample@meta.data[subsample@cell.names, 'nUMI'])
  
  #Run ProcessData
  data_for_knn = ProcessData(subsample)
  
  #Now, compute all pairs of separabilities for all clusters
  clusters_compute <- unique(cluster_ids)
  subsample_cells <- subsample@cell.names
  subsample_cells_store[[as.character(n_cell)]] <- subsample_cells
  
  for(i in 1:length(clusters_compute)) {
    for(j in 1:i) {
      if (j<i) {
        cells_1 <- intersect(subsample_cells,names(cluster_ids[which(cluster_ids==clusters_compute[i])]))
        cells_2 <- intersect(subsample_cells,names(cluster_ids[which(cluster_ids==clusters_compute[j])]))
        cluster_sep <- ComputeSeparability(data_for_knn,cells_1, cells_2, k = 20)
        tp = get_tp(data_for_knn = data_for_knn, cells_1 = cells_1, cells_2 = cells_2)
        # k-means
        kdist = get_kmdist(data_for_knn = data_for_knn, cells_1 = cells_1, cells_2 = cells_2)
        cluster_comparison <- paste(clusters_compute[i],clusters_compute[j],sep=":")  
        row_n <- nrow(cluster_separabilities)+1
        cluster_separabilities[row_n,"ncell"] <- downsample_cells; 
        cluster_separabilities[row_n,"separability"] <- cluster_sep
        cluster_separabilities[row_n,"clusters"] <- paste(clusters_compute[i],clusters_compute[j],sep=":")
        # compute cluster balance
        cluster_separabilities[row_n,"ncell_cluster1"] <- length(cells_1)
        cluster_separabilities[row_n,"ncell_cluster2"] <- length(cells_2)
        # compute mean coverage
        cluster_separabilities[row_n,"mean_nUMI"] <- mean_nUMI
        # k-means features
        cluster_separabilities[row_n, "tp"] <- tp
        cluster_separabilities[row_n, "kdist"] <- kdist
      }
    }
  }
  if (downsample_cells%%1000 == 0) print(downsample_cells)
}

print(head(cluster_separabilities))
print(tail(cluster_separabilities))

### plot 
all_data <- filter(cluster_separabilities, !is.na(separability))
test <- filter(all_data, clusters %in% c("gamma:delta","gamma:acinar","quiescent_stellate:activated_stellate","ductal:delta","mast:macrophage")) 
ggplot(test,aes(x=ncell, y=separability, col=clusters))+geom_point() + geom_line()

object <- SetIdent(object, names(cluster_ids),cluster_ids)
avg_cluster <- AverageExpression(object,return.seurat = T,show.progress = F)
# par(mfrow=c(2,2))
# CellPlot(avg_cluster,"gamma","delta")
# CellPlot(avg_cluster,"gamma","acinar")
# CellPlot(avg_cluster,"quiescent_stellate","activated_stellate")
# CellPlot(avg_cluster,"ductal","delta")
# par(mfrow=c(1,1))

### save data
data_fn = paste0(out_prefix, '_data1.RData')
save.image(data_fn)

### num DE genes and log fold change
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

### save data
data_fn = paste0(out_prefix, '_data2.RData')
save.image(data_fn)


### plot separability vs cluster balance : all clusters together
ggplot(all_data,aes(x=cluster_balance, y=separability, col=ncell))+geom_point()

### prediction evaluator
evaluate_predictions <- function(predicted, actual){
  predicted[predicted>1] = 1
  pearson = cor.test(actual, predicted)
  return(list(r=as.numeric(pearson$estimate), p=pearson$p.value))
}

plot_lm_predictions <- function(predicted, actual, r, title){
  plot(predicted, actual, xlab="Predicted separation", ylab="Actual separation",pch=16,cex=0.2)
  abline(a = 0, b=1)
  title(title)
  legend( 'bottomright', sprintf('r: %.3f', r), bg = rgb(1,1,1,0.5))
}

### given model
provided_model <- lm(separability ~ num_de_genes + frequency_cluster1 + frequency_cluster2 + ncell ,data = all_data)
coef(provided_model)
drop1(provided_model)
provided_model_performance = evaluate_predictions(provided_model$fitted.values, all_data$separability)
plot_lm_predictions(provided_model$fitted.values, 
                    all_data$separability, 
                    provided_model_performance$r, 
                    "separability ~ num_de_genes + frequency_cluster1 \n+ frequency_cluster2 + ncell")

### model1
model1 <- lm(separability ~ num_de_genes + frequency_cluster1 + frequency_cluster2 + 
               ncell + ncell_big_cluster + ncell_small_cluster + 
               mean_nUMI + mean_fold_change + max_fold_change + 
               median_fold_change + var_fold_change + mad_fold_change + 
               cluster_balance, 
             data = all_data)
coef(model1)
drop1(model1)
model1_performance = evaluate_predictions(model1$fitted.values, all_data$separability)
plot_lm_predictions(model1$fitted.values, 
                    all_data$separability, 
                    model1_performance$r, 
                    "separability ~ num_de_genes + freq_c1 + freq_c2 + ncell + 
                    ncell_big + ncell_small + mean_nUMI + mean_fold + max_fold + 
                    median_fold + var_fold + mad_fold + cluster_balance")


### model2
model2 <- lm(separability ~ num_de_genes + frequency_cluster1 + frequency_cluster2 + 
               ncell + ncell_big_cluster + ncell_small_cluster + 
               mean_nUMI + mean_fold_change + max_fold_change + 
               median_fold_change + var_fold_change + mad_fold_change + 
               cluster_balance + ncell * cluster_balance, 
             data = all_data)
coef(model2)
drop1(model2)
model2_performance = evaluate_predictions(model2$fitted.values, all_data$separability)
plot_lm_predictions(model2$fitted.values, 
                    all_data$separability, 
                    model2_performance$r, 
                    "separability ~ num_de_genes + freq_c1 + freq_c2 + ncell + 
                    ncell_big + ncell_small + mean_nUMI + mean_fold + max_fold + 
                    median_fold + var_fold + mad_fold + cluster_balance + ncell * cluster_balance")
summary(model2)

### model3
model3 <- lm(separability ~ num_de_genes + frequency_cluster1 + frequency_cluster2 + 
               ncell + ncell_big_cluster + ncell_small_cluster + 
               median_fold_change + var_fold_change +
               cluster_balance + ncell * cluster_balance, 
             data = all_data)
coef(model3)
drop1(model3)
model3_performance = evaluate_predictions(model3$fitted.values, all_data$separability)
plot_lm_predictions(model3$fitted.values, 
                    all_data$separability, 
                    model3_performance$r, 
                    "separability ~ num_de_genes + freq_c1 + freq_c2 + ncell + 
                    ncell_big + ncell_small + median_fold + var_fold + 
                    cluster_balance + ncell * cluster_balance")

### model4
model4 <- lm(separability ~ num_de_genes + frequency_cluster1 + frequency_cluster2 + 
               log2(ncell) + log2(ncell_big_cluster) + log2(ncell_small_cluster) + 
               median_fold_change + var_fold_change +
               cluster_balance + log2(ncell) * cluster_balance, 
             data = all_data)
coef(model4)
drop1(model4)
model4_performance = evaluate_predictions(model4$fitted.values, all_data$separability)
plot_lm_predictions(model4$fitted.values, 
                    all_data$separability, 
                    model4_performance$r, 
                    "separability ~ num_de_genes + freq_c1 + freq_c2 + 
                    log2(ncell) + log2(ncell_big) + log2(ncell_small) + 
                    median_fold + var_fold + cluster_balance + log2(ncell) * cluster_balance")

### model5
model5 <- lm(separability ~ num_de_genes + log2(ncell) + log2(ncell_big_cluster) + log2(ncell_small_cluster) + 
               median_fold_change + var_fold_change +
               cluster_balance + log2(ncell) * cluster_balance, 
             data = all_data)
coef(model5)
drop1(model5)
model5_performance = evaluate_predictions(model5$fitted.values, all_data$separability)
plot_lm_predictions(model5$fitted.values, 
                    all_data$separability, 
                    model5_performance$r, 
                    "separability ~ num_de_genes + log2(ncell) + log2(ncell_big) +
                    log2(ncell_small) + median_fold + var_fold +
                    cluster_balance + log2(ncell) * cluster_balance")

### model6
model6 <- lm(separability ~ log2(ncell_big_cluster) + cluster_balance + log2(ncell) : cluster_balance,
             data = all_data)
coef(model6)
drop1(model6)
model6_performance = evaluate_predictions(model6$fitted.values, all_data$separability)
plot_lm_predictions(model6$fitted.values, 
                    all_data$separability, 
                    model6_performance$r, 
                    "separability ~ log2(ncell_big_cluster) + cluster_balance
                    + log2(ncell) : cluster_balance")
summary(model6)


### save data
data_fn = paste0(out_prefix, '_data3.RData')
save.image(data_fn)

dev.off()
