### this script creates CAT-plots to visualize hub agreements between tissues

library(argparser)
library(AUC)

args <- arg_parser("program");
args <- add_argument(args, '-hub',
                     help='hub files, semi-colon (;) separated.',
                     #default='/Users/ashissaha/github/gtex_twn/results/quic_scale_free_combined/filtered_conflicts/hubs/WholeBlood_hub_exp_iso.txt;/Users/ashissaha/github/gtex_twn/results/quic_scale_free_combined/filtered_conflicts/hubs/Muscle-Skeletal_hub_exp_iso.txt;/Users/ashissaha/github/gtex_twn/results/quic_scale_free_combined/filtered_conflicts/hubs/Heart-LeftVentricle_hub_exp_iso.txt;/Users/ashissaha/github/gtex_twn/results/quic_scale_free_combined/filtered_conflicts/hubs/Skin-SunExposed_Lowerleg__hub_exp_iso.txt;/Users/ashissaha/github/gtex_twn/results/quic_scale_free_combined/filtered_conflicts/hubs/Skin-NotSunExposed_Suprapubic__hub_exp_iso.txt')
                     default='/scratch1/battle-fs1/ashis/results/gtex_twn/rsem/results_15000/quic_scale_free_combined/selected/filtered_conflicts/hubs/WholeBlood_hub_exp_iso.txt;/scratch1/battle-fs1/ashis/results/gtex_twn/rsem/results_15000/quic_scale_free_combined/selected/filtered_conflicts/hubs/Muscle-Skeletal_hub_exp_iso.txt;/scratch1/battle-fs1/ashis/results/gtex_twn/rsem/results_15000/quic_scale_free_combined/selected/filtered_conflicts/hubs/Skin-SunExposed_Lowerleg__hub_exp_iso.txt;/scratch1/battle-fs1/ashis/results/gtex_twn/rsem/results_15000/quic_scale_free_combined/selected/filtered_conflicts/hubs/Skin-NotSunExposed_Suprapubic__hub_exp_iso.txt;/scratch1/battle-fs1/ashis/results/gtex_twn/rsem/results_15000/quic_scale_free_combined/selected/filtered_conflicts/hubs/Heart-LeftVentricle_hub_exp_iso.txt')
args <- add_argument(args, '-f', '--features',
                     help='feature files, TE/IR data files, semi-colon (;) separated. needed only for feature names',
                     default='/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/rsem/gene_rpkm_15000/WholeBlood.txt;/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/rsem/gene_rpkm_15000/Muscle-Skeletal.txt;/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/rsem/gene_rpkm_15000/Skin-SunExposed_Lowerleg_.txt;/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/rsem/gene_rpkm_15000/Skin-NotSunExposed_Suprapubic_.txt;/scratch0/battle-fs1/GTEx_Analysis_2015-01-12/processedData/data_used_for_twn/rsem/gene_rpkm_15000/Heart-LeftVentricle.txt')
args <- add_argument(args, "-annot",
                     help="NA for TE-TE or TE-IR hubs, transcript annotation file for IR-TE or IR-IR hubs",
                     default=NA)
#default="/scratch1/battle-fs1/ashis/progdata/gtex_hidden_factor/rsem/annot/transcript_annot.txt")
args <- add_argument(args, '-l', '--labels',
                     help='setting labels, semi-colon (;) separated.',
                     default='WholeBlood;Muscle-Skeletal;Skin-SunExposed_Lowerleg_;Skin-NotSunExposed_Suprapubic_;Heart-LeftVentricle')
args <- add_argument(args, '-o', '--output',
                     help='output directory.',
                     default='results/catplots/te_ir')


argv = parse_args(args)
hub_fn_input <- argv$hub
feature_fn_input <- argv$f
label_input <- argv$l
annot_fn <- argv$annot
method <- argv$method
out_prefix <- argv$o

##### process hub file and label input
parts <- strsplit(hub_fn_input, ';')[[1]]
valid_parts <- as.vector(sapply(parts, function(p) nchar(p)>0))
hub_fns <- parts[valid_parts]

parts <- strsplit(label_input, ';')[[1]]
valid_parts <- as.vector(sapply(parts, function(p) nchar(p)>0))
labels <- parts[valid_parts]

parts <- strsplit(feature_fn_input, ';')[[1]]
valid_parts <- as.vector(sapply(parts, function(p) nchar(p)>0))
feature_fns <- parts[valid_parts]

if ((length(hub_fns) != length(labels)) ||(length(hub_fns) != length(feature_fns) )){
  stop('hubs, features, and labels parameters must have same length.')
}

n_settings <- length(hub_fns)

if (n_settings < 2){
  stop('at least 2 hubs needed for overlap analysis.')
}

##### read and process hubs
hub_list <- lapply(hub_fns, FUN = read.table, header = F, sep = '\t', quote = "", col.names=c('hub','deg','connections'), colClasses=c('character', 'integer','character'), comment.char = "", stringsAsFactors = F)

for(i in 1:length(hub_list)){
  hub <- hub_list[[i]]
  hub_list[[i]] <- hub$hub
}

##### read and process features
if(!is.na(annot_fn)){
  annot <- read.table(annot_fn, sep='\t', header=T, quote="", comment.char = "", colClasses = c(chr='character'), stringsAsFactors=F)
  trans2gene <- new.env(hash = T)
  
  tmp <- apply(annot, 1, function(row){
    trans2gene[[row[1]]] <<- row[['sym']]
    return(NA)
  })
}

feature_list <- lapply(feature_fns, function(fn){
  fh <- file(fn)
  line <- readLines(fh, n=1)
  cols <- strsplit(line, '\t')[[1]][-1]
  close(fh)
  if(is.na(annot_fn)){
    return(cols)
  } else{
    genes <- unique(unlist(lapply(cols, function(col) trans2gene[[col]])))
    return(genes)
  }
})

##### create pairwise catpots
overlap_data <- matrix(1, nrow = n_settings, ncol = n_settings)
dummyROCobj = roc(1:5/5, as.factor(c(0,0,1,1,1))) # required for AUC-like computation

for (i in 2:n_settings){
  for (j in 1:(i-1)){
    hubs1 <- hub_list[[i]]
    hubs2 <- hub_list[[j]]
    common_features <- intersect(feature_list[[i]], feature_list[[j]])
    hubs1 <- hubs1[hubs1 %in% common_features]
    hubs2 <- hubs2[hubs2 %in% common_features]
    
    n_max_rank <- min(length(hubs1), length(hubs2))
    prop <- sapply(1:n_max_rank, function(n){
      length(intersect(hubs1[1:n], hubs2[1:n])) / n
    })
    
    n_common_features = length(common_features)
    out_fn <- paste0(out_prefix, '-', labels[i], '-', labels[j],'.png')
    png(out_fn,  width = 800, height = 600)
    plot(c(0, 1:n_max_rank, n_common_features), c(0,prop,1), pch=19, xlab = 'size of hubs', ylab='proportion in common', cex.lab=1.5, cex.axis=1.5)
    lines(c(0, 1:n_max_rank, n_common_features), c(0,prop,1))
    lines(c(0, n_common_features), c(0,n_common_features)/n_common_features, lty=2)
    dev.off()
    
    dummyROCobj$cutoffs = c(0, 1:n_max_rank/n_common_features, 1)
    dummyROCobj$fpr = c(0, 1:n_max_rank/n_common_features,1)
    dummyROCobj$tpr=c(0, prop, 1)
    
    overlap_data[i,j] = overlap_data[j,i] = auc(dummyROCobj)
  }
}

colnames(overlap_data) <- labels
rownames(overlap_data) <- labels

### save overlap data
write.table(overlap_data, file = paste0(out_prefix, '.txt'), sep='\t', col.names = NA, quote = F)

