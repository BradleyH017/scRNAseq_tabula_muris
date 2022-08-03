##### Bradley Jun 2022
##### DE between clusters to test merging
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_env

# Set up
library('Seurat')
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_colon")
load("results/objects/seur_initial_dimred_cluster.Rds")

# Now set up the options for DE
library(optparse)
option_list = list(
       make_option(c("-c", "--cond1"), action = "store", default = NA, type ="character",
               help="What is the first cluster to compare"),
        make_option(c("-d", "--cond2"), action = "store", default = NA, type ="character",
                help="What is the second cluster to compare")
       )
opt = parse_args(OptionParser(option_list=option_list))

print("Option for cluster 1 is")
opt$cond1
print("Option for cluster 2 is")
opt$cond2


print("Doing MAST")
DefaultAssay(seur) <- "RNA"
test <- FindMarkers(seur, ident.1 = opt$cond1, ident.2 = opt$cond2, test.use = "MAST", slot="data")
test$cond1 <- rep(opt$cond1, nrow(test))
test$cond2 <- rep(opt$cond2, nrow(test))
head(test)
write.table(test, paste("results/tables/pairwiseDEonSeuratClusters/MAST/MAST_", opt$cond1, "vs", opt$cond2, ".txt",sep = ""), sep = "\t")
print("Done MAST")

print("Doing ROC")
DefaultAssay(seur) <- "log2TP10K"
test <- FindMarkers(seur, ident.1 = opt$cond1, ident.2 = opt$cond2, test.use = "roc", min.pct = 0.1, logfc.threshold = 0.25, slot="data")
test$cond1 <- rep(opt$cond1, nrow(test))
test$cond2 <- rep(opt$cond2, nrow(test))
head(test)
write.table(test, paste("results/tables/pairwiseDEonSeuratClusters/roc/roc_", opt$cond1, "vs", opt$cond2, ".txt", sep = ""), sep = "\t")
print("Done ROC")
