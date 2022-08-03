##### Bradley Jun 2022
##### Cluster merging after DE results
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_env

# Set up
library('Seurat')
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_colon")
load("results/objects/seur_initial_dimred_cluster.Rds")

# Load the DE results
resdir <- "results/tables/pairwiseDEonSeuratClusters/roc"
res <- list.files(resdir)
res <- paste(resdir, res, sep = "/")
roc <- lapply(res, read.table, sep = '\t')
roc_raw <- do.call(rbind, roc)


# Reduce to the DE genes by given cut off
like_smil = T

if(like_smil == T){
        for(x in seq(1:length(roc))){
                roc[[x]]$gene <- rownames(roc[[x]])
                roc[[x]] <-  roc[[x]][roc[[x]]$power > 0.2,]
                }
        } else { for(x in seq(1:length(roc))){
                roc[[x]]$gene <- rownames(roc[[x]])
                roc[[x]] <-  roc[[x]][roc[[x]]$power > 0.45,];
                roc[[x]] <- roc[[x]][(abs(roc[[x]]$pct.1 - roc[[x]]$pct.2) > 0.3),];
        }
}


# Exploring the distribution of total DE genes based on this
counts <- vector("list", length = length(roc))
for(x in seq(1:length(roc))){
        counts[[x]] <- nrow(roc[[x]])
}
counts <- do.call(rbind, counts)


hist(counts, main = "Number of DE genes between clusters after filtering")


cuts <- data.frame(cutoff = seq(from = 0, to = 200, by = 10), nmerge = rep(0, length(seq(from = 0, to = 200, by = 10))))
for(r in seq(1:nrow(cuts))){
        cuts[r,2] <- length(counts[counts > cuts[r,1]])
}

# Generate a matrix, colour by the number of DE genes based on this subsetting
mat <- matrix(nrow = length(levels(factor(seur@meta.data$seurat_clusters))), ncol = length(levels(factor(seur@meta.data$seurat_clusters))))

dimnames = seq(1:(length(levels(factor(seur@meta.data$seurat_clusters)))))
dimnames = dimnames - 1

colnames(mat) <- dimnames
rownames(mat) <- dimnames


# Fill the matrix
roc_df <- do.call(rbind, roc)
roc_df$cond_pair <- rep(0,nrow(roc_df))
for(r in seq(1:nrow(roc_df))){
        if(as.numeric(roc_df$cond1[r]) > as.numeric(roc_df$cond2[r])){
                roc_df$cond_pair[r] = paste(roc_df$cond1[r], roc_df$cond2[r], sep = ".")
        } else {
                roc_df$cond_pair[r] = paste(roc_df$cond2[r], roc_df$cond1[r], sep = ".")
        }
}


for(c in seq(1:ncol(mat))){
        for(r in seq(1:nrow(mat))){
                if(as.numeric(colnames(mat)[c]) == as.numeric(rownames(mat)[r])){
                        mat[r,c] <- NA
                } else {
                        if(as.numeric(colnames(mat)[c]) > as.numeric(rownames(mat)[r])) {
                                mat[r,c] <- nrow(roc_df[roc_df$cond1 == colnames(mat)[c] & roc_df$cond2 == rownames(mat)[r],]);
                        } else {
                                mat[r,c] <- nrow(roc_df[roc_df$cond1 == colnames(mat)[r] & roc_df$cond2 == rownames(mat)[c],]);
                        }
                }
        }
}

# Using a cut off of 30 DEGs at this threshold - we only merge clusters 0 and 8.
seur@meta.data$seurat_clusters <- as.character(seur@meta.data$seurat_clusters)
seur@meta.data$seurat_clusters[seur@meta.data$seurat_clusters == "8"] <- "0"

# Recalculate the markers
new_markers <- FindAllMarkers(seur, assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the new markers and the object
write.csv(new_markers, "results/tables/cluster_markers_after_merging.csv")
write.csv(new_markers[new_markers$cluster == "15",], "results/tables/cluster15_markers_after_merging.csv")
save(seur, file="results/objects/seur_dimred_post_merging.Rds")
