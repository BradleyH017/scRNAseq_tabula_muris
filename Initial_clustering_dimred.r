##### Bradley Jun 2022
##### Initial clustering and dimensionality reduction of scRNAseq from tabula muris
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_env

# Set up
library('Seurat')
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_colon")
load("data/Preprocessing/raw_seur_pre-processed.Rds")
seur <- seur_final
rm(seur_final)

# Initial dimred without batch correction to Check
# Not regressing MT% as doesn't seem to be an issue
seur <- SCTransform(seur, verbose = FALSE)
seur <- FindVariableFeatures(seur)
seur <- RunPCA(seur, verbose = FALSE)
# Check loadings/scree etc
library(patchwork)
pdf(file="results/plots/PCA_QC.pdf")
DimPlot(seur, reduction="pca") + ElbowPlot(seur) + VizDimLoadings(seur, dims = 1:2, reduction = "pca")
dev.off()
# Now run UMAP
seur <- RunUMAP(seur, reduction = "pca", dims = 1:20)
DimPlot(seur, reduction="umap")

# now run clustering (requires different r-matrix package)
# conda install -c conda-forge r-matrix=1.3.2
seur <- FindNeighbors(seur, dims = 1:20, reduction="pca", k.param=20)
seur <- FindClusters(seur, verbose = FALSE)
cl <- DimPlot(seur, label = TRUE) + NoLegend()
smp <- DimPlot(seur, group.by="Sample",label = TRUE) + NoLegend()

seur@meta.data$Sex <- unlist(strsplit(seur@meta.data$Sample, "\\_"))[c(F,F,T)]
sx <- DimPlot(seur, group.by="Sex", label = TRUE) + NoLegend()
pdf(file="results/plots/UMAP_samples_clusters.pdf")
cl + smp + sx
dev.off()


# Looking at our C11orf compound vs WT DE results on the UMAP
de <- c("Colca2","Sh2d7" , "Trpm5" , "C11orf53" ,"Sh2d6","Pou2f3","Gnat3","Itprid1","Ccdc129",
          "Zfp579","Avil","Rgs13","Igkv4-91","Alox5","Fgf7","Trim38","St18","Cirbp","Pik3r5",
          "Nebl","Pcsk4","Igkv3-1","Tgfbr3l","Hmx2","Ly6g6f","Adgrg1","Aspn","Krt12","Izumo4",
          "Spib","Zfp692","Cav2","2410002F23Rik","Ccdc9","Samd91","Iglv1","IL25")
de <- de[de %in% rownames(seur@assays$RNA)]
plots <- vector("list", length = length(de))
for(g in seq_along(de)){
  p <- FeaturePlot(seur,features=de[g]) + NoLegend()
  plots[[g]] <- p
}
pdf(file="results/plots/compound_de_umap.pdf", width=20, height=16)
wrap_plots(plots, nrow = 5, ncol = 6)
dev.off()

pdf(file="results/plots/cis_umap.pdf", width=20, height=16)
FeaturePlot(seur, features=cis)
dev.off()

# Before calculating the markers of each cluster, generate the log2p10k Matrix
counts_list <- vector("list", length = ceiling(ncol(seur@assays$RNA@counts)/1000))
for(x in seq(1:(length(counts_list)-1))){
       counts_list[[x]] <- as.data.frame(seur@assays$RNA@counts[,(((1000*x)-999):(1000*x))])
}
counts_list[[length(counts_list)]] <- seur@assays$RNA@counts[,((floor(ncol(seur@assays$RNA@counts)/1000)*1000)+1):ncol(seur@assays$RNA@counts)]

# Now perform the logcpm normalisation of each 1000. These are too large to do this in R, so will save and do this from the command line
library('edgeR')
for(x in seq(1:length(counts_list))){
        y <- DGEList(counts = counts_list[[x]]);
        cpm <- edgeR::cpm(y, log=F, normalized.lib.sizes=F);
        cpm <- cpm/1000; #So that is per 10k
        logcpm <- log2(cpm +1);
        counts_list[[x]] <- logcpm;
        counts_list[[x]] <- as(counts_list[[x]], "sparseMatrix");
}
log2TPM <- do.call(cbind, counts_list)
seur[['log2TP10K']] <- CreateAssayObject(counts = log2TPM)
# Performing marker analysis on this doesn't seem to perform too well. We get genes like MALAT1 and spike ins as cluster markers.
# May be because of the different technology (plate not droplet), PCR step etc
# From: https://github.com/satijalab/seurat/issues/2180 . It is suggested to use log-normalised dataset
# Can calculate this with seurat
DefaultAssay(seur) <- "RNA"
seur <- NormalizeData(seur, verbose = FALSE, normalization.method = "LogNormalize")
markers <- FindAllMarkers(seur, assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Calculate the markers of each cluster
DefaultAssay(seur) <- "RNA"
markers <- FindAllMarkers(seur, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(markers, "results/tables/cluster_markers.csv")
write.csv(markers[markers$cluster == 15,], "results/tables/cluster15_markers.csv")


# save
save(seur, file = "results/objects/seur_initial_dimred_cluster.Rds")


# Save a list of comparisons for DE across these (merging)
clusters <- levels(factor(seur@meta.data$seurat_clusters))
comparisons <- expand.grid(clusters, clusters)
comparisons <- comparisons[comparisons$Var1 != comparisons$Var2,]
comparisons$comb <- ""
for(r in seq(1:nrow(comparisons))){
       if(as.numeric(comparisons$Var1)[r] < as.numeric(comparisons$Var2)[r]){
               comparisons[r, 3] = paste(comparisons$Var1[r], comparisons$Var2[r], sep = ".")
       } else {
               comparisons[r, 3] = paste(comparisons$Var2[r], comparisons$Var1[r],sep = ".")
       }
}

comp <- as.data.frame(comparisons)
comp <- comp[!duplicated(comp$comb),-3]

# To submi this as an array job, want to do this once per test. So save the list as seperate files for each condition and pass it as an argument
write.csv(comp[,1], "results/tables/comparisons_for_MAST_test_using_seurat_clusters_cond1.txt", row.names = F, quote=F)
write.csv(comp[,2], "results/tables/comparisons_for_MAST_test_using_seurat_clusters_cond2.txt", row.names = F, quote=F)
## Need to go and remove the header for these > *_final.txt
