##### Bradley Jun 2022
##### Interesting cluster analysis
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_env

# Set up
library('Seurat')
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_colon")
load("results/objects/seur_dimred_post_merging.Rds")

# Plot C53 and Colca2 nice for Thesis
# C53 = "1810046K07Rik"
# Colca2 = "Gm684"
library(ggplot2)
p1 <- FeaturePlot(seur, features = c("1810046K07Rik")) + ggtitle("C11orf53")
p2 <- FeaturePlot(seur, features = c("Gm684")) + ggtitle("Colca2")
library(patchwork)
ppi=300
png(file="../Fig_plots/Thesis/Chapter3/C53_Colca2_umap.png", width=5*ppi, height=8*ppi, res=ppi)
p1/p2
dev.off()

##### Annotating clusters by enrichment of markers (external dataset)
##### Went to pangola db and searched for tuft cells
##### Picked the one which used mouse colon (SRA653146)
generate_external_markers = FALSE
if(generate_external_markers == T) {
  library(Matrix)
  load("data/Pangola_gene_sets/SRA653146_SRS2874272.sparse.RData")
  cell_anno <- read.table("data/Pangola_gene_sets/SRA653146_SRS2874272.clusters.txt", sep = " ")
  library(Seurat)
  pg_seur <- CreateSeuratObject(counts=sm)
  cell_anno_add <- cell_anno[,2]
  names(cell_anno_add) <- cell_anno[,1]
  pg_seur <- AddMetaData(pg_seur, cell_anno_add, col.name="cluster")
  DefaultAssay(seur) <- "RNA"
  pg_seur <- NormalizeData(pg_seur, verbose = FALSE, normalization.method = "LogNormalize")
  # Generalise to annotation - see https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA653146&srs=SRS2874272
  pg_seur@meta.data$gen_cluster <- ifelse(pg_seur@meta.data$cluster %in% c(0,3,4,6,7,8,13,16,17), "Enterocyte", "")
  pg_seur@meta.data$gen_cluster <- ifelse(pg_seur@meta.data$cluster %in% c(1,2), "Epithelial_cell", pg_seur@meta.data$gen_cluster)
  pg_seur@meta.data$gen_cluster <- ifelse(pg_seur@meta.data$cluster %in% c(5,9,10,11,12,15,19,20,21), "Goblet", pg_seur@meta.data$gen_cluster)
  pg_seur@meta.data$gen_cluster <- ifelse(pg_seur@meta.data$cluster %in% c(14), "Tuft", pg_seur@meta.data$gen_cluster)
  pg_seur@meta.data$gen_cluster <- ifelse(pg_seur@meta.data$cluster %in% c(18,22), "Enteroendocrine", pg_seur@meta.data$gen_cluster)
  # Subset non anno cells
  good_cells <- rownames(pg_seur@meta.data[!is.na(pg_seur@meta.data$cluster),])
  pg_seur <- subset(pg_seur, cells = good_cells)
  Idents(pg_seur) <- "gen_cluster"
  pg_markers <- FindAllMarkers(pg_seur, assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(pg_markers, "data/Pangola_gene_sets/general_cluster_markers.csv")
} else {
  pg_markers <- read.csv("data/Pangola_gene_sets/general_cluster_markers.csv", row.names=1)
  pg_markers <- pg_markers[pg_markers$p_val_adj < 0.05,]
}

# reformat for enrichment
for(r in 1:nrow(pg_markers)){
  if(length(unlist(strsplit(pg_markers$gene[r], "\\-"))) == 2){
    pg_markers$gene[r] <- unlist(strsplit(pg_markers$gene[r], "\\-"))[c(T,F)]
  } else {
    temp <- unlist(strsplit(pg_markers$gene[r], "\\-"))[c(T,T,F)]
    pg_markers$gene[r] <- paste(temp, collapse = "-")
  }
}

pg_clusters <- levels(factor(pg_markers$cluster))
pg_sigs <- vector("list", length = length(pg_clusters))
names(pg_sigs) <- pg_clusters
for(c in seq_along(pg_clusters)){
  pg_sigs[[c]] <- pg_markers[pg_markers$cluster == pg_clusters[c],]$gene
}

# Load in the cluster markers and generate the gene List
markers_all <- read.csv("results/tables/cluster_markers_after_merging.csv", row.names=1)

# Generate the geneList list to loop through the markers
markers <- vector("list", length = length(levels(factor(markers_all$cluster))))
names(markers) <- levels(factor(markers_all$cluster))
for(x in seq(1:length(markers))){
        markers[[x]] <- markers_all[markers_all$cluster == names(markers)[x],]
}

# Now turn this into a ranked gene list (by logFC)
geneList <- vector("list", length = length(markers))
names(geneList) <- names(markers)
for(x in seq(1:length(markers))){
        geneList[[x]] <- markers[[x]]$avg_logFC;
        names(geneList[[x]]) <- markers[[x]]$gene
}

# Now do enrichment for FDR and nom sig and manually BH correct
library(fgsea)
library(ggplot2)
fgsea_list <- vector("list", length= length(geneList))
names(fgsea_list) <- names(geneList)
for(x in seq(1:length(fgsea_list))){
  # Get very small pvals
        fgsea_list[[x]] <- fgseaMultilevel(pg_sigs, geneList[[x]], minSize = 0, maxSize = 500, scoreType = "pos", eps=0)
        fgsea_list[[x]]$cluster <- rep(names(geneList)[x], nrow(fgsea_list[[x]]))
}

fgseaRes <- do.call(rbind, fgsea_list)
# Manually adjusting p-vals
fgseaRes$BH_padj <- p.adjust(fgseaRes$pval, method="fdr")
fgseaRes <- fgseaRes[fgseaRes$BH_padj < 0.05,]
fgseaRes <- fgseaRes[,-8]
write.csv(fgseaRes, "results/tables/enrichment_results/fgsea_markers_vs_pg_database_anno.csv")

# Re-annotate the clusters using these (only most significant, and then by proximity)
fgseaRes <- fgseaRes[fgseaRes$BH_padj < 0.01,]
seur@meta.data$gen_cluster <- ifelse(seur@meta.data$seurat_clusters %in% c(1,2,3,5,7,9,16), "Enterocyte", seur@meta.data$seurat_clusters)
seur@meta.data$gen_cluster <- ifelse(seur@meta.data$seurat_clusters %in% c(4,6,10,12,13), "Goblet", seur@meta.data$gen_cluster)
seur@meta.data$gen_cluster <- ifelse(seur@meta.data$seurat_clusters %in% c(15), "Tuft", seur@meta.data$gen_cluster)
seur@meta.data$gen_cluster <- ifelse(seur@meta.data$seurat_clusters %in% c(0), "Epithelial", seur@meta.data$gen_cluster)
seur@meta.data$gen_cluster <- ifelse(seur@meta.data$seurat_clusters %in% c(11), "Enterocyte", seur@meta.data$gen_cluster)

# Regerenerate the markers of my clusters (now annotated to Tuft etc)
Idents(seur) <- "gen_cluster"
gen_my_markers = F
if(gen_my_markers == T){
  markers_all <- FindAllMarkers(seur, slot="data", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
  write.csv(markers_all, "results/tables/gen_cluster_markers.csv")
} else {
  markers_all <- read.csv("results/tables/gen_cluster_markers.csv", row.names=1)
}
markers <- vector("list", length = length(levels(factor(markers_all$cluster))))
names(markers) <- levels(factor(markers_all$cluster))
for(x in seq(1:length(markers))){
        markers[[x]] <- markers_all[markers_all$cluster == names(markers)[x],]
}
geneList <- vector("list", length = length(markers))
names(geneList) <- names(markers)
for(x in seq(1:length(markers))){
        geneList[[x]] <- markers[[x]]$avg_logFC;
        names(geneList[[x]]) <- markers[[x]]$gene
}

######## Now do the enrichment against the interesting clusters
# Load the C11orf53 correlated genes from the C53 and WT WGCNA
gi <- read.csv("../Nextflow_alignment/WGCNA/C53_wt_MAD_5000/tables/geneInfo.csv")
gi_cor <- gi[gi$GS.C11orf53 > 0.7 & gi$p.GS.C11orf53 < 0.01,]

# Turn this into a list for enrichment
path_list <- list()
path_list[['C11orf53_correlated_genes_WGCNA']] <- gi_cor$X

# Add the compound DE results
path_list[['C11orfCompound_DE_genes']] <- c("Gm684","Sh2d7","Trpm5","1810046K07Rik","Sh2d6","Pou2f3","Gnat3","Ccdc129","Itprid1",
  "Zfp579","Avil","Rgs13","Igkv4-91","Alox5","Fgf7","Trim38","St18","Cirbp","Pik3r5","Nebl","Pcsk4","Igkv3-1","Tgfbr3l","Hmx2"
  ,"Ly6g6f","Adgrg1","Aspn","Krt12","Izumo4","Spib","Zfp692","Cav2","2410002F23Rik","Ccdc9","Samd91","Iglv1","IL25")


# Now do enrichment for FDR and nom sig and manually BH correct
library(fgsea)
library(ggplot2)
fgsea_list <- vector("list", length= length(geneList))
names(fgsea_list) <- names(geneList)
for(x in seq(1:length(fgsea_list))){
        fgsea_list[[x]] <- fgseaMultilevel(path_list, geneList[[x]], minSize = 0, maxSize = 500, scoreType = "pos")
        fgsea_list[[x]]$cluster <- rep(names(geneList)[x], nrow(fgsea_list[[x]]))
}

fgseaRes <- do.call(rbind, fgsea_list)
# Manually adjusting p-vals
fgseaRes$BH_padj <- p.adjust(fgseaRes$pval, method="fdr")
fgseaRes <- fgseaRes[fgseaRes$BH_padj < 0.05,]

enr_plots <- vector("list", length = length(path_list))
enr_labels <- vector("list", length = length(path_list))
for(p in seq_along(path_list)){
  enr_res = fgseaRes[fgseaRes$pathway == names(path_list)[p],]
  enr_label = paste0("NES=", signif(enr_res$NES, 2), "\n", "FDR=", as.numeric(signif(enr_res$BH_padj,2)))
  temp <- plotEnrichment(path_list[[p]], geneList[['Tuft']]) + theme_bw() + ggtitle(paste("Enrichment of", names(path_list)[p])) +
                          xlab("Rank") + ylab("Enrichment score in Tuft cell markers") + theme(plot.title=element_text(face="bold", size=14)) +
                          annotate(geom="text", x=270, y=0.50, label=enr_label, color="black", size = 4, hjust=1)
  enr_plots[[p]] <- temp
}

# Visualing this at the single cell level
# conda install -c bioconda bioconductor-escape
library(escape)
GS <- c(GSEABase::GeneSet(path_list[[1]]), GSEABase::GeneSet(path_list[[2]]))
GS[[1]]@setName <- names(path_list)[1]
GS[[2]]@setName <- names(path_list)[2]
ES.seurat <- enrichIt(obj = seur,
                      gene.sets = GS,
                      groups = 1000, cores = 2)
seur <- Seurat::AddMetaData(seur, ES.seurat)


# PLotting this
# seurat colour palette is horrible and misleading as bins
# USing archr palette
# Sys.setenv(CONDA_BUILD_SYSROOT="/")
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
library(magrittr)
library(RColorBrewer)
library(purrr)

# Extract archR palette
col.ls <- ArchRPalettes
to_plot <- c("C11orf53_correlated_genes_WGCNA", "C11orfCompound_DE_genes")
cell_plots <- vector("list", length = 2)
for(p in seq_along(to_plot)){
  temp <- FeaturePlot(seur, features = to_plot[p], cols=col.ls[['horizonExtra']]) +
    theme(panel.border = element_rect(color = "black", fill = NA,  size = 0.5)) + ggtitle("")
  cell_plots[[p]] <- temp
}

library(patchwork)
ppi=300
png(file="results/plots/C11orf_signature_enrichment_fgsea_umap.png", width=15*ppi, height=8*ppi, res=ppi)
DimPlot(seur, group.by="gen_cluster", label=T) + NoLegend() + (enr_plots[[1]] + cell_plots[[1]]) / (enr_plots[[2]] + cell_plots[[2]])
dev.off()

# Save this for thesis too
png(file="../Fig_plots/Thesis/Chapter3/C11orf_signature_enrichment_fgsea_umap.png", width=15*ppi, height=8*ppi, res=ppi)
DimPlot(seur, group.by="gen_cluster", label=T) + NoLegend() + ggtitle("") + (enr_plots[[1]] + cell_plots[[1]]) / (enr_plots[[2]] + cell_plots[[2]])
dev.off()


# Also do UMAPs for C11orf53 and COLCA2 only
png(file="../Fig_plots/Thesis/Chapter3/C53_Colca2_umap.png", width=9*ppi, height=4*ppi, res=ppi)
FeaturePlot(seur, features="1810046K07Rik") + ggtitle("C11orf53") + (FeaturePlot(seur, features="Gm684") + ggtitle("Colca2")) + plot_layout(ncol=2)
dev.off()


#### SCRAP


# Extract the horizon Extra, rename the cols
test_cols <- col.ls[['horizonExtra']]
#test_cols <- test_cols[order(as.numeric(names(test_cols)))]
max_range_to_plot = range(sapply(seur@meta.data[,c("C11orf53_correlated_genes_WGCNA", "C11orfCompound_DE_genes")], range))
# Divide the range into specific groups, so matches the actual NES vals when plotted
spread <- max(max_range_to_plot)-min(max_range_to_plot)
names(test_cols) <- signif(seq(from=min(max_range_to_plot), to = max(max_range_to_plot), by=spread/(length(test_cols)-1)),2)


FeaturePlot(seur, features = c("C11orf53_correlated_genes_WGCNA", "C11orfCompound_DE_genes"), cols=test_cols)





FeaturePlot(seur, features = c("C11orf53_correlated_genes_WGCNA", "C11orfCompound_DE_genes"), cols=CustomPalette(low="blue", high="red", mid="white"))

p.ls <- col.ls %>% imap( ~ {
  FeaturePlot(object = seur, features = c("C11orf53_correlated_genes_WGCNA", "C11orfCompound_DE_genes"))
})


FeaturePlot_scCustom(seur, reduction="umap", features = c("C11orf53_correlated_genes_WGCNA", "C11orfCompound_DE_genes"))

&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


horizonExtra
  1         4         6         8         3         9         7         5
"#000436" "#021EA9" "#1632FB" "#6E34FC" "#C732D5" "#FD619D" "#FF9965" "#FFD32B"
  2
"#FFFC5A"
