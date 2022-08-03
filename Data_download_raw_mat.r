# Bradley May 2022
# Looking at tabula muris scRNAseq
# Downaloded using wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109774/suppl/GSE109774_Colon.tar.gz

# Set up
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)

# looking at the data initially
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_colon/data/Colon")

# Make the dataframe
names <- list.files()
flist <- vector("list", length = length(names))
for(f in seq_along(flist)){
  flist[[f]] <- read.csv(names[f], header=F, row.names=1)
  colnames(flist[[f]]) <- gsub(".csv", "", names[f])
}

raw <- do.call(cbind, flist)
write.csv(raw, "../raw_counts.csv")

# Check interseting genes
int <- c("1810046K07Rik", "2010007H06Rik", "Gm684", "Trpm5", "Dclk1", "Pou2f3")

raw_int <- as.data.frame(t(raw))
raw_int <- raw_int[,colnames(raw_int) %in% int]
