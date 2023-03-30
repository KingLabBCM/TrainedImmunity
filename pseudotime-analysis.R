library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(ggbeeswarm)
library(ggthemes)
library(pals)
library(slingshot)
library(destiny)
library(gam)

# Setup -------------------------------------------------------------------

## Load data
load('data/full.combined-labeled.Rdata')

## Reduce to WT experiments
fc.wt = subset(full.combined, subset = group %in% c('wt.naive','wt.infected'))

## Transform to SingleCellExperiment object
sce.fc = as.SingleCellExperiment(fc.wt, assay = 'integrated')
rm(full.combined)

col.pal.celltype = c(RColorBrewer::brewer.pal(12,"Paired"),'grey')[c(1:3,6,5,4,7:13)]
col.pal.clusters = cols25(22)

mydir = 'output/plots/pseudotime'


# Slingshot pseudotime ----------------------------------------------------

dim.red = 'PCA'

# First look at Slingshot lineages
lin1 <- getLineages(sce.fc, clusterLabels = 'celltype',
                    start.clus = 'LT-HSC',
                    reducedDim = dim.red)
png(paste0(mydir, "/slingshot_lineage_celltype_",dim.red,".png"))
plot(reducedDims(sce.fc)[[dim.red]], col = col.pal.celltype[as.numeric(sce.fc$celltype)], asp = 1, pch = 16,
     main = 'Pseudotime lineages (Slingshot)')
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
legend(
  "bottomleft",
  legend=levels(sce.fc$celltype),
  col=col.pal.celltype,
  cex=0.7,pch=19, ncol=1)
dev.off()

## Slingshot using cell type
colData(sce.fc)$celltype_chr <- as.character(colData(sce.fc)$celltype)
colData(sce.fc)$celltype_chr[which(colData(sce.fc)$celltype_chr == 'Stage I Neutrophil')] <- 'Neutrophils'
sce.fc <- slingshot(sce.fc, clusterLabels = 'celltype_chr',
                    reducedDim = dim.red,
                    start.clus = 'LT-HSC',
                    end.clus = c('Bcell','Macrophages monocytes','Neutrophils','NK cell', 'Erythrocytes/MKE'),
                    approx_points = 100)
save(sce.fc, file = paste0("output/pseudotime-slingshot-celltype_",dim.red,".Rdata"))

# Load saved Slingshot output
load(paste0("output/pseudotime-slingshot-celltype_",dim.red,".Rdata"))
sce.fc@colData = sce.fc@colData[,-(16:22)] # remove previous results -- this file saved both
celltypes = levels(sce.fc$celltype)
sce.fc$celltype <- factor(sce.fc$celltype, labels = c(celltypes[1:3],'IA-HSCs',celltypes[5:13]))
celltypes = levels(sce.fc$celltype)

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- colorRampPalette(RColorBrewer::brewer.pal(11,'Spectral')[-6])(100)
sling.pt = sce.fc@colData[,grep('^slingPseudotime', colnames(sce.fc@colData), value = TRUE)]
pdf(paste0(mydir, "/slingshot_avg_pseudotime_celltype_",dim.red,".pdf"))
plot(reducedDims(sce.fc)[[dim.red]], col = colors[cut(apply(sling.pt,1,function(x) mean(x, na.rm = TRUE)),breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce.fc), lwd=2, col = 'black')
dev.off()

# Plot UMAP colored by Slingshot pseudotime
pdf(paste0(mydir, "/slingshot_avg_pseudotime_celltype_PCA-to-UMAP.pdf"))
plot(reducedDims(sce.fc)[['UMAP']], col = colors[cut(apply(sling.pt,1,function(x) mean(x, na.rm = TRUE)),breaks=100)], pch=16, asp = 1)
dev.off()

# Plot PC1 vs PC2 colored by cell type
pdf(paste0(mydir, "/slingshot_pseudotime_celltype_by-celltype.pdf"))
plot(reducedDims(sce.fc)[[dim.red]], col = col.pal.celltype[sce.fc$celltype], pch=16, asp = 1)
points(reducedDims(sce.fc)[[dim.red]][which(sce.fc$celltype == 'IA-HSCs'),], pch=16, asp = 1, col = col.pal.celltype[4])
lines(SlingshotDataSet(sce.fc), lwd=2, col = 'black')
legend(
    "bottomleft",
    legend=levels(sce.fc$celltype),
    col=col.pal.celltype,
    cex=0.7,pch=19, ncol=1)
dev.off()
