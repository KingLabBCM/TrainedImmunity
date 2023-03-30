library(Seurat)
library(patchwork)
library(cowplot)
library(tidyverse)
library(scCATCH)

samp.type = 'full.combined'

# Cluster data ------------------------------------------------------------

load('data/full.combined-pre-cc-scaled.Rdata')

# Relevel groups for prettier plotting
full.combined$orig.ident = factor(full.combined$orig.ident, 
                                  levels = c('wt.naive','wt.infected','ko.naive','ko.infected'))

# Perform integrated analysis
DefaultAssay(full.combined) <- "integrated"

# Scale data
full.combined <- ScaleData(full.combined, features = rownames(full.combined))

# Find highly variable features
full.combined <- FindVariableFeatures(full.combined, selection.method = "disp", nfeatures = 4000)

# Run PCA 
full.combined <- RunPCA(full.combined, features = VariableFeatures(full.combined))
DimPlot(full.combined, group.by = 'Phase', reduction = 'pca')
DimPlot(full.combined, group.by = 'orig.ident', reduction = 'pca')

# Elbow plot - approximate number of dimensions
ElbowPlot(full.combined, ndims = 50, reduction = 'pca')
ggsave(filename = paste0('output/plots/elbow-',samp.type,'.png'), device = 'png', width = 4, height = 3)

# UMAP
n_dim=30
full.combined <- RunUMAP(full.combined, dims = 1:n_dim, 
                         metric = 'cosine',
                         angular.rp.forest = TRUE,
                         n.neighbors = 50, min.dist = 0.5, spread = 1.5)
DimPlot(full.combined, reduction = "umap")

# Cluster cells
full.combined <- FindNeighbors(full.combined, dims = 1:n_dim)
full.combined <- FindClusters(full.combined, 
                              method = 'igraph',
                              algorithm = 2,
                              resolution = 0.4)

# Save clustered data
saveRDS(full.combined, file = paste0('data/full.combined-clustered.rds'))

# Plot clusters
DimPlot(full.combined, reduction = "umap")
ggsave(filename = paste0('output/plots/clusters-',samp.type,'.png'), device = 'png',
       height = 5, width=6)

# Visualization
p1 <- DimPlot(full.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(full.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave(filename = paste0('output/plots/clusters-type-',samp.type,'.png'), device = 'png',
       height = 5, width=11)

DimPlot(full.combined, reduction = "umap", split.by = "orig.ident", ncol = 2)
ggsave(filename = paste0('output/plots/clusters-type-faceted-',samp.type,'.png'), device = 'png',
       height = 6, width=7)

# Plot cell type markers --------------------------------------------------

full.combined = readRDS('data/full.combined-clustered.rds')

DefaultAssay(full.combined) <- "RNA"

load('data/Gene_lists_scRNA-seq.Rdata')
type.combs = all.markers %>% 
  select(-Gene) %>% 
  unique()

for (i in 1:nrow(type.combs)) {
  cur.type = type.combs$Type[i]
  cur.subtype = type.combs$Subtype[i]
  mrkrs = all.markers %>% filter(Type == cur.type, Subtype == cur.subtype) %>% pull(Gene) %>% unique()
  
  print(
    DotPlot(full.combined, features = mrkrs, dot.scale = 8) + 
      RotatedAxis() + 
      ggtitle(cur.type, cur.subtype)
  )
  ggsave(filename = paste0('output/plots/dot-gene-lists-',samp.type,'-',cur.type,'-',cur.subtype,'.png'),
         device = 'png', width = 8, height = 10)
}

# Differential expression -------------------------------------------------

full.combined = readRDS('data/full.combined-clustered.rds')
DefaultAssay(full.combined) <- "RNA"

# Find cluster markers
full.combined.markers <- FindAllMarkers(full.combined, 
                                        only.pos = FALSE, 
                                        min.pct = 0.25, 
                                        logfc.threshold = 0.25, 
                                        test.use = 'wilcox')
saveRDS(full.combined.markers, file = paste0('output/all-markers-',samp.type,'.rds'))

mrkrs = full.combined.markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_logFC, n=3)

DoHeatmap(full.combined, features = mrkrs$gene, assay = 'integrated') + NoLegend()
ggsave(filename = paste0('output/plots/markers-heatmap-',samp.type,'.png'), device = 'png', width = 11, height = 7)

# Run scCATCH -------------------------------------------------------------

clu_markers <- findmarkergenes(object = full.combined,
                               species = 'Mouse',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = 'Bone marrow',
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = 'Bone marrow')

# Annotate cell types -----------------------------------------------------

full.combined = readRDS('data/full.combined-clustered.rds')
load('data/Gene_lists_scRNA-seq.Rdata')

all.markers = all.markers %>% group_by(Gene) %>% mutate(occurrences = n())
load('output/annotations-full.combined.Rdata')

de.markers = list()
scc.markers = list()
for (i in 0:21) {
  print(paste('Cluster',i))
  
  # Estimate label by DE results
  df = mrkrs %>% filter(cluster == i & avg_logFC >= 0.7) %>% select(gene, avg_logFC) %>% rename(Gene = gene)
  de.markers[[paste('Cluster',i)]] = all.markers %>% left_join(df) %>% na.omit() %>% arrange(-avg_logFC)
  
  # estimate label by scCATCH results
  df = clu_markers$clu_markers %>% filter(cluster == i) %>% select(gene, avg_logfc) %>% rename(Gene = gene)
  if (nrow(df) > 0) scc.markers[[paste('Cluster',i)]] = all.markers %>% left_join(df) %>% na.omit()  %>% arrange(-avg_logfc)
}
openxlsx::write.xlsx(de.markers, 'output/celltype-labeling-DE.xlsx')
openxlsx::write.xlsx(scc.markers, 'output/celltype-labeling-scCATCH.xlsx')

# Label predicted cell types
full.combined <- RenameIdents(full.combined,
                              `0` = "CD41+ HSC",
                              `1` = "LT-HSC", 
                              `2` = "MPP3/GMP",
                              `3` = "MPP3/GMP", 
                              `4` = "Macrophages monocytes", 
                              `5` = "ST-HSC",
                              `6` = "Neutrophils", 
                              `7` = "Stage I Neutrophil",
                              `8` = "Pre B cell", 
                              `9` = "Macrophages monocytes", 
                              `10` = "Pre B cell",
                              `11` = "Erythrocytes/MkE", 
                              `12` = "Macrophages monocytes", 
                              `13` = "MPP3/GMP", 
                              `14` = "Pre B cell", 
                              `15` = "Neutrophils",  
                              `16` = "B cell", 
                              `17` = "B cell", 
                              `18` = "Stage I Neutrophil", 
                              `19` = "Pro B cell", 
                              `20` = "IA-HSC", 
                              `21` = "NK cell") 

# Create ident variables
full.combined$celltype.group <- paste(Idents(full.combined), full.combined$orig.ident, sep = "_")
full.combined$group = full.combined$orig.ident
full.combined$celltype <- Idents(full.combined)

# Format celltype data for plots
celltype.ord = c("LT-HSC","CD41+ HSC","ST-HSC", "IA-HSC",
                 "Erythrocytes/MkE",
                 "MPP3/GMP", "Macrophages monocytes",
                 "Stage I Neutrophil","Neutrophils",
                 "Pro B cell","Pre B cell","B cell",
                 "NK cell"    
)
full.combined$celltype <- factor(full.combined$celltype, levels = celltype.ord)
Idents(full.combined) <- "celltype"

# Save labeled data
save(full.combined, file = 'data/full.combined-labeled.Rdata')

# Create smaller version of data
fc.lite = DietSeurat(
  full.combined,
  counts = FALSE,
  data = TRUE,
  scale.data = TRUE,
  features = NULL,
  assays = 'integrated',
  dimreducs = 'umap',
  graphs = NULL
)

save(fc.lite, file = 'data/full.combined-labeled_LITE.Rdata')

# Cleanup -----------------------------------------------------------------
rm(list=ls())
