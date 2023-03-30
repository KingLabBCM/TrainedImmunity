library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)


exp.groups <- c('wt.naive',
               'wt.infected',
               'ko.naive',
               'ko.infected')

out.dir = 'output/qc/'

# Compare QC stats across groups ------------------------------------------

# load data
wt.naive <- readRDS('data/wt.naive-base.rds')
wt.infected <- readRDS('data/wt.infected-base.rds')
ko.naive <- readRDS('data/ko.naive-base.rds')
ko.infected <- readRDS('data/ko.infected-base.rds')

### Compare QC stats

# Create metadata dataframe
metadata <- rbind(wt.naive@meta.data, wt.infected@meta.data,
                  ko.naive@meta.data, ko.infected@meta.data)

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(sample = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA) %>% 
  separate(sample, into = c('type','treatment'), sep = '\\.', remove = FALSE) %>% 
  mutate(type = factor(type, levels = c('wt','ko'))) %>% 
  mutate(treatment = factor(treatment, levels = c('naive','infected')))

metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") +
  scale_fill_discrete('')
ggsave(filename = paste0(out.dir, 'barplot-ncells.png'), height = 3, width = 4)

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)  +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-nUMIs.png'), height = 3, width = 5)


metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  facet_wrap(~ type) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-nUMIs-byType.png'), height = 3, width = 6)


# Visualize the distribution of genes detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-nGenes.png'), height = 3, width = 5)

metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500) +
  facet_wrap(~ type) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-nGenes-byType.png'), height = 3, width = 6)

metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500) +
  facet_wrap(~ treatment) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-nGenes-byTmt.png'), height = 3, width = 6)


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt / 100)) + 
  geom_point() + 
  scale_colour_gradient('Percent mt',low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample)
ggsave(filename = paste0(out.dir, 'points-nGenes-vs-nUMIs.png'), height = 5, width = 6)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-pctmt.png'), height = 3, width = 5)

metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  xlim(c(0,11)) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-pctmt-red.png'), height = 3, width = 5)


metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  xlim(c(0,11)) +
  facet_wrap(~ type)  +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-pctmt-red-byType.png'), height = 3, width = 6)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-complexity.png'), height = 3, width = 5)

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  xlim(c(0.78, 0.97)) +
  geom_vline(xintercept = 0.8) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-complexity-red.png'), height = 3, width = 5)

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  xlim(c(0.78, 0.97)) +
  geom_vline(xintercept = 0.8) +
  facet_wrap(~ type) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-complexity-red-byType.png'), height = 3, width = 6)

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  xlim(c(0.78, 0.97)) +
  geom_vline(xintercept = 0.8) +
  facet_wrap(~ treatment) +
  scale_fill_discrete('') + scale_color_discrete('')
ggsave(filename = paste0(out.dir, 'density-complexity-red-byTmt.png'), height = 3, width = 6)


# Verify cell cycle scaling results ---------------------------------------

load('data/m.cc.genes.Rdata')

find.genes.pca <- function(gene.list, pcs = 1:10, top_n = 10) {
  all.genes <- c()
  for (pc in pcs) {
    top.pos = names(sort(Loadings(grp.obj[['pca']])[,pc], decreasing = TRUE)[1:top_n])
    top.neg = names(sort(Loadings(grp.obj[['pca']])[,pc])[1:top_n])
    
    genes = gene.list[which(gene.list %in% c(top.pos, top.neg))]
    if (length(genes)) cat(paste0('PC',pc,': ', paste(genes,collapse=', ')),'\n')
    all.genes = c(all.genes, genes)
  }
  return(all.genes)
}

### Plot default scaled data

for (grp in exp.groups) {
  # Load data
  grp.obj <- readRDS(paste0('data/',grp,'-filtered-unscaled.rds'))
  
  # Use default linear scaling
  grp.obj <- ScaleData(grp.obj, features = rownames(grp.obj))
  grp.obj <- RunPCA(grp.obj, features = VariableFeatures(grp.obj))
  
  # Calculate cell cycle scores and phase assignments
  grp.obj <- CellCycleScoring(grp.obj, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
  
  # Find genes associated with cell cycle in PCs
  pca.genes <- find.genes.pca(m.cc.genes)

  # Visualize the distribution of cell cycle markers across
  RidgePlot(grp.obj, features = pca.genes, ncol = 4)
  ggsave(filename = paste0(out.dir,'cc-markers-',grp,'.png'), height = 10, width = 11)
  
  # Run PCA on cell cycle genes to verify cells separate entirely by phase
  grp.obj <- RunPCA(grp.obj, features = m.cc.genes)
  print(DimPlot(grp.obj))
  ggsave(filename = paste0(out.dir,'default-scaled-',grp,'.png'), height = 5, width = 6)
  
  ### Plot cell cycle scaled data
  
  # Load data
  grp.obj <- readRDS(paste0('data/',grp,'-filtered-cc-scaled.rds'))
  
  # Plot cell cycles
  grp.obj <- RunPCA(grp.obj, features = m.cc.genes)
  print(DimPlot(grp.obj, group.by = 'Phase'))
  ggsave(filename = paste0(out.dir,'cc-scaled-',grp,'.png'), height = 5, width = 6)
  
  ### Plot cell cycle diff scaled data
  
  # Load data
  grp.obj <- readRDS(paste0('data/',grp,'-filtered-cc-diff-scaled.rds'))
  
  # Plot cell cycles
  grp.obj <- RunPCA(grp.obj, features = m.cc.genes)
  print(DimPlot(grp.obj, group.by = 'Phase'))
  ggsave(filename = paste0(out.dir,'cc-diff-scaled-',grp,'.png'), height = 5, width = 6)
}





