library(Seurat)
library(patchwork)
library(cowplot)
library(tidyverse)

# Set up data for plotting ------------------------------------------------

# Load labeled data
load('data/full.combined-labeled_LITE.Rdata')

# Create ident variables
fc.lite$celltype.group <- paste(Idents(fc.lite), fc.lite$orig.ident, sep = "_")
fc.lite$group = fc.lite$orig.ident
fc.lite$celltype <- Idents(fc.lite)

# Format celltype data for plots
celltype.ord = c("LT-HSC","CD41+ HSC","ST-HSC", "IA-HSC",
                 "Erythrocytes/MkE",
                 "MPP3/GMP", "Macrophages monocytes",
                 "Stage I Neutrophil","Neutrophils",
                 "Pro B cell","Pre B cell","B cell",
                 "NK cell"
)
fc.lite$celltype <- factor(fc.lite$celltype, levels = celltype.ord)
Idents(fc.lite) <- "celltype"

# Filter to WT only
fc.lite.wt <- subset(fc.lite, subset = orig.ident %in% c('wt.naive','wt.infected'))
fc.lite.wt$orig.ident <- factor(fc.lite.wt$orig.ident, labels = c('WT naive', 'WT infected'))

celltypes = levels(fc.lite.wt$celltype)

plt.genes = c('Batf2', 'Cxcl9', 'Arhgap9','Ccl5', 'Nr4a1','Tnfrsf1a')


# Set up raincloud plot function =-----------------------------------------
library(raincloudplots)
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}


# Raincloud plots ---------------------------------------------------------

# Format data
rc.df = data.frame(orig.ident = fc.lite.wt$orig.ident, 
                 celltype = fc.lite.wt$celltype, 
                 t(as.matrix(fc.lite.wt@assays$integrated@data))) 

# Set plotting parameters
par.adj = 1
par.bw = 0.1

# Plot selected genes
for (g in plt.genes) {
  print(g)
  
  # HSCs
  p1 = ggplot(rc.df %>% filter(celltype %in% celltypes[1:3]),
              aes(x=orig.ident,y=get(g), fill = orig.ident))+
    geom_flat_violin(position = position_nudge(x = .2, y = 0),
                     adjust=par.adj, bw=par.bw, scale = 'width', alpha=0.6) +
    geom_point(aes(color = orig.ident),
               position = position_jitter(width = .15, height = 0),
               size = 1, alpha = 0.6, show.legend = FALSE) +
    facet_wrap(~celltype, nrow = 1, strip.position = 'bottom',
               labeller = labeller(celltype = label_wrap_gen(width=15, multi_line = TRUE))) +
    ylab('Expression')+xlab('Cell type')+
    theme_cowplot(font_size = 12)+
    guides(color = FALSE)+
    scale_color_brewer(palette = "Set1", direction = -1) +
    scale_fill_brewer("",palette = "Set1", direction = -1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          strip.text = element_text(size=10)) +
    ggtitle(g)
  print(p1 + guides(fill=FALSE))
  ggsave(filename = paste0('output/plots/wt-only/paper/',g,'-raincloud_HSCs.pdf'), device = 'pdf',
         height = 4, width=5)

  # Differentiated
  ggplot(rc.df %>% filter(celltype %in% celltypes[c(7:9,12)]),
         aes(x=orig.ident,y=get(g), fill = orig.ident))+
    geom_flat_violin(position = position_nudge(x = .2, y = 0),
                     adjust=par.adj, bw=par.bw, scale = 'width', alpha=0.6) +
    geom_point(aes(color = orig.ident),
               position = position_jitter(width = .15, height = 0),
               size = 1, alpha = 0.6, show.legend = FALSE) +
    facet_wrap(~celltype, nrow = 1, strip.position = 'bottom',
               labeller = labeller(celltype = label_wrap_gen(width=15, multi_line = TRUE))) +
    ylab('Expression')+xlab('Cell type')+
    theme_cowplot(font_size = 12)+
    guides(color = FALSE, fill = FALSE)+
    scale_color_brewer(palette = "Set1", direction = -1) +
    scale_fill_brewer("",palette = "Set1", direction = -1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          strip.text = element_text(size=8)) +
    ggtitle(g)
  ggsave(filename = paste0('output/plots/wt-only/paper/',g,'-raincloud_Diff.pdf'), device = 'pdf',
         height = 4, width=5)
  
  # All cell types
  ggplot(rc.df %>% filter(celltype %in% celltypes[c(1:3,6:9,12)]), 
         aes(x=orig.ident,y=get(g), fill = orig.ident))+
    geom_flat_violin(position = position_nudge(x = .2, y = 0),
                     adjust=par.adj, bw=par.bw, scale = 'width', alpha=0.6) +
    geom_point(aes(color = orig.ident),
               position = position_jitter(width = .15, height = 0), 
               size = 1, alpha = 0.6, show.legend = FALSE) +
    facet_wrap(~celltype, nrow = 1, strip.position = 'bottom', 
               labeller = labeller(celltype = label_wrap_gen(width=15, multi_line = TRUE))) +
    ylab('Expression')+xlab('Cell type')+
    theme_cowplot(font_size = 16)+
    guides(color = FALSE, fill = FALSE)+
    scale_color_brewer(palette = "Set1", direction = -1) +
    scale_fill_brewer("",palette = "Set1", direction = -1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          strip.text = element_text(size=15)) +
    ggtitle(g)
  ggsave(filename = paste0('output/plots/wt-only/paper/',g,'-raincloud_full.pdf'), device = 'pdf',
         height = 4, width=10)
}

ggpubr::as_ggplot(ggpubr::get_legend(p1))
ggsave(filename = paste0('output/plots/wt-only/paper/legend.pdf'), device = 'pdf',
       height = 1, width=2)


# Feature plots -----------------------------------------------------------

pt.col = 'red'
           
# Plotted individually
for (g in plt.genes) {
  g.vals = FetchData(fc.lite.wt, vars = g)
  cs = scale_color_gradientn(colors = c('gray95',pals::stepped2()[20:17]), 
                             limits = c(min(g.vals), max(g.vals)))
  
  FeaturePlot(subset(fc.lite.wt, subset = group == 'wt.naive'),
                     features = g, cols = c('gray95','red'),
              keep.scale = 'feature',
              pt.size = 0.5, order = TRUE) + cs
  ggsave(filename = paste0('output/plots/wt-only/paper/',g,'-feature_wt-naive.pdf'), device = 'pdf', width = 4, height = 3.5)

  FeaturePlot(subset(fc.lite.wt, subset = group == 'wt.infected'),
              features = g, cols = c('gray95','red'),
              keep.scale = 'feature',
              pt.size = 0.5, order = TRUE) + cs
  ggsave(filename = paste0('output/plots/wt-only/paper/',g,'-feature_wt-infected.pdf'), device = 'pdf', width = 4, height = 3.5)
}

## UMAP plot with HSC/Bcells in red
umap.cols = c(RColorBrewer::brewer.pal(12,"Paired"),'grey')[c(1:3,6,5,4,7:13)]
DimPlot(fc.lite.wt, reduction = "umap", label = TRUE, group.by = 'celltype',
        cols = umap.cols, 
        label.size = 3) + 
  NoLegend()
ggsave(filename = paste0('output/plots/wt-only/paper/umap-celltype-labeled.pdf'), device = 'pdf',
       height = 5, width=5.5)

# Split by treatment
DimPlot(fc.lite.wt, reduction = "umap", label = TRUE, group.by = 'celltype',
        cols = umap.cols, split.by = 'group',
        label.size = 3) + 
  NoLegend()
ggsave(filename = paste0('output/plots/wt-only/paper/umap-celltype-labeled_by-group.pdf'), device = 'pdf',
       height = 5, width=9)


# Save legend
p1 = DimPlot(fc.lite.wt, reduction = "umap", label = TRUE, group.by = 'celltype',
             cols = umap.cols) 
ggpubr::as_ggplot(ggpubr::get_legend(p1))
ggsave(filename = paste0('output/plots/wt-only/paper/celltype-legend.pdf'), device = 'pdf',
       height = 3, width=3)

