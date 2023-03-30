library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)

n.var.features = 4000
cr.data.dir = "data/cellranger/"

# Set up cell cycle data --------------------------------------------------

# Cell cycle markers, from Tirosh et al, 2015
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Convert from human to murine genes
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- c(convertHumanGeneList(s.genes),'Atad2','Pold3') # Manually convert missing genes
m.g2m.genes <- c(convertHumanGeneList(g2m.genes),'Hn1','Fam64a') # Synonyms: Hn1 = Jpt1, Fam64a = Pimreg
m.g2m.genes <- m.g2m.genes[-grep('Pimreg',m.g2m.genes)]
m.cc.genes <- c(m.s.genes, m.g2m.genes)

# Save murine cell cycle gene lists
save(m.s.genes,m.g2m.genes,m.cc.genes, file='data/m.cc.genes.Rdata')


# Create data objects -----------------------------------------------------

data.dirs <- c('wt.naive' = '17422',
               'wt.infected' = '17423',
               'ko.naive' = '17424',
               'ko.infected' = '17425')
exp.groups <- names(data.dirs)

# Iterate through experimental groups to create objects
for (grp in exp.groups) {

  # Load raw data
  grp.data <- Read10X(data.dir = paste0(cr.data.dir,
                                        data.dirs[grp],
                                        "/outs/filtered_feature_bc_matrix/"))

  # Initialize the Seurat object with the raw data
    # Features in at least 10 cells
    # Cells with at least 200 features
  grp.data <- CreateSeuratObject(counts = grp.data,
                                 project = grp,
                                 min.cells = 10, min.features = 200)

  # Add some QC stats for filtering cells
  grp.data[["percent.mt"]] <- PercentageFeatureSet(grp.data, pattern = "^mt-")
  grp.data[["log10GenesPerUMI"]] <- log10(grp.data$nFeature_RNA) / log10(grp.data$nCount_RNA)

  saveRDS(grp.data, file = paste0('data/',grp,'-base.rds'))

  # Subset data and visually verify filtering
  grp.data <- subset(grp.data, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 10 & nCount_RNA < 60000)

  # Normalize data
  grp.data <- NormalizeData(grp.data, normalization.method = "LogNormalize", scale.factor = 10000)

  # Find highly variable features
  grp.data <- FindVariableFeatures(grp.data, selection.method = "disp", nfeatures = n.var.features)

  # Save filtered data
  saveRDS(grp.data, file = paste0('data/',grp,'-filtered-unscaled.rds'))

  # Calculate cell cycle scores and phase assignments
  grp.data <- CellCycleScoring(grp.data, s.features = m.s.genes, g2m.features = m.g2m.genes)

  # Scale by cell cycle scores
  grp.data <- ScaleData(grp.data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(grp.data))

  # Find highly variable features
  grp.data <- FindVariableFeatures(grp.data, selection.method = "disp", nfeatures = n.var.features)

  # Save scaled data
  saveRDS(grp.data, file = paste0('data/',grp,'-filtered-cc-scaled.rds'))
}

rm(grp.data)

# Integrate pre-scaled data -----------------------------------------------

# Load WT naive
wt.naive = readRDS('data/wt.naive-filtered-cc-scaled.rds')
wt.infected = readRDS('data/wt.infected-filtered-cc-scaled.rds')
ko.naive = readRDS('data/ko.naive-filtered-cc-scaled.rds')
ko.infected = readRDS('data/ko.infected-filtered-cc-scaled.rds')

# Combine into list object
full.list <- list('wt.naive' = wt.naive, 'wt.infected' = wt.infected, 'ko.naive' = ko.naive, 'ko.infected' = ko.infected)
rm(wt.naive,wt.infected,ko.naive,ko.infected) # remove original large objects

# Get variable features
full.list <- lapply(X = full.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "disp", nfeatures = n.var.features)
})

# Perform integration
full.anchors <- FindIntegrationAnchors(object.list = full.list, anchor.features = n.var.features, dims = 1:50)

full.combined <- IntegrateData(anchorset = full.anchors,  dims = 1:50)
save(full.combined, file = 'data/full.combined-pre-cc-scaled.Rdata')

