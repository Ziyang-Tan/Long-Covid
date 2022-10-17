library(dplyr)
library(readr)
source("src/clone_expansion_plots.R")
library(ggplot2)
library(Seurat)
library(SeuratDisk)

sc_data <- read_csv(file = 'data/scTCR_data_merge.csv') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv')
obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
cell_type <- read_csv('data/cell_types.csv')
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type

# top clone phenotype
top_clone_id <- get_top_expansion_id(clone_exp, 10)
data <- sc_data %>% filter(clone_id %in% top_clone_id)
#obj <- ScaleData(obj)

obj_exp <- subset(obj, subset = unique_index %in% data$unique_index)
obj_exp <- ScaleData(obj_exp)

FeaturePlot(obj_exp, features = c('CD3E', 'CD4', 'CD8A', 'ITGAE'))

tmp <- obj_exp[[]] %>% left_join(data, by = 'unique_index')
obj_exp$clone_id <- tmp$clone_id
Idents(obj_exp) <- 'clone_id'
all_features <- FindAllMarkers(obj_exp)
top_features <- all_features %>% 
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC)
DimPlot(obj_exp, group.by = 'clone_id')

tmp <- GetAssayData(obj_exp, slot = 'scale.data')
d <- hclust(dist(tmp[unique(top_features$gene),], method = 'euclidean'), method = 'average')
DoHeatmap(obj_exp, group.by = 'clone_id', features = c(d$labels[d$order]))
ggsave('figures/top_clone_phenotype.pdf', width = 15, height = 10, units = 'in')

# search for CD8 T with high CD39 and CD11b and their clone expansions

top_clone_id <- get_top_expansion_id(clone_exp %>% filter(grepl('ISAC99', Sample_Name)), 20)
data <- sc_data %>% filter(clone_id %in% top_clone_id)
obj_exp <- subset(obj, subset = unique_index %in% data$unique_index)
FeaturePlot(obj_exp, features = c('CD3E', 'CD8A', 'ENTPD1', 'ITGAM'))
FeaturePlot(obj, features = c('CD3E', 'CD8A', 'ENTPD1', 'ITGAM'))
ggsave('figures/all_cells_CD39_CD11b.pdf')

tmp <- obj_exp[[]] %>% left_join(data, by = 'unique_index')
obj_exp$clone_id <- tmp$clone_id
Idents(obj_exp) <- 'clone_id'
DoHeatmap(obj_exp, slot='data', group.by = 'clone_id', features = c('CD3E', 'CD8A', 'ENTPD1', 'ITGAM')) + 
  scale_fill_gradientn(colors = c("black", "yellow"))
ggsave('figures/ISAC99_CD39_CD11b_top_clones.pdf', width = 15, height = 5, units = 'in')


################################################################################
# gene DE analysis, highly expanded vs other cells
library(EnhancedVolcano)
for (cur_sample in unique(sc_data$Sample_Name)){
  for (cur_subpop in c('CD4T', 'CD8T', 'gdT')){
    cur_subpop <- 'CD8T'
    cur_sample <- 'ISAC99_6'
    
    obj_sub <- subset(obj, subset = Sample_Name == cur_sample & cell_type == cur_subpop & TCR_Paired_Chains)
    high_expand_id <- (clone_exp %>% filter(clone_count >= 10, Sample_Name == cur_sample) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
    high_expand_sc <- sc_data %>% filter(Sample_Name == cur_sample, cell_type == cur_subpop, clone_id %in% high_expand_id) %>% 
      tibble::add_column(highly_expand = 'highly_expanded') %>%
      select(unique_index, highly_expand)
    
    if (dim(high_expand_sc)[1] < 3) next # no highly expanded clones in this sample x cell type
    
    tmp <- obj_sub[[]] %>% left_join(high_expand_sc, by = 'unique_index') %>%
      tidyr::replace_na(list(highly_expand='others'))
    tmp2 <- tmp$highly_expand
    names(tmp2) <- tmp$unique_index
    obj_sub$highly_expand <- tmp2
    Idents(obj_sub) <- 'highly_expand'
    res <- FindMarkers(obj_sub, ident.1 = 'highly_expanded', ident.2 = 'others',
                       logfc.threshold = 0)
    EnhancedVolcano(res, x='avg_log2FC', y='p_val', lab = rownames(res), 
                    pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
                    selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
                    title = paste0(cur_subpop, ' of ', cur_sample),
                    subtitle = 'highly_expanded cells vs others')
    ggsave(paste0('figures/gene_DE/', cur_sample, '_', cur_subpop, '.pdf'))
    
  }
}




