library(dplyr)
library(readr)
source("src/clone_expansion_plots.R")
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(EnhancedVolcano)
library(rstatix)

sc_data <- read_csv(file = 'data/scTCR_data_merge.csv') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv')
obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
cell_type <- read_csv('data/cell_types.csv')
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type
virus_specific <- read_csv('data/virus_specific_clone_id_CD8.csv')
cov2_specific <- read_csv('data/COVID_specific_clones_CD8T.csv')

# top clone phenotype
# top_clone_id <- get_top_expansion_id(clone_exp, 10)
# data <- sc_data %>% filter(clone_id %in% top_clone_id)
#obj <- ScaleData(obj)

# obj_exp <- subset(obj, subset = unique_index %in% data$unique_index)
# obj_exp <- ScaleData(obj_exp)
# 
# FeaturePlot(obj_exp, features = c('CD3E', 'CD4', 'CD8A', 'ITGAE'))
# 
# tmp <- obj_exp[[]] %>% left_join(data, by = 'unique_index')
# obj_exp$clone_id <- tmp$clone_id
# Idents(obj_exp) <- 'clone_id'
# all_features <- FindAllMarkers(obj_exp)
# top_features <- all_features %>% 
#   group_by(cluster) %>% 
#   slice_max(n=10, order_by = avg_log2FC)
# DimPlot(obj_exp, group.by = 'clone_id')
# 
# tmp <- GetAssayData(obj_exp, slot = 'scale.data')
# d <- hclust(dist(tmp[unique(top_features$gene),], method = 'euclidean'), method = 'average')
# DoHeatmap(obj_exp, group.by = 'clone_id', features = c(d$labels[d$order]))
# ggsave('figures/top_clone_phenotype.pdf', width = 15, height = 10, units = 'in')

################################################################################
# gene DE analysis, highly expanded vs other cells
# in all samples
cur_subpop = 'CD8T'
obj_sub <- subset(obj, subset = cell_type == cur_subpop & TCR_Paired_Chains)
obj_sub <- subset(obj, subset = cell_type == cur_subpop & 
                    TCR_Paired_Chains & 
                    unique_index %in% (sc_data %>% filter(!(clone_id %in% virus_specific)))$unique_index) # those are all non virus specific cells
high_expand_id <- (clone_exp %>% filter(clone_count >= 20) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
# expansion threshold for normal DE is usually 6
high_expand_sc <- sc_data %>% filter(cell_type == cur_subpop, clone_id %in% high_expand_id, !(clone_id %in% virus_specific)) %>% 
  tibble::add_column(highly_expand = 'highly_expanded') %>%
  select(unique_index, highly_expand)
tmp <- obj_sub[[]] %>% left_join(high_expand_sc, by = 'unique_index') %>%
  tidyr::replace_na(list(highly_expand='others'))
tmp2 <- tmp$highly_expand
names(tmp2) <- tmp$unique_index
obj_sub$highly_expand <- tmp2
Idents(obj_sub) <- 'highly_expand'
res <- FindMarkers(obj_sub, ident.1 = 'highly_expanded', ident.2 = 'others',
                   logfc.threshold = 0)
EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-04, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = res %>% filter(abs(avg_log2FC) > 1 | p_val < 1e-04) %>% rownames(),
                title = paste0(cur_subpop, ' of all samples'),
                subtitle = 'highly_expanded(>20) cells vs others')
ggsave('figures/gene_DE/highly_expanded(>20) cells vs others_virus excluded.pdf', width = 14, height = 14)

# DE cov2 specific
cur_subpop = 'CD8T'
expand_id <- (clone_exp %>% filter(clone_count >= 2) %>% select(clone_id) %>% unique())[[1]]
obj_sub <- subset(obj, subset = cell_type == cur_subpop & 
                    TCR_Paired_Chains & 
                    unique_index %in% (sc_data %>% filter(!(clone_id %in% virus_specific) & clone_id %in% expand_id))$unique_index,) # those are all non virus specific cells
cov2_sc <- sc_data %>% filter(cell_type == cur_subpop, clone_id %in% cov2_specific$clone_id) %>% 
  tibble::add_column(cov2 = 'SARS-CoV-2 specific') %>%
  select(unique_index, cov2)
tmp <- obj_sub[[]] %>% left_join(cov2_sc, by = 'unique_index') %>%
  tidyr::replace_na(list(cov2='Other CD8+ T cells'))
tmp2 <- tmp$cov2
names(tmp2) <- tmp$unique_index
obj_sub$cov2 <- tmp2
Idents(obj_sub) <- 'cov2'
res <- FindMarkers(obj_sub, ident.1 = 'SARS-CoV-2 specific', ident.2 = 'Other CD8+ T cells',
                   logfc.threshold = 0)
EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-04, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = res %>% filter(abs(avg_log2FC) > 1 | p_val < 1e-04) %>% rownames(),
                title = paste0(cur_subpop, ' of all samples'),
                xlim=c(-2,2),
                col= 'black',
                arrowheads=F,
                subtitle = 'highly_expanded(>20) cells vs others') + 
  theme_minimal()
ggsave('figures/gene_DE/CoV2 expansion vs other expansion.pdf', width = 7, height = 7)







# specific cluster vs the others expanded (>1)
tp <- read_csv(paste0('data/chosen_meta_cluster_', cur_subpop, '.csv'))
expand_id <- (clone_exp %>% filter(clone_count >= 2) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
obj_sub <- subset(obj, subset = cell_type == cur_subpop & 
                    TCR_Paired_Chains & 
                    unique_index %in% (sc_data %>% filter(clone_id %in% expand_id))$unique_index) # those are all the expanded cells
chosen_sc <- sc_data %>% filter(cell_type == cur_subpop, clone_id %in% tp$clone_id) %>% 
  tibble::add_column(chosen = 'chosen_cluster') %>%
  select(unique_index, chosen)
tmp <- obj_sub[[]] %>% left_join(chosen_sc, by = 'unique_index') %>%
  tidyr::replace_na(list(chosen='others'))
tmp2 <- tmp$chosen
names(tmp2) <- tmp$unique_index
obj_sub$chosen <- tmp2
Idents(obj_sub) <- 'chosen'
res <- FindMarkers(obj_sub, ident.1 = 'chosen_cluster', ident.2 = 'others',
                   logfc.threshold = 0)
EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-04, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = res %>% filter(abs(avg_log2FC) > 1 | p_val < 1e-04) %>% rownames(),
                title = paste0(cur_subpop, ' of all samples'),
                subtitle = 'cells in chosen_cluster vs others')
ggsave('figures/gene_DE/cells in chosen GLIPH cluster vs others CD8T.pdf', width = 14, height = 14)

VlnPlot(object = obj_sub, features = c('RORC', 'KLRB1', 'KLRG1', 'LGALS1', 'IFNG', 'ITGAM'))
ggsave('figures/gene_DE/cells in chosen GLIPH cluster vs others CD8T_some features.pdf')



# per sample
# for (cur_sample in unique(sc_data$Sample_Name)){
#   for (cur_subpop in c('CD4T', 'CD8T', 'gdT')){
#     #cur_subpop <- 'CD8T'
#     #cur_sample <- 'ISAC99_6'
#     obj_sub <- subset(obj, subset = Sample_Name == cur_sample & cell_type == cur_subpop & TCR_Paired_Chains)
#     high_expand_id <- (clone_exp %>% filter(clone_count >= 10, Sample_Name == cur_sample) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
#     high_expand_sc <- sc_data %>% filter(Sample_Name == cur_sample, cell_type == cur_subpop, clone_id %in% high_expand_id) %>% 
#       tibble::add_column(highly_expand = 'highly_expanded') %>%
#       select(unique_index, highly_expand)
#     
#     if (dim(high_expand_sc)[1] < 3) next # no highly expanded clones in this sample x cell type
#     
#     tmp <- obj_sub[[]] %>% left_join(high_expand_sc, by = 'unique_index') %>%
#       tidyr::replace_na(list(highly_expand='others'))
#     tmp2 <- tmp$highly_expand
#     names(tmp2) <- tmp$unique_index
#     obj_sub$highly_expand <- tmp2
#     Idents(obj_sub) <- 'highly_expand'
#     res <- FindMarkers(obj_sub, ident.1 = 'highly_expanded', ident.2 = 'others',
#                        logfc.threshold = 0)
#     EnhancedVolcano(res, x='avg_log2FC', y='p_val', lab = rownames(res), 
#                     pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
#                     selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
#                     title = paste0(cur_subpop, ' of ', cur_sample),
#                     subtitle = 'highly_expanded cells vs others')
#     ggsave(paste0('figures/gene_DE/', cur_sample, '_', cur_subpop, '.pdf'))
#     
#   }
# }

# the abundance of the chosen cluster in each sample
cur_subpop = 'CD8T'
tp <- read_csv(paste0('data/chosen_meta_cluster_', cur_subpop, '.csv'))
pacsr <- read_csv('data/PACSR_tmp.csv') %>% mutate(name=as.character(name))
sc_data <- read_csv(file = 'data/scTCR_data_merge.csv') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
sample_cell_count <- sc_data %>%
  filter(cell_type == cur_subpop, !Sample_Name %in% c('Hepatitis_Italy', 'LC3_443', NA)) %>%
  group_by(Sample_Name) %>%
  summarise(cell_count = n())
tb <- sc_data %>%
  filter(cell_type == cur_subpop) %>%
  # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
  group_by(CDR3_concat, Sample_Name, clone_id) %>%
  summarise(clone_count = n()) %>%
  filter(clone_count>1) %>%
  filter(!Sample_Name %in% c('Hepatitis_Italy', 'LC3_443', NA),
         clone_id %in% tp$clone_id) %>%
  group_by(Sample_Name) %>%
  summarise(count = sum(clone_count)) %>% 
  mutate(name = gsub('LC.+_', '', Sample_Name)) %>%
  left_join(pacsr, by='name') %>%
  left_join(sample_cell_count, by='Sample_Name') %>%
  mutate(ratio = count/cell_count)

ggplot(tb, aes(x=pacsr, y=ratio)) + geom_boxplot() + geom_jitter() + ggtitle(paste0('p=',(tb %>% t_test(ratio ~ pacsr))[,'p']))
ggsave('figures/raito of chosen cluster_PACSR_CD8T.pdf')


# some chosen gene expression
# KIR
FeaturePlot(obj_sub, c('KIR2DL4', 'KIR2DL5B', 'KIR3DL1'))
ggsave(filename = 'figures/KIRpos/all CD8T.pdf')
# DimPlot(obj_sub, group.by='highly_expand')

obj_sub_highexp <- subset(obj_sub, subset = highly_expand == 'highly_expanded')
FeaturePlot(obj_sub_highexp, c('KIR2DL4', 'KIR2DL5B', 'KIR3DL1'))
ggsave(filename = 'figures/KIRpos/highly expanded(>20) CD8T.pdf')



