library(dplyr)
library(readr)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(EnhancedVolcano)
library(rstatix)

data_lc <- read_csv('/Users/tan/Long-Covid/data/GLIPH_results/longCovid_df.csv')
data_conv <- read_csv('/Users/tan/Long-Covid/data/GLIPH_results/convCovid_df.csv')

sc_data <- read_csv(file = 'data/scTCR_data_merge.csv') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
cell_type <- read_csv('data/cell_types.csv')
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type
clone_exp <- read_csv(file = 'data/clone_expansion.csv')

get_cdr3aa_set <- function(d, group) {
  d %>% 
    filter(condition == group,
           !is.na(cdr3a)) %>% 
    mutate(cdr3b_trim = sub('^C', '', cdr3b),
           CDR3aa_concat = paste0(cdr3a, '_', cdr3b_trim)) %>%
    select(CDR3aa_concat) %>%
    distinct() %>%
    unlist()
}

COVID_specific_clones <- sc_data %>% 
  filter(cell_type == 'CD8T',
         CDR3aa_concat %in% union(get_cdr3aa_set(data_lc, 'LongCOVID'), get_cdr3aa_set(data_conv, 'Convalescent'))
  )

write_csv(COVID_specific_clones, 'data/COVID_specific_clones_CD8T.csv')


COVID_specific <- rbind(
  sc_data %>% filter(CDR3aa_concat %in% get_cdr3aa_set(data_lc, 'LongCOVID')) %>% select(unique_index) %>% mutate(COVID_specific = 'LongCOVID'),
  sc_data %>% filter(CDR3aa_concat %in% get_cdr3aa_set(data_conv, 'Convalescent')) %>% select(unique_index) %>% mutate(COVID_specific = 'Convalescent')
)

rownames(COVID_specific) <- COVID_specific$unique_index

data_conv %>% group_by(condition) %>% tally
