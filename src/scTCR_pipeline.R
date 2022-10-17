library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(circlize)

proj_id <- c('P27352_1001', 'P27352_1002')
sample_list = list(LongCovid1=c("LC1_415", "LC1_393", "LC1_898", "LC1_176", "LC1_407", "LC1_109"),
                   LongCovid2=c("LC2_168", "LC2_427", "LC2_342", "LC2_156", "LC2_104", "LC2_383"))
# load data
source('src/load_BD_scTCR.R')
glob_path <- '/Users/tan/OneDrive - KI.SE/TCR_processed_data/single cell'
raw_tcr <- lapply(proj_id, BD_load_VDJ, dir_path = glob_path) %>% do.call(what = rbind)
sample_tag <- lapply(proj_id, BD_load_sample_tag, dir_path = glob_path) %>% do.call(what = rbind)

cell_type <- read_csv('data/cell_types.csv')

raw_tcr_merge <- left_join(raw_tcr, sample_tag, by = 'unique_index') %>%
  left_join(cell_type, by = 'unique_index') %>%
  mutate(Sample_Name = case_when(
    is.na(Sample_Name) ~ proj_id,
    TRUE ~ Sample_Name
  ))

# summarize 
df <- raw_tcr_merge %>% filter(!Sample_Name %in% c('Multiplet', 'Undetermined'))
table_summary <- df %>% group_by(proj_id) %>% summarise(demultiplexed = n()) %>% left_join(
  df %>% filter(at_least_one_chain) %>% group_by(proj_id) %>% summarise(at_least_one_CDR3 = n())
) %>% left_join(
  df %>% filter(TCR_Paired_Chains) %>% group_by(proj_id) %>% summarise(paired = n())
)

data <- raw_tcr_merge %>%
  filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
  filter(!is.na(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant)) %>%
  filter(!Sample_Name %in% c('Multiplet', 'Undetermined')) %>%
  mutate(CDR3_concat = paste0(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, '_', 
                              TCR_Beta_Delta_CDR3_Nucleotide_Dominant),
         CDR3aa_concat = paste0(TCR_Alpha_Gamma_CDR3_Translation_Dominant, '_', 
                                TCR_Beta_Delta_CDR3_Translation_Dominant)) %>%
  mutate(clone_id = as.character(as.numeric(as.factor(CDR3_concat))))
write_csv(data, file = 'data/scTCR_data_merge.csv')
clone_id_map <- data %>% select(CDR3_concat, clone_id) %>% unique()
clone_exp <- data %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')
#inner_join(data %>% select(clone_id, CDR3aa_concat) %>% unique(), by='clone_id')
write_csv(clone_exp, file = 'data/clone_expansion.csv')

# donut chart 
source("src/clone_expansion_plots.R")
for (patient in c('LongCovid1', 'LongCovid2')){
  for (sub_name in c('CD4T', 'CD8T', 'gdT', 'others')){
    fig_dir <- file.path('figures', patient)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    clone_exp_sub <- data %>%
      filter(cell_type == sub_name) %>%
      # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
      group_by(CDR3_concat, Sample_Name) %>%
      summarise(clone_count = n()) %>%
      ungroup() %>%
      inner_join(clone_id_map, by='CDR3_concat')
    g_list1 <- lapply(sample_list[[patient]],function(x){clone_expansion_donut(x, clone_exp_sub)})
    ggarrange(plotlist = g_list1, ncol = 3, nrow = 2) %>%
      ggexport(filename = file.path(fig_dir, paste0(patient, '_clone_expansion_', sub_name, '.pdf')), 
               width = 20, height = 10)
  }
}

# chord diagram vb
source('src/chord_plots.R')
data_vb <- data %>% 
  mutate(Vb = sub('\\*.*$', '', TCR_Beta_Delta_V_gene_Dominant))
for (patient in c('LongCovid1', 'LongCovid2')){
  for (sub_name in c('CD4T', 'CD8T', 'gdT', 'others')){
    fig_dir <- file.path('figures', patient)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    data_vb_sub <- data_vb %>% filter(cell_type == sub_name)
    exp_id <- data_vb_sub %>%
      group_by(CDR3_concat, Sample_Name) %>%
      summarise(clone_count = n()) %>%
      ungroup() %>%
      inner_join(clone_id_map, by='CDR3_concat') %>%
      filter(clone_count > 2) %>% select(clone_id) %>% unlist()
    data_vb_sub <- data_vb_sub %>% mutate(clone_id = if_else(clone_id %in% exp_id, paste0('TCR', clone_id), 'unique'))
    
    pdf(file = file.path(fig_dir, paste0(patient, '_vb_chord_diagram_', sub_name, '.pdf')),
        width = 20, height =10)
    par(mfrow = c(2, 3))
    lapply(sample_list[[patient]],function(x){clone_vb_chord_plot(data_vb_sub, x)})
    dev.off()
  }
}













