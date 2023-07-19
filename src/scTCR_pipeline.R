library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(circlize)
library(rstatix)
source("src/clone_expansion_plots.R")
source('src/load_BD_scTCR.R')


proj_id <- c('P27352_1001', 'P27352_1002', 'P27470_1001', 'P27470_1002', 'P27566_1001', 'P28564_1001', 'P28564_1002')
sample_list = list(LongCovid1=c("LC1_415", "LC1_393", "LC1_898", "LC1_176", "LC1_407", "LC1_109"),
                   LongCovid2=c("LC2_168", "LC2_427", "LC2_342", "LC2_156", "LC2_104", "LC2_383"),
                   LongCovid3=c("LC3_937", "LC3_503", "LC3_443", "LC3_413", "LC3_358", "LC3_140"),
                   LongCovid4=c("LC4_891", "LC4_159", "LC4_323", "LC4_364", "LC4_528", "LC4_443", "LC4_941"),
                   LongCovid5=c("LC5_331", "LC5_395", "LC5_101", "LC5_944", "LC5_324"),
                   Convalescent1=c("CP76", "CP76b", "CP77", "CP95", "CP96", "CP97"),
                   Convalescent2=c("CP98", "CP163", "CP164b", "CP168", "CP172", "CP176"))
# load data
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
clone_exp_sub <- data %>%
  filter(cell_type == 'CD4T') %>%
  # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')
write_csv(clone_exp_sub, file = 'data/clone_expansion_CD4T.csv')
clone_exp_sub <- data %>%
  filter(cell_type == 'CD8T') %>%
  # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')
write_csv(clone_exp_sub, file = 'data/clone_expansion_CD8T.csv')



# donut chart 
for (patient in c('LongCovid1', 'LongCovid2', 'LongCovid3', 'LongCovid4', 'LongCovid5', 'Convalescent1', 'Convalescent2')){
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
  mutate(V = sub('\\*.*$', '', TCR_Beta_Delta_V_gene_Dominant))
for (patient in c('LongCovid1', 'LongCovid2', 'LongCovid3', 'LongCovid4', 'LongCovid5', 'Convalescent1', 'Convalescent2')){
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


data_va <- data %>% 
  mutate(V = sub('\\*.*$', '', TCR_Alpha_Gamma_V_gene_Dominant))
for (patient in c('LongCovid1', 'LongCovid2', 'LongCovid3', 'LongCovid4', 'LongCovid5', 'Convalescent1', 'Convalescent2')){
  for (sub_name in c('CD4T', 'CD8T', 'gdT', 'others')){
    fig_dir <- file.path('figures', patient)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    data_va_sub <- data_va %>% filter(cell_type == sub_name)
    exp_id <- data_va_sub %>%
      group_by(CDR3_concat, Sample_Name) %>%
      summarise(clone_count = n()) %>%
      ungroup() %>%
      inner_join(clone_id_map, by='CDR3_concat') %>%
      filter(clone_count > 2) %>% select(clone_id) %>% unlist()
    data_va_sub <- data_va_sub %>% mutate(clone_id = if_else(clone_id %in% exp_id, paste0('TCR', clone_id), 'unique'))
    
    pdf(file = file.path(fig_dir, paste0(patient, '_va_chord_diagram_', sub_name, '.pdf')),
        width = 20, height =10)
    par(mfrow = c(2, 3))
    lapply(sample_list[[patient]],function(x){clone_vb_chord_plot(data_va_sub, x)})
    dev.off()
  }
}

# public clone

public_id <- data %>% 
  filter(cell_type == 'CD8T', !Sample_Name %in% c('Hepatitis_Italy', 'LC3_443')) %>% 
  group_by(CDR3_concat, Sample_Name) %>% 
  summarise(clone_count=n()) %>%
  ungroup() %>%
  left_join(clone_id_map, by='CDR3_concat') %>%
  #filter(clone_count>1) %>%
  group_by(clone_id) %>%
  tally() %>% 
  filter(n>1)
# no public clone at all

# stats comparison
pacsr <- read_csv('data/PACSR_tmp.csv') %>% mutate(name=as.character(name))
sub_name <- 'CD8T'
clone_exp_sub <- data %>%
  filter(cell_type == sub_name) %>%
  # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')
tb <- lapply(unlist(sample_list),function(x){
  df <- get_clone_expansion_table(x, clone_exp_sub)
  return(df %>% filter(label=='Unique'))
}) %>% do.call(what = rbind)

tb <- tb %>% filter(!Sample_Name %in% c('Hepatitis_Italy', 'LC3_443')) %>% 
  mutate(name = gsub('LC.+_', '', Sample_Name)) %>% 
  left_join(pacsr, by='name') %>% 
  mutate(expanded_fraction=1-fraction)

write_csv(tb, 'data/CD8T expanded fraction of 29 LongCOVID patient.csv')

ggplot(tb, aes(x=pacsr, y=expanded_fraction)) + geom_boxplot() + geom_jitter() + ggtitle(paste0('p=',(tb %>% t_test(expanded_fraction ~ pacsr))[,'p']))
ggsave('figures/expansion_PACSR_gdT.pdf')

# connect to cortisol results
d_merge <- read_delim('data/LC_meta_SweBelg.csv', delim = ';') %>%
  left_join(read_excel('data/acth cortisol results.xlsx') %>% rename(sample_id.ACTHcortisol = sample), by='sample_id.ACTHcortisol') %>%
  left_join(tb %>% rename(subject_id = name), by='subject_id')

d_merge %>% filter(!is.na(expanded_fraction)) %>% View()

ggplot(d_merge %>% filter(!is.na(expanded_fraction)), aes(x=`age`, y=expanded_fraction, color=pacsr)) + 
  geom_point() 
+ geom_smooth(method='lm') +stat_cor(method = "spearman")

# clonal expansion ratio comparison

group_info <- tibble::enframe(sample_list) %>% 
  tidyr::unnest(cols = 'value') %>% 
  mutate(group = gsub('.{1}$', '', name)) %>%
  rename(Sample_Name = value)

sub_name <- 'CD8T'

clone_exp_sub <- data %>%
  filter(cell_type == sub_name) %>%
  # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')

d <- lapply(unique(clone_exp_sub$Sample_Name), function(x) {
  get_clone_expansion_table(x, clone_exp_sub) %>% dplyr::filter(label %in% c('Unique'))
}) %>%
  do.call(rbind, .) %>%
  dplyr::filter(Sample_Name != 'Hepatitis_Italy') %>% 
  group_by(Sample_Name) %>%
  summarise(fraction = sum(fraction)) %>%
  left_join(group_info, by='Sample_Name')

d.stats <- d %>% t_test(fraction ~ group, p.adjust.method = 'fdr') %>% add_xy_position(x = 'group')
ggplot(d, aes(x=group, y=fraction, color=group)) + geom_point() + 
  stat_pvalue_manual(d.stats, tip.length = 0.01, label='p') + 
  ylab('fraction of Unique') + ggtitle('CD8T')


# chord diagram vb21.3 (TRBV11-2) combine all
data <- read_csv(file = 'data/scTCR_data_merge.csv') %>% filter(Sample_Name != 'Hepatitis_Italy')
data_vb <- data %>% 
  mutate(V = sub('\\*.*$', '', TCR_Beta_Delta_V_gene_Dominant))
data_vb_sub <- data_vb %>% filter(cell_type == 'CD4T')

clone_id_map <- data %>% select(CDR3_concat, clone_id) %>% unique()
exp_id <- data_vb_sub %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat') %>%
  filter(clone_count > 1) %>% select(clone_id) %>% unlist()
data_vb_sub <- data_vb_sub %>% mutate(
  clone_id = if_else(clone_id %in% exp_id, paste0('TCR', clone_id), 'unique')) %>%
  mutate(group = case_when(
    grepl('LC', Sample_Name) ~ 'LongCOVID',
    grepl('CP', Sample_Name) ~ 'Convalescent',
    TRUE ~ 'others'
  ))

# separate LongCOVID/convalescent ------
pdf(file = 'figures/chordDiagram_TRBV_TCR_groups_CD8T.pdf',
    width = 20, height =10)
par(mfrow = c(1, 2))

df <- data_vb_sub %>% 
  filter(
    clone_id != 'unique',
    group == 'LongCOVID') %>%
  group_by(clone_id, V) %>%
  tally()
chord_order <- c(df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F),
                 df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F))
circos.clear()
chordDiagram(df, 
             annotationTrack = "grid", 
             order = chord_order,
             small.gap = 0,
             preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
title('LongCOVID')

df <- data_vb_sub %>% 
  filter(
    clone_id != 'unique',
    group == 'Convalescent') %>%
  group_by(clone_id, V) %>%
  tally()
chord_order <- c(df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F),
                 df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F))
circos.clear()
chordDiagram(df, 
             annotationTrack = "grid", 
             order = chord_order,
             small.gap = 0,
             preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
title('Convalescent')

dev.off()

# highlight TRBV11-2 ---------
library(RColorBrewer)
pdf(file = 'figures/chordDiagram_TRBV_TCR_groups_CD4T_highlightTRBV11-2.pdf',
    width = 20, height =10)
par(mfrow = c(1, 2))
df <- data_vb_sub %>% 
  filter(
    clone_id != 'unique',
    group == 'LongCOVID') %>%
  mutate(V = if_else(V=='TRBV11-2', V, 'Others')) %>%
  group_by(V, clone_id) %>%
  tally()

chord_order <- c(df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F), 
                 df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F)
                 )
circos.clear()
grid.col = colorRampPalette(brewer.pal(10, "Blues"))(length(chord_order))
names(grid.col) = chord_order
grid.col['Others'] = 'darkgrey'
grid.col['TRBV11-2'] = 'red'
chordDiagram(df, 
             grid.col = grid.col,
             annotationTrack = "grid", 
             order = chord_order,
             small.gap = 0,
             preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
title('LongCOVID')

df <- data_vb_sub %>% 
  filter(
    clone_id != 'unique',
    group == 'Convalescent') %>%
  mutate(V = if_else(V=='TRBV11-2', V, 'Others')) %>%
  group_by(V, clone_id) %>%
  tally()

chord_order <- c(df %>% group_by(V) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(V) %>% unlist(use.names = F), 
                 df %>% group_by(clone_id) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% select(clone_id) %>% unlist(use.names = F)
)
circos.clear()
grid.col = colorRampPalette(brewer.pal(10, "Blues"))(length(chord_order))
names(grid.col) = chord_order
grid.col['Others'] = 'darkgrey'
grid.col['TRBV11-2'] = 'red'
chordDiagram(df, 
             grid.col = grid.col,
             annotationTrack = "grid", 
             order = chord_order,
             small.gap = 0,
             preAllocateTracks = list(track.height = max(strwidth(unique(c(df[[1]],df[[2]]))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

title('Convalescent')
dev.off()
