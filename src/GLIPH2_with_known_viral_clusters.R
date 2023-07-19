library(dplyr)
library(readr)
library(igraph)
source('src/GLIPH_related_plots.R')
source('src/clone_expansion_plots.R')

subpop <- 'CD8'

data_TRA_raw <- read_csv(paste0('data/', subpop, '_TRA_cluster_annotated.csv')) %>% filter(sid != 'Hepatitis_Italy') %>% mutate(index = as.character(index))
data_TRB_raw <- read_csv(paste0('data/', subpop, '_TRB_cluster_annotated.csv')) %>% filter(sid != 'Hepatitis_Italy') %>% mutate(index = as.character(index))

data_TRAB <- rbind(data_TRA_raw, data_TRB_raw) %>%
  mutate(CDR3aa_concat = paste0(cdr3a, '_', cdr3b),
         Sample_Name = sid,
         cluster_id = paste0(index, '_', stringr::str_extract(pattern, '^.{4}')),
         virus_specific = !is.na(annotation))

clone_exp <- read_csv(file = 'data/clone_expansion.csv') %>%
  filter(Sample_Name != 'Hepatitis_Italy') %>%
  mutate(clone_id = as.character(clone_id),
         patient_id = sub('_.*$', '', Sample_Name))
data_scTCR <- read_csv('data/scTCR_data_merge.csv') %>% filter(Sample_Name != 'Hepatitis_Italy')
# top_clones <- read_csv('data/top_clones_all_timepoints.csv') %>%
#   filter(sub_name == 'CD8T')
nt2aa <- data_scTCR %>% select(CDR3_concat, CDR3aa_concat) %>% distinct()

# write virus specificity 
virus_specific_list <- left_join(
  data_scTCR %>% select(CDR3aa_concat, clone_id) %>% distinct(),
  data_TRAB %>% select(CDR3aa_concat, virus_specific) %>% distinct(), by = 'CDR3aa_concat') %>%
  filter(virus_specific) %>%
  select(clone_id) %>% 
  distinct()

write_csv(virus_specific_list, paste0('data/virus_specific_clone_id_', subpop, '.csv'))

# SARS CoV2 specificity (only TRB)
data_lc <- read_csv('/Users/tan/Long-Covid/data/GLIPH_results/longCovid_df.csv')
data_conv <- read_csv('/Users/tan/Long-Covid/data/GLIPH_results/convCovid_df.csv')

tmp <- rbind(data_lc %>% filter(condition == 'LongCOVID') %>% mutate(CoV2 = 'SARS-CoV2'),
             data_conv %>% filter(condition == 'Convalescent') %>% mutate(CoV2 = 'SARS-CoV2')) %>%    
  mutate(cdr3b = sub('^C', '', cdr3b))
# merge two GLIPH results
# data_TRB <- data_TRB %>% mutate(index = paste0('1_', index)) %>% select(index, sid, cdr3a, cdr3b, annotation) %>%
#   full_join(tmp %>% mutate(index = paste0('2_',index)) %>%select(index, sid, cdr3a, cdr3b, CoV2), by=c('cdr3a', 'cdr3b')) %>% distinct()

clone_types <- rbind(data_TRB_raw %>% select(cdr3a, cdr3b),
                     tmp %>% select(cdr3a, cdr3b)) %>% distinct()

data_TRB <- clone_types %>% 
  left_join(data_TRB_raw %>% mutate(index = paste0('1_', index)) %>% select(index, sid, cdr3a, cdr3b, annotation), by=c('cdr3a', 'cdr3b')) %>%
  left_join(tmp %>% mutate(index = paste0('2_',index)) %>% select(index, sid, cdr3a, cdr3b, CoV2), by=c('cdr3a', 'cdr3b')) %>%
  mutate(index = if_else(is.na(index.x), index.y, index.x),
         sid = if_else(is.na(sid.x), sid.y, sid.x)) %>%
  select(cdr3a, cdr3b, index, sid, annotation, CoV2) %>% distinct()


data_TRA <- data_TRA_raw %>% select(index, sid, cdr3a, cdr3b, annotation) %>% distinct()

# parse annotation
data_TRA <- data_TRA %>% 
  mutate(CMV = if_else(grepl('CMV', annotation), 'CMV', ''),
         EBV = if_else(grepl('EBV', annotation), 'EBV', '') ,
         HomoSapiens = if_else(grepl('HomoSapiens', annotation), 'HomoSapiens', ''),
         InfluenzaA = if_else(grepl('InfluenzaA', annotation), 'InfluenzaA', ''),
         YFV = if_else(grepl('YFV', annotation), 'YFV', ''),
         tuberculosis = if_else(grepl('tuberculosis', annotation), 'tuberculosis', ''),
         HTLV = if_else(grepl('HTLV-1', annotation), 'HTLV-1', ''),
         HIV = if_else(grepl('HIV-1', annotation), 'HIV-1', ''),
         Influenza = if_else(grepl('Influenza', annotation), 'Influenza', ''),
         DENV = if_else(grepl('DENV', annotation), 'DENV', '' ),
         CoV2 = ''
  ) %>%
  mutate(# CDR3a = paste(CMV, EBV, HomoSapiens, InfluenzaA, YFV, tuberculosis, HTLV, HIV, Influenza, DENV, sep=':'), 
    CDR3a = paste(CMV, EBV, InfluenzaA, Influenza, HomoSapiens,CoV2, sep=':'), 
    CDR3a = gsub('(:)+',':',CDR3a),
    CDR3a = gsub('(^:)|(:$)', '', CDR3a),
    CDR3b = '')
data_TRB <- data_TRB %>% 
  mutate(CMV = if_else(grepl('CMV', annotation), 'CMV', ''),
         EBV = if_else(grepl('EBV', annotation), 'EBV', '') ,
         HomoSapiens = if_else(grepl('HomoSapiens', annotation), 'HomoSapiens', ''),
         InfluenzaA = if_else(grepl('InfluenzaA', annotation), 'InfluenzaA', ''),
         YFV = if_else(grepl('YFV', annotation), 'YFV', ''),
         tuberculosis = if_else(grepl('tuberculosis', annotation), 'tuberculosis', ''),
         HTLV = if_else(grepl('HTLV-1', annotation), 'HTLV-1', ''),
         HIV = if_else(grepl('HIV-1', annotation), 'HIV-1', ''),
         Influenza = if_else(grepl('Influenza', annotation), 'Influenza', ''),
         DENV = if_else(grepl('DENV', annotation), 'DENV', '' ),
         CoV2 = if_else(is.na(CoV2), '', CoV2)
  ) %>%
  mutate(# CDR3b = paste(CMV, EBV, HomoSapiens, InfluenzaA, YFV, tuberculosis, HTLV, HIV, Influenza, DENV, sep=':'), 
    CDR3b = paste(CMV, EBV, InfluenzaA, Influenza, HomoSapiens, CoV2, sep=':'), 
    CDR3b = gsub('(:)+',':',CDR3b),
    CDR3b = gsub('(^:)|(:$)', '', CDR3b),
    CDR3a = '')

# network of all samples
pdf(paste0('figures/GLIPH_network/GLIPH_with_CoV2_clusters_clonecount', 3, '_', subpop, '.pdf'))
plot.gliph(data_TRA, data_TRB, clone_exp, nt2aa, clone_thre=3)
title('all samples')
dev.off()

tp <- data_TRAB %>%
  filter(cluster_id %in% c('6_TRAV', '15_TRAV', '16_TRAV', '17_TRAV', '276_TRBV', '287_TRBV', '365_TRBV')) %>%
  select(CDR3aa_concat) %>%
  distinct() %>%
  left_join(nt2aa) %>%
  left_join(clone_exp)

write_csv(tp, paste0('data/chosen_meta_cluster_', subpop, 'T.csv'))

# pdf(paste0('figures/GLIPH_with_known_viral_clusters_clonecount', 3, '_with_label_', subpop, '.pdf'))
# plot.gliph(data_TRA, data_TRB, clone_exp, nt2aa, clone_thre=3, labels = tp$clone_id)
# title('all Long covid samples')
# dev.off()

# network per sample
for (name_i in unique(clone_exp$Sample_Name)){
  pdf(paste0('figures/GLIPH_network/expanded clones/GLIPH_with_CoV2_clusters_clonecount', 1, '_', name_i, '_', subpop, '.pdf'))
  plot.gliph(data_TRA %>% filter(sid==name_i), 
             data_TRB %>% filter(sid==name_i), 
             clone_exp %>% filter(Sample_Name==name_i), 
             nt2aa, 
             clone_thre=1,
             #labels = tp$clone_id
             #labels = filter(top_clones, patient==name_i)$clone_id
  )
  title(main=name_i)
  dev.off()
}

data_TRB %>% filter(!grepl('LC', sid))

# clonal expansion and changes on GLIPH cluster level

# cluster_exp <- data_TRAB %>% 
#   group_by(pattern, Sample_Name, cluster_id) %>% 
#   summarise(clone_count=sum(frequency)) %>%
#   rename(clone_id = cluster_id,
#          CDR3_concat = pattern)
#mutate(patient = sub('_.*$', '', Sample_Name))

# g_list <- lapply(unique(data_TRAB$sid), clone_expansion_alluvium, cluster_exp, top_mod = 'large', n=10)
# ggarrange(plotlist = g_list, ncol = 2, nrow = 2) %>%
#   ggexport(filename = 'figures/GLIPH_cluster_changes_CD8T.pdf')


# analysis of graph
# comp <- components(g)
# clusters <- groups(comp)
# big_clusters <- clusters[comp$csize > 3]
# 
# 
# tmp <- data_scTCR %>% filter(clone_id %in% clusters$`5`)
# tmp %>% group_by(Sample_Name, clone_id) %>% tally() %>% View()
# 
# data_TRAB %>% filter(CDR3aa_concat %in% unique(tmp$CDR3aa_concat)) %>% View()


