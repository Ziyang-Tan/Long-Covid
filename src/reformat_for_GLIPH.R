library(dplyr)
library(readr)

data <- read_csv(file = 'data/scTCR_data_merge.csv')
group_info <- read_delim(file= 'data/group_info.csv', delim = ';')

clone_id_map <- data %>% select(CDR3_concat, clone_id) %>% unique()

sub_name = 'CD8T'
data_sub <- data %>% filter(cell_type == sub_name)
clone_exp_sub <- data %>%
  filter(cell_type == sub_name) %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')

df <- left_join(
  clone_exp_sub %>% select(clone_id, clone_count, Sample_Name),
  data_sub %>% select(clone_id, TCR_Alpha_Gamma_V_gene_Dominant, TCR_Alpha_Gamma_J_gene_Dominant, TCR_Alpha_Gamma_CDR3_Translation_Dominant, 
                      TCR_Beta_Delta_V_gene_Dominant, TCR_Beta_Delta_J_gene_Dominant, TCR_Beta_Delta_CDR3_Translation_Dominant) %>% distinct(clone_id, .keep_all=TRUE),
  by='clone_id'
) %>% filter(!is.na(TCR_Beta_Delta_CDR3_Translation_Dominant)) %>%
  left_join(group_info, by='Sample_Name') %>%
  filter(Group != 'PACS-R' & !is.na(Group)) %>%
  rename(sid = Sample_Name, 
         condition = Group,
         CDR3b = TCR_Beta_Delta_CDR3_Translation_Dominant,
         Vb = TCR_Beta_Delta_V_gene_Dominant,
         Jb = TCR_Beta_Delta_J_gene_Dominant,
         CDR3a = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
         Va = TCR_Alpha_Gamma_V_gene_Dominant,
         Ja = TCR_Alpha_Gamma_J_gene_Dominant,
         Frequency = clone_count) %>%
  select(sid, condition, CDR3b, Vb, Jb, CDR3a, Va, Ja, Frequency)

write_csv(df, file = 'data/data_LC/GLIPH_input_CD8T.csv')
