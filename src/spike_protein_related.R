library(dplyr)
library(readr)
library(Seurat)
library(SeuratDisk)
library(scatterpie)
library(tidyr)
library(wesanderson)
source("src/clone_expansion_plots.R")

meta <- read_csv2('/Users/tan/Long-Covid/data/LC_meta_SweBelgRoma_ext.csv') %>%
  filter(cohort %in% c('Long covid_solna (UppCov)', 'Pan')) %>% 
  rename(Sample_Name = subject_id)
clone_exp <- read_csv(file = 'data/clone_expansion_CD8T.csv') %>% 
  filter(Sample_Name != 'Hepatitis_Italy') %>%
  mutate(clone_id = as.character(clone_id),
         Sample_Name = sub('(LC\\d_)', '', Sample_Name))
COVID_specific_clones <- read_csv('data/COVID_specific_clones_CD8T.csv')

## CoV2 specific expansion

tmp <- lapply(unique(clone_exp$Sample_Name), function(name){
  tb <- get_clone_expansion_table(name, clone_exp) %>%
    mutate(label = case_when(
      clone_id %in% COVID_specific_clones$clone_id ~ label,
      TRUE ~ 'Others'
    )) %>%
    group_by(Sample_Name, label) %>% summarise(value = sum(fraction)) %>%
    ungroup()
}) %>% do.call(rbind, .)

# shannon diversity
# diversity <- clone_exp %>% group_by(Sample_Name) %>% summarise(shannon = vegan::diversity(clone_count),
#                                                   lnS = log(n())) %>%
#   mutate(shannon_equitability = shannon/lnS)
# tmp <- clone_exp %>% filter(Sample_Name == '101') %>% 
#   mutate(p = clone_count/sum(clone_count),
#          lnp = log(p),
#          plnp = p*lnp)


data <- left_join(tmp, meta %>% select(Sample_Name, days_since_diseaseCoV1, `simoa_conc.`, cohort))

tmp <- data %>% pivot_wider(names_from = label, values_from = value, values_fill = 0) %>%
  filter(!is.na(`simoa_conc.`)) %>%
  mutate(spike_conc = as.numeric(`simoa_conc.`), 
         log2_spike_conc = log2(spike_conc + 1),
         #cohort = jitter(10*as.numeric(as.factor(cohort)), amount = 3),
         CoV2_specific_expansion_fraction = 1-Others
         )

# regression
ggplot(tmp, aes(x=CoV2_specific_expansion_fraction, y=log2_spike_conc)) + 
  geom_point() +
  geom_smooth(method = 'lm', se=FALSE) +
  ggpubr::stat_cor(method = "pearson") + 
  theme_bw() + 
  theme(panel.background = element_blank())
ggsave('figures/clonal_expansion_spike_CD8M_CoV2specific_linear.pdf', width = 8, height = 4.5)

# non-linear regression
glm_model <- glm(log2_spike_conc ~ CoV2_specific_expansion_fraction, family = poisson(link = "log"), data = tmp)
summary(glm_model)
glm_model$model$fitted <- predict(glm_model, type = "response")
ggplot(glm_model$model) + 
  geom_point(aes(CoV2_specific_expansion_fraction, log2_spike_conc)) +
  geom_line(aes(CoV2_specific_expansion_fraction, fitted)) +
  theme_bw() + 
  theme(panel.background = element_blank()) +
  annotate('text', label = 'p of intercept < 2e-16 ***, 
           p of CoV2_specific_expansion_fraction = 0.00163 **', x=0.15, y=9)
ggsave('figures/clonal_expansion_spike_CD8M_CoV2specific_glm.pdf', width = 8, height = 4.5)


# scatterpie
# cols = c('#FFFFFF', wes_palette("Rushmore1", 7, type = "continuous"))
# ggplot() + 
#   geom_scatterpie(aes(x=CoV2_specific_expansion_fraction*100, y=log2_spike_conc, r=0.5), data = tmp,
#                   cols=c('Others', '2', '3', '4', '5', '6-9', '10-19', '>=20')) + 
#   coord_fixed() + 
#   scale_fill_manual(values = cols) +
#   theme_bw() + 
#   theme(panel.background = element_blank())
# ggsave('figures/clonal_expansion_spike_CD8M_CoV2specific.pdf')

## not CoV2 specific expansion

tmp2 <- lapply(unique(clone_exp$Sample_Name), function(name){
  tb <- get_clone_expansion_table(name, clone_exp) %>%
    mutate(label = case_when(
      clone_id %in% COVID_specific_clones$clone_id ~ 'CoV2 specific', 
      label == 'Unique' ~ 'Unique',
      TRUE ~ 'non CoV2 specific'
    )) %>%
    group_by(Sample_Name, label) %>% summarise(value = sum(fraction)) %>%
    ungroup()
}) %>% do.call(rbind, .)

data <- left_join(tmp2, meta %>% select(Sample_Name, days_since_diseaseCoV1, `simoa_conc.`, cohort))

tmp2 <- data %>% pivot_wider(names_from = label, values_from = value, values_fill = 0) %>%
  filter(!is.na(`simoa_conc.`)) %>%
  mutate(spike_conc = as.numeric(`simoa_conc.`), 
         log2_spike_conc = log2(spike_conc + 1),
         #cohort = jitter(10*as.numeric(as.factor(cohort)), amount = 3),
  )

ggplot(tmp2, aes(x=`non CoV2 specific`, y=log2_spike_conc)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = 0) + 
  ggpubr::stat_cor(method = "pearson") +
  theme_bw() + 
  theme(panel.background = element_blank())
ggsave('figures/clonal_expansion_spike_CD8M_nonCoV2specific_linear.pdf', width = 8, height = 4.5)

write_csv(tmp2, 'data/CoV2_specific_clonal_expansion_fraction_GLIPH2.csv')

## all expansion

# tmp <- lapply(unique(clone_exp$Sample_Name), function(name){
#   tb <- get_clone_expansion_table(name, clone_exp) %>%
#     group_by(Sample_Name, label) %>% summarise(value = sum(fraction)) %>%
#     ungroup()
# }) %>% do.call(rbind, .)
# 
# data <- left_join(tmp, meta %>% select(Sample_Name, days_since_diseaseCoV1, `simoa_conc.`, cohort))
# 
# tmp <- data %>% pivot_wider(names_from = label, values_from = value, values_fill = 0) %>%
#   filter(!is.na(`simoa_conc.`)) %>%
#   mutate(log2_spike_conc = log2(as.numeric(`simoa_conc.`)),
#          #cohort = jitter(10*as.numeric(as.factor(cohort)), amount = 3),
#          clonal_expansion_fraction = 1-Unique
#   )
# 
# ggplot(tmp, aes(x=clonal_expansion_fraction, y=log2_spike_conc)) + geom_point()
# type 1
# ggplot() + 
#   geom_scatterpie(aes(x=days_since_diseaseCoV1, y=log2_spike_conc, r=0.5), data = tmp,
#                   cols=c('Unique', '2', '3', '4', '5', '6-9', '10-19', '>=20')) + coord_fixed() + 
#   scale_x_continuous(
#     trans=scales::trans_new(name = '10log2_trans', transform = function(x){10*log2(x)}, inverse = function(x){2^(x/10)})
#     )
# type 2
# cols = c('#FFFFFF', wes_palette("Rushmore1", 7, type = "continuous"))
# ggplot() + 
#   geom_scatterpie(aes(x=clonal_expansion_fraction*100, y=log2_spike_conc, r=0.5), data = tmp,
#                   cols=c('Unique', '2', '3', '4', '5', '6-9', '10-19', '>=20')) + 
#   coord_fixed() + 
#   scale_fill_manual(values = cols) +
#   theme_bw() + 
#   theme(panel.background = element_blank())
# 
# ggsave('figures/clonal_expansion_spike_CD8M_all_expansion.pdf')





