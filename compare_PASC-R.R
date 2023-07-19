library(dplyr)
library(readr)
library(readxl)

meta1 <- read_excel('/Users/tan/Long-Covid/data/UppCov_Plasma_received_metadata 220714.xlsx')
meta2 <- read_excel('/Users/tan/Long-Covid/data/Brodinu_Sthlm_Patients_metadata_221118_PBedit.xlsx', sheet = 1)
tcr_ID <- read_delim('/Users/tan/Long-Covid/data/LC_sampleID.csv', delim = ';')

r1 <- meta1 %>% filter(`PACS /PACS recovered (PACS-R)` == 'PACS-R') %>% select(`Subject/Patient ID`)
r2 <- meta2 %>% filter(recoverd_PASC_at_sampling == 'yes') %>% select(`Record ID`)

intersect(unlist(r1), unlist(r2))

tmp <- meta1 %>% filter(`Subject/Patient ID` %in% tcr_ID$x)


intersect(unlist(r1), tcr_ID$x)
intersect(unlist(r2), tcr_ID$x)

pacsr_id <- union(intersect(unlist(r1), tcr_ID$x), intersect(unlist(r2), tcr_ID$x))

tibble(name=unique(tcr_ID$x),
       pacsr=if_else(name %in% pacsr_id, 'PACS-R', 'PACS')) %>%
  write_csv('data/PACSR_tmp.csv')
