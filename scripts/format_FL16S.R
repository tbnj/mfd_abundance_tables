### Setup env
library(tidyverse)

setwd("/mfd_abundance_tables")

## Import OTU table
operon <- data.table::fread("data/MFD_operon_fSSU_vs_MFG_ssu_database_v1.1_NR987_MFD_only_OTUtab.txt", sep = "\t", header = TRUE) %>%
  rename(OTU = 1) %>%
  rename_with(., ~str_remove(., "_pool"),
              starts_with("MFD")) %>%
  # select(-MFD10339) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

UMI <- data.table::fread("data/MFD_UMI_fSSU_vs_MFG_ssu_database_v1.1_NR987_MFD_only_OTUtab.txt", sep = "\t", header = TRUE) %>%
  rename(OTU = 1) %>%
  rename_with(., ~str_remove(., "_pool"),
              starts_with("MFD")) %>%
  # select(-MFD10339) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

samples <- c(colnames(UMI), colnames(operon)) %>%
  str_subset(., "MFD") %>%
  unique() %>%
  as.data.frame()

OTUs <- c(UMI$OTU, operon$OTU) %>%
  unique()

## Import OTU taxonomy and format
taxonomy <- data.table::fread("data/MFG_ssu_database_v1.2_NR987_taxonomy.txt",
                               sep = "\t", header = FALSE, col.names = c("OTU")) %>%
  separate(OTU, into = c('OTU', 'Taxonomy'), sep = ';tax=', extra = 'merge', fill = 'right') %>%
  mutate(across(OTU, ~str_remove_all(., '>')),
         across(Taxonomy, ~str_remove_all(., '[a-zA-Z]:|;'))) %>%
  filter(OTU %in% OTUs) %>%
  separate(Taxonomy, into = paste0('tax', 1:7), sep = ',', extra = 'merge', fill = 'right') %>%
  rename(Kingdom = tax1,
         Phylum = tax2,
         Class = tax3,
         Order = tax4,
         Family = tax5,
         Genus = tax6,
         Species = tax7) %>%
  mutate(across(Phylum, ~str_replace(., "MFD_p_179760|MFD_p_312341", "Patescibacteria")))

## Combine tables
df.UMI <- UMI %>%
  left_join(taxonomy, by = "OTU") %>%
  select(OTU, starts_with("MFD"), Kingdom:Species) %>%
  filter(rowSums(across(where(is.integer))) != 0) %>%
  column_to_rownames(var = "OTU") %>%
  filter(Kingdom == "Bacteria") %>%
  mutate(across(Kingdom:Species, ~na_if(., ""))) %>%
  filter(!Kingdom == "Eukaryota",
         !Order == "Chloroplast",
         !Family == "Mitochondria")

df.operon <- operon %>%
  left_join(taxonomy, by = "OTU") %>%
  select(OTU, starts_with("MFD"), Kingdom:Species) %>%
  filter(rowSums(across(where(is.integer))) != 0) %>%
  column_to_rownames(var = "OTU") %>%
  filter(Kingdom == "Bacteria") %>%
  mutate(across(Kingdom:Species, ~na_if(., ""))) %>%
  filter(!Kingdom == "Eukaryota",
         !Order == "Chloroplast",
         !Family == "Mitochondria")


## Write file to output directory
data.table::fwrite(df.UMI, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_16S_UMI_OTU.csv"))

data.table::fwrite(df.operon, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_16S_operon_OTU.csv"))
