### Setup env
library(dplyr)
library(tibble)
library(stringr) 
library(tidyr)

setwd("/mfd_abundance_tables")

## Import metadata file
fl16s.metadata <- data.table::fread("data/2024-01-25_OTU_minimal_metadata.csv", header = TRUE, sep = ",")

## Number of distinct samples
metadata.samples <- fl16s.metadata %>%
  select(fieldsample_barcode) %>%
  filter(str_detect(fieldsample_barcode, "MFD")) %>%
  distinct()

## Import OTU table
df.fl16s <- data.table::fread("data/2024-02-13_MFD_OTU_table.txt", sep = "\t", header = TRUE) %>%
  rename(OTU = 1)

## Extract sample IDs
mfd.samples <- colnames(df.fl16s) %>%
  str_subset(., "MFD")

## Filter sample IDs, remove empty rows
df.fl16s <- df.fl16s %>%
  select(OTU, all_of(mfd.samples)) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Import OTU taxonomy and format
fl16s.tax <- data.table::fread("data/2023-12-07_MFD_OTU_taxonomy.sintax", 
                               sep = "\t", header = FALSE, col.names = c("OTU", "Tax", "Strand", "Taxonomy")) %>%
  select(OTU, Taxonomy) %>%
  mutate(across(Taxonomy, ~str_replace(., 'tax=[a-zA-Z]:', 'd:')),
         across(Taxonomy, ~str_remove_all(., '[a-zA-Z]:'))) %>%
  separate(Taxonomy, into = paste0('tax', 1:7), sep = ',', extra = 'merge', fill = 'right') %>%
  rename(Kingdom = tax1,
         Phylum = tax2,
         Class = tax3,
         Order = tax4,
         Family = tax5,
         Genus = tax6,
         Species = tax7)

## Create long format
sample.df.long <- df.fl16s %>%
  pivot_longer(!OTU, names_to = "Sample", values_to = "Counts")

### Check total sequence counts from UMI and OPERON data respectively
## UMI
UMI <- sample.df.long %>%
  filter(str_detect(Sample, "Biowide_|MFDlight_")) %>% # Biowide and MFDlight are the UMI sequencing projects
  filter(!Counts == 0) %>%
  left_join(fl16s.tax)

UMI %>%
  pull(Counts) %>%
  sum()

## Unwanted taxa
OPERON %>%
  filter(Order == "Chloroplast" | Family == "Mitochondria" | Kingdom == "Eukaryota" | Kingdom == "Archaea") %>%
  pull(Counts) %>%
  sum()

## Without unwanted taxa
UMI %>%
  filter(!Kingdom == "Eukaryota",
         !Kingdom == "Archaea",
         !Order == "Chloroplast",
         !Family == "Mitochondria") %>%
  pull(Counts) %>%
  sum()

## Total UMI
107350+6200841

## Check
UMI %>%
  filter(!Kingdom == "Bacteria") %>%
  filter(!Counts == 0) %>%
  pull(Kingdom)

## OPERON
OPERON <- sample.df.long %>%
  filter(!str_detect(Sample, "Biowide_|MFDlight_")) %>% # Biowide and MFDlight are the UMI sequencing projects
  left_join(fl16s.tax)

OPERON %>%
  pull(Counts) %>%
  sum()

## Unwanted taxa
OPERON %>%
  filter(Order == "Chloroplast" | Family == "Mitochondria" | Kingdom == "Eukaryota" | Kingdom == "Archaea") %>%
  pull(Counts) %>%
  sum()

## Without unwanted taxa
OPERON %>%
  filter(!Kingdom == "Eukaryota",
         !Kingdom == "Archaea",
         !Order == "Chloroplast",
         !Family == "Mitochondria") %>%
  pull(Counts) %>%
  sum()

54475+14563108

OPERON %>%
  filter(!Kingdom == "Bacteria") %>%
  filter(!Counts == 0) %>%
  pull(Kingdom)

### Format data
## Select only MFD projects
sample.df.long.mfd <- sample.df.long %>%
  mutate(across(Sample, ~str_remove(., "Biowide_|MFDlight_|MFDextra.*_")))

## Distinct samples
distinct.df.mfd <- sample.df.long.mfd %>% select(Sample) %>% distinct()

## Group and summarise 
grouped.df.mfd <- sample.df.long.mfd %>%
  group_by(OTU, Sample) %>%
  summarise(Combined = sum(Counts))

## Wide format
sample.df.mfd <- grouped.df.mfd %>%
  group_by(OTU) %>%
  pivot_wider(names_from = "Sample", values_from = "Combined") %>%
  ungroup()

## Join with OTU taxonomy
otu.df.mfd <- sample.df.mfd %>%
  left_join(fl16s.tax, by = "OTU") %>%
  select(OTU, starts_with("MFD"), Kingdom:Species)

### Remove empty rows, remove unwanted taxa
otu.df_sub <- otu.df.mfd %>%
  filter(rowSums(across(where(is.integer))) != 0) %>%
  column_to_rownames(var = "OTU") %>%
  filter(Kingdom == "Bacteria") %>%
  mutate(across(Kingdom:Species, ~na_if(., ""))) %>%
  filter(!Kingdom == "Eukaryota",
         !Order == "Chloroplast",
         !Family == "Mitochondria")

## Write file to output directory
data.table::fwrite(otu.df_sub, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_OTU.csv"))

## Clean
rm(sample.df.long.mfd, distinct.df.mfd, grouped.df.mfd, sample.df.mfd, otu.df.mfd,
   df.fl16s, fl16s.metadata, mfd.samples)
gc()


### Aggregate to each taxonomic level and write to output directory

## Phylum
otu.df_sub.phylum <- otu.df_sub %>%
  mutate(across(Phylum, ~replace_na(., "Unclassified"))) %>%
  group_by(Phylum) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  select(-c(Phylum), Phylum) %>%
  mutate(Class = NA,
         Order = NA,
         Family = NA,
         Genus = NA,
         Species = NA) %>%
  left_join(fl16s.tax %>% select(Kingdom, Phylum) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.phylum, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_phylum.csv"))

rm(otu.df_sub.phylum)
gc()

## Class
otu.df_sub.class <- otu.df_sub %>%
  mutate(across(Class, ~replace_na(., "Unclassified"))) %>%
  group_by(Class) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Order = NA,
         Family = NA,
         Genus = NA,
         Species = NA) %>%
  left_join(fl16s.tax %>% select(Kingdom, Phylum, Class) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Phylum, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.class, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_class.csv"))

rm(otu.df_sub.class)
gc()

## Order
otu.df_sub.order <- otu.df_sub %>%
  mutate(across(Order, ~replace_na(., "Unclassified"))) %>%
  group_by(Order) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Family = NA,
         Genus = NA,
         Species = NA) %>%
  left_join(fl16s.tax %>% select(Kingdom, Phylum, Class, Order) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Class, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.order, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_order.csv"))

rm(otu.df_sub.order)
gc()

## Family
otu.df_sub.family <- otu.df_sub %>%
  mutate(across(Family, ~replace_na(., "Unclassified"))) %>%
  group_by(Family) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Genus = NA,
         Species = NA) %>%
  left_join(fl16s.tax %>% select(Kingdom, Phylum, Class, Order, Family) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Family, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.family, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_family.csv"))

rm(otu.df_sub.family)
gc()

## Genus
otu.df_sub.genus <- otu.df_sub %>%
  mutate(across(Genus, ~replace_na(., "Unclassified"))) %>%
  group_by(Genus) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Species = NA) %>%
  left_join(fl16s.tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Genus, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.genus, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_genus.csv"))

rm(otu.df_sub.genus)
gc()

## Species
otu.df_sub.species <- otu.df_sub %>%
  mutate(across(Species, ~replace_na(., "Unclassified"))) %>%
  group_by(Species) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  left_join(fl16s.tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Species, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.species, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_FL16S_aggregated_species.csv"))

rm(otu.df_sub.species)
gc()

rm(list=ls())
gc()

