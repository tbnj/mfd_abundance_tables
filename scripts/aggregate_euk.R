### Setup env
library(dplyr)
library(tibble)
library(stringr) 
library(tidyr)

setwd("/mfd_abundance_tables")

### Import metadata files
seq.metadata <- readr::read_csv("data/2023-10-11_samples_minimal_metadata_collapsed.csv")

sample.metadata <- readxl::read_excel("data/2024-02-13_mfd_db.xlsx") %>%
  filter(!project_id %in% c("P04_1", "P12_4", "P20_1"))

## Combine metadata and filter for unused projects
comb.metadata <- sample.metadata %>%
  left_join(seq.metadata, by = "fieldsample_barcode") %>%
  filter(!is.na(before_total_reads)) %>%
  mutate(across(mfd_sampletype:mfd_hab3, ~str_to_title(.)))

## Write combined metadata file to output directory
## Uncomment if arcbac-version of script was not run beforehand
#data.table::fwrite(comb.metadata, sep = ",", row.names = FALSE, col.names = TRUE,
#                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_combined_metadata.csv"))

## Sanity check
intersect(seq.metadata$fieldsample_barcode, sample.metadata$fieldsample_barcode) %>% length()
intersect(sample.metadata$fieldsample_barcode, seq.metadata$fieldsample_barcode) %>% length()

## Extract sample IDs
samples <- comb.metadata %>%
  pull(fieldsample_barcode)

## Import full observational table
otu.df <- data.table::fread("data/2024-03-07_MFD_euk_shallow_release.csv")

### Filter sample IDs, remove empty rows, remove unwanted taxa
##!! IMPORTANT column OTU is only named so, as it is needed by ampvis2 						 !!##
##!! Each entry does not represent a classical OTU, but can be thought of as a taxonomic bin !!##
otu.df_sub <- otu.df %>%
  select(OTU, any_of(samples), Kingdom:Species) %>%
  filter(rowSums(across(where(is.integer))) != 0) %>%
  column_to_rownames(var = "OTU") %>%
  mutate(across(Kingdom:Species, ~na_if(., "")))

## Select taxonomy
tax <- otu.df_sub %>% 
  select(Kingdom:Species)

rm(comb.metadata, otu.df, sample.metadata, seq.metadata)
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
  left_join(tax %>% select(Kingdom, Phylum) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.phylum, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_euk_aggregated_phylum.csv"))

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
  left_join(tax %>% select(Kingdom, Phylum, Class) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Phylum, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.class, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_euk_aggregated_class.csv"))

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
  left_join(tax %>% select(Kingdom, Phylum, Class, Order) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Class, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.order, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_euk_aggregated_order.csv"))

rm(otu.df_sub.order)
gc()

## Family
otu.df_sub.family <- otu.df_sub %>%
  mutate(across(Family, ~replace_na(., "Unclassified"))) %>%
  group_by(Family) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Genus = NA,
         Species = NA) %>%
  left_join(tax %>% select(Kingdom, Phylum, Class, Order, Family) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Family, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.family, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_euk_aggregated_family.csv"))

rm(otu.df_sub.family)
gc()

## Genus
otu.df_sub.genus <- otu.df_sub %>%
  mutate(across(Genus, ~replace_na(., "Unclassified"))) %>%
  group_by(Genus) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Species = NA) %>%
  left_join(tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus) %>% distinct()) %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Kingdom:Genus, ~replace_na(., "Unclassified")))

data.table::fwrite(otu.df_sub.genus, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_euk_aggregated_genus.csv"))

rm(otu.df_sub.genus)
gc()

## Species - a little different than the other taxonomic levels
## Filter for entries unclassified at the Species level, and overwrite as "Unclassified"
spec.acc <- otu.df_sub %>%
  filter(is.na(Species)) %>%
  mutate(across(Species, ~replace_na(., "Unclassified"))) %>%
  group_by(Species) %>%
  summarise(across(where(is.integer), ~sum(.))) %>%
  mutate(Kingdom = "Unclassified",
         Phylum = "Unclassified",
         Class = "Unclassified",
         Order = "Unclassified",
         Family = "Unclassified",
         Genus = "Unclassified") %>%
  select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species), Kingdom, Phylum, Class, Order, Family, Genus, Species)

## Remove entries not classified at the Species level, and add the overwritten version back
otu.df_sub.species <- otu.df_sub %>%
  filter(!is.na(Species)) %>%
  rbind(spec.acc)

data.table::fwrite(otu.df_sub.species, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_euk_aggregated_species.csv"))

rm(otu.df_sub.species)
gc()

rm(list=ls())
gc()
