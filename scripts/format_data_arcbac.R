### Setup env
library(tidyverse)
library(ampvis2)
library(vegan)

setwd("/mfd_abudance_tables")

### Load ampvis-formatted data
load("release/2025-02-13_MFD_ampvis_arcbac_data.RData")


### Filtered data
## Extract dataframe
arcbac.data <- mfd.ampvis.arcbac$abund %>%
  cbind(mfd.ampvis.arcbac$tax) %>%
  select(starts_with("MFD"), Kingdom:Genus)

## Change to relative abundance
arcbac.data.rel <- arcbac.data %>%
  mutate(across(where(is.numeric), ~./sum(.)*100))


### Rarefied data
## Extract dataframe
arcbac.data.ra <- mfd.ampvis.arcbac.ra$abund %>%
  cbind(mfd.ampvis.arcbac.ra$tax) %>%
  select(starts_with("MFD"), Kingdom:Genus)

## Change to relative abundance
arcbac.data.ra.rel <- arcbac.data.ra %>%
  mutate(across(where(is.numeric), ~./sum(.)*100))

### Write outputs to disk
data.table::fwrite(arcbac.data, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_arcbac_genus_count.csv"))
data.table::fwrite(arcbac.data.rel, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_arcbac_genus_rel.csv"))

data.table::fwrite(arcbac.data.ra, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_arcbac_genus_rarefaction_count.csv"))
data.table::fwrite(arcbac.data.ra.rel, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_arcbac_genus_rarefaction_rel.csv"))

rm(list=ls())
gc()
