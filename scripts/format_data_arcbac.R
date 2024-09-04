### Setup env
library(tidyverse)
library(ampvis2)
library(vegan)

setwd("/mfd_abudance_tables")

### Load ampvis-formatted data
load("output/2024-03-07_MFD-ampvis-arcbac-data.RData")


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

# Write outputs to disk
data.table::fwrite(arcbac.data, "output/2024-03-07_arcbac-count.csv")
data.table::fwrite(arcbac.data.rel, "output/2024-03-07_arcbac-rel.csv")

data.table::fwrite(arcbac.data.ra, "output/2024-03-07_arcbac-rarefaction-count.csv")
data.table::fwrite(arcbac.data.ra.rel, "output/2024-03-07_arcbac-rarefaction-rel.csv")

