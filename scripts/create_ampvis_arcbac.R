### Setup env
library(ampvis2)
library(tidyverse)

setwd("/mfd_abudance_tables")

### Import sample metadata
combined.metadata <- data.table::fread('/output/2024-03-07_combined_metadata.csv',
                                       sep = ",", na.strings = "") %>%
  relocate(fieldsample_barcode, .before = project_id)


### Import aggregated observational table
genus.table.arcbac <- data.table::fread('/output/2024-03-07_MFD_arcbac_aggregated_genus.csv', 
                                        sep = ",", header = TRUE)


### Filter data
## Find the 16S rRNA gene count for each sample
counts.arcbac <- genus.table.arcbac %>%
  select(where(is.numeric)) %>%
  colSums() %>%
  sort()

## Create dataframe
df.counts.arcbac <- counts.arcbac %>%
  as.data.frame() %>%
  rename(count = 1) %>%
  mutate(rank = seq(1, length(counts.arcbac), 1))

## Filter sample IDs based on 16S counts
filter.arcbac <- counts.arcbac[counts.arcbac >= 1000]

## Extract sample IDs
samples.filt <- names(filter.arcbac)

## Filter observational table and remove empty rows
genus.table.arcbac.filt <- genus.table.arcbac %>%
  select(any_of(samples.filt), Kingdom:Species) %>%
  filter(rowSums(across(where(is.numeric)))!=0)


### Create ampvis object
## Select table
data.arcbac.filt <- genus.table.arcbac.filt %>%
  select(where(is.numeric))

## Select taxonomy
tax.arcbac.filt <- genus.table.arcbac.filt %>%
  select(Kingdom:Species)

## Filter metadata 
metadata.filt <- combined.metadata %>%
  filter(fieldsample_barcode %in% samples.filt)

## Load filtered data into an ampvis object
mfd.ampvis.arcbac <- amp_load(
  otutable = genus.table.arcbac.filt,
  metadata = metadata.filt)


### Perform a single random subsampling without replacement
## Extract rarefy level as the smallest number of reads in the filtered dataset 
rarefy <- min(colSums(mfd.ampvis.arcbac[["abund"]]))

## Use built in rarefy function from ampvis2
set.seed(123) # needed as this is a random process
mfd.ampvis.arcbac.ra <- mfd.ampvis.arcbac %>%
  amp_subset_samples(., rarefy = rarefy, normalise = FALSE)

## Clean and save image to disk
rm(list=setdiff(ls(), c("metadata.filt", "mfd.ampvis.arcbac", "mfd.ampvis.arcbac.ra")))

save.image('release/2024-03-07_MFD-ampvis-arcbac-data.RData')
