### Setup env
library(ampvis2)
library(tidyverse)

setwd("/mfd_abudance_tables")

### Import sample metadata
combined.metadata <- data.table::fread('/output/2024-03-07_combined_metadata.csv',
                                       sep = ",", na.strings = "") %>%
  relocate(fieldsample_barcode, .before = project_id)


### Import aggregated observational table
genus.table.euk <- data.table::fread('/output/2024-03-07_MFD_euk_aggregated_genus.csv', 
                                        sep = ",", header = TRUE)


### Filter data
## Find the 16S rRNA gene count for each sample
counts.euk <- genus.table.euk %>%
  select(where(is.numeric)) %>%
  colSums() %>%
  sort()

## Create dataframe
df.counts.euk <- counts.euk %>%
  as.data.frame() %>%
  rename(count = 1) %>%
  mutate(rank = seq(1, length(counts.euk), 1))

## Filter sample IDs based on 16S counts
filter.euk <- counts.euk[counts.euk >= 50]

## Extract sample IDs
samples.filt <- names(filter.euk)

## Filter observational table and remove empty rows
genus.table.euk.filt <- genus.table.euk %>%
  select(any_of(samples.filt), Kingdom:Species) %>%
  filter(rowSums(across(where(is.numeric)))!=0)


### Create ampvis object
## Select table
data.euk.filt <- genus.table.euk.filt %>%
  select(where(is.numeric))

## Select taxonomy
tax.euk.filt <- genus.table.euk.filt %>%
  select(Kingdom:Species)

## Filter metadata 
metadata.filt <- combined.metadata %>%
  filter(fieldsample_barcode %in% samples.filt)

## Load filtered data into an ampvis object
mfd.ampvis.euk <- amp_load(
  otutable = genus.table.euk.filt,
  metadata = metadata.filt)


### Perform a single random subsampling without replacement
## Extract rarefy level as the smallest number of reads in the filtered dataset 
rarefy <- min(colSums(mfd.ampvis.euk[["abund"]]))

## Use built in rarefy function from ampvis2
set.seed(123) # needed as this is a random process
mfd.ampvis.euk.ra <- mfd.ampvis.euk %>%
  amp_subset_samples(., rarefy = rarefy, normalise = FALSE)

## Clean and save image to disk
rm(list=setdiff(ls(), c("metadata.filt", "mfd.ampvis.euk", "mfd.ampvis.euk.ra")))

save.image('output/2024-03-07_mfd-ampvis-euk-data.RData')
