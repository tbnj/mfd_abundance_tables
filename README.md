# mfd_abundance_tables
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 
The scripts are used to aggregate, transform the data into a format compatible with [ampvis2](https://kasperskytte.github.io/ampvis2/articles/ampvis2.html), and finally extract observational tables used in the analysis included in the manuscript. 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### Amplicon 16S data 
`scripts/aggregate_FL16s.R` aggregates data from 16S sequences based on the Nanopore-UMI and PacBio operon sequencing. 

### Metagenomic 16S data 
`scripts/aggregate_arcbac.R` aggreate data from 16S fragments dervied from the metagenomes.


`scripts/create_ampvis_arcbac.R` imports the genus-aggregated table, filters and performs a single subsampling without replacement. 


`scripts/format_data_arcbac.R` formats observational tables in both count and relative abundance formats. 

### Metagenomic 18S data 
`scripts/aggregate_euk.R` aggreate data from 18S fragments dervied from the metagenomes. 


`scripts/create_ampvis_euk.R` imports the genus-aggregated table, filters and performs a single subsampling without replacement. 

## Data
The scripts rely on data files available from the MFD Zenodo [repo](https://zenodo.org/records/12605769) and the MFD [github](https://github.com/cmc-aau/mfd_metadata), from where the original output files are also available. 
