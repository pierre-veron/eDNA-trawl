# Contact 
__Pierre Veron__ 

📧 pierre [dot] veron [dot] 2017 [at] polytechnique [dot] org

# Project 
This project is the source and code associated with the publication:

> Pierre Veron, Romane Rozanski, Virginie Marques, Stéphane Joost, Marie Emilie Deschez, Verena M Trenkel, Pascal Lorance, Alice Valentini, Andrea Polanco F., Loïc Pellissier, David Eme, Camille Albouy, **Environmental DNA complements scientific trawling in surveys of marine fish biodiversity**, _ICES Journal of Marine Science_, 2023. DOI: [10.1093/icesjms/fsad139](https://doi.org/10.1093/icesjms/fsad139)

In this project, we compare two sampling methods for fishes (eDNA and scientific trawling) under the scope of multi-component indices of diversity. The data are taken from the EVHOE 2019 measurement campaign in the Bay of Biscay (Northern Atlantic). We compared eDNA metabarcoding and trawling for assessing:
* taxonomic diversity
* functional diversity
* phylogenetic diversity
* abundance. 
 
This repository allows to reproduce the analyses and the figures in this publication, and the same data. 
The whole analysis is perfored in `R` version `4.1.2`.

# How to use it
The script `main.R` calls the used R packages, loads the data, runs the analyses and draws the figures. This script calls all the other scripts in the right order. 


# Project architecture

```
├───data
│   └───icons
│       └───char
├───diversity_indices
├───figures
├───output
│   ├───trees
│   ├───trees_merged
│   ├───trees_null
│   └───trees_null_merged
├───R
│   ├───functions
│   │   ├───fonctions_source_f
│   │   └───fonctions_source_p
│   ├───graphics
│   └───process_data
└───renv

```
* `main.R` the main script that calls the other scripts
* `renv.lock` information on the environments and versions
* `eDNA_trawling.Rproj` metadata for the project (for RStudio)
* `renv/`:directory loaded by `renv` for restoring the environment (see above)
* `data/` the data associated with the project:
   * `eDNA_2019_reads.csv` the number of reads from the metabarcoding analysis, per station and taxa,
   * `eDNA_metadata.csv` information about the conditions of sampling for eDNA (position, depth...)
   * `GISdata.rds` background of the maps used in the figures
   * `incidence_matrix.csv` presence or absence of each taxon in each station and each sampling method
   * `incidence_matrix_eDNA.csv` incidence matrix only for eDNA
   * `incidence_matrix_trawl.csv` incidence matrix only for trawling
   * `name_stations.csv` name of the stations
   * `reftax.csv` the classification of each taxa detected by the methods
   * `reftax_species.csv` list of all Atlantic species belonging to the genus detected by the methods and their classification
   * `traits.csv` functional traits used to generate the functional space
   * `traits_cat.csv` list of the types of traits
   * `trawl_2019_abundance.csv` the number of individuals caught by the trawl per taxa and site
   * `trawl_metadata.csv` position of each trawling sample
   * `trees.rds` phylogenetic trees used for the phylogenetic indices of diversity. 
* `diversity_indices/` the formatted output of the scripts:
   * `alpha.csv`, `alpha_merged.csv`, `alpha_ses.csv`, `alpha_merged_ses.csv`: alpha indices of diversity on each of the sampling methods and on the merged communities (trawling+eDNA) and their SES values
   * `phylo_signal_alpha` and `phylo_signal_gamma.csv`: result from the analysis of phylogenetic signal of clustering on the trees. 
* `figures/` contains the figures generated by the scripts and used in the manuscript
* `output/` contains the raw output from the script (in RDS format). If the file contained in this directory are cleaned, the program will re-run the analyses (it might take a long time). Please note that this directory should contain at least the folders `trees`, `trees_merged`, `trees_null` and `trees_null_merged`. Once the analyses are run, their results are stored here and accessed by the script `R/process_data/diversity_indices.R` to avoid computing several times the same indices. 
* `R/` contains the scripts used to load the data, run the analyses, store the results and draw the figures. 
