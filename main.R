### ----------------------------------------------------------------------------
#   This is the main script associated with the publication
#   Veron et al., eDNA complements scientific trawling to survey multiple 
#      components of marine fish biodiversity
#   
#   It calls the different scripts to run the analyses and produce the 
#      figures used in this publication. 
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
### ----------------------------------------------------------------------------

# Initialisation ---------------------------------------------------------------
rm(list = ls()) # automatically remove all variables
source("R/load_packages.R") 
source("R/functions/useful_functions.R", encoding = "utf8") # import useful functions 
source("R/load_data.R")



# Analysis ---------------------------------------------------------------------

source("R/functions/PIPELINE_TAXO_V1_26-03-2022.R") # load the pipelines
source("R/functions/PIPELINE_fonctio_V7_23-03-22.R")
source("R/functions/PIPELINE_phylo_V5_24-03-2022.R")
source("R/process_data/diversity_indices.R") # run the pipelines 
source("R/process_data/indices_export.R") # export diversity indices to a csv file
source("R/process_data/cluster_phylo.R") # phylogenetic signal


# Graphics ---------------------------------------------------------------------

# Define colors for the graphs and functions called to draw the plots
source("R/graphics/colors.R")
source("R/functions/functions_plot_pcoa.R")
source("R/functions/phylogenetic_functions.R")

# Draw figures (each script is independent)
source("R/graphics/1_gamma_div.R", encoding = "utf8")
source("R/graphics/2_species_accumulation_curve.R", encoding = "utf8")
source("R/graphics/3_alpha_div.R", encoding = "utf8")
source("R/graphics/4_alpha_div_merged.R", encoding = "utf8")
source("R/graphics/5_alpha_selected.R")
source("R/graphics/6_S6_PCoA.R", encoding = "utf8")

source("R/graphics/S1_order_by_method.R", encoding = "utf8")
source("R/graphics/S2_traits_by_method.R", encoding = "utf8")
source("R/graphics/S3_community_matrix.R", encoding = "utf8")
source("R/graphics/S4_Venn_diagrams.R", encoding = "utf8")
source("R/graphics/S5_functional_space.R", encoding = "utf8")
source("R/graphics/S6_S7_S8_abundance.R", encoding = "utf8")
