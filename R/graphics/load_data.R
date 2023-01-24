### ----------------------------------------------------------------------------
#   Loads the data. 
#    
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
#   Created 2023/01
### ----------------------------------------------------------------------------

Incidence_matrix <- read.csv("data/incidence_matrix.csv", row.names = 1)

Incidence_matrix_merged <- as_tibble(Incidence_matrix, rownames = "site")
Incidence_matrix_merged$site <- gsub("t", "", gsub("e", "", Incidence_matrix_merged$site))

Incidence_matrix_merged <- Incidence_matrix_merged %>% 
  group_by(site) %>%
  summarise_each(list(max)) # we take maximum : if species detected in at least
# one method, we set to 1

Incidence_matrix_merged <- as.data.frame(Incidence_matrix_merged)
rownames(Incidence_matrix_merged) <- Incidence_matrix_merged$site
Incidence_matrix_merged <- Incidence_matrix_merged[,-1]

Incidence_matrix_eDNA <- read.csv("data/incidence_matrix_eDNA.csv", row.names = 1)
Incidence_matrix_trawl<- read.csv("data/incidence_matrix_trawl.csv", row.names = 1)

Reftax <- read.csv("data/reftax.csv", row.names = 1)
Reftax_species <- read.csv("data/reftax_species.csv", row.names = 1)
Traits <- read.csv("data/traits.csv", row.names = 1)

for (trt in c("Env_2", "Body_shape", "mode", "fert", "care")) {
  Traits[, trt] <- as.factor(Traits[, trt])
}
rm(trt)
Traits$Depth.Range <- as.numeric(Traits$Depth.Range) 
Traits_cat <- read.csv("data/traits_cat.csv", row.names = 1)

Trees <- readRDS("data/trees.rds")
EDNA_2019_reads <- read.csv("data/eDNA_2019_reads.csv", row.names = 1)
EDNA_metadata <- read.csv("data/eDNA_metadata.csv", row.names = 1)
TRAWL_metadata <- read.csv("data/trawl_metadata.csv", row.names = 1)

Trawl_2019_abundance <- read.csv("data/trawl_2019_abundance.csv", row.names = 1)

Name_stations_tab <- read.csv("data/name_stations.csv")
Name_stations <- Name_stations_tab[, 2]
names(Name_stations) <- Name_stations_tab[, 1]
rm(Name_stations_tab)


GISdata <- readRDS("data/GISdata.rds")