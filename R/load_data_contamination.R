### ----------------------------------------------------------------------------
#    Loads the data used in this project
### ----------------------------------------------------------------------------

# Metadata
EDNA_metadata <- read.csv("data/eDNA_metadata.csv", row.names = 1)
TRAWL_metadata <- read.csv("data/trawl_metadata.csv", row.names = 1)
Name_stations_df <- read.csv("data/name_stations.csv", row.names = 1)
Name_stations <- Name_stations_df[, 1]
names(Name_stations) <- rownames(Name_stations_df)
rm(Name_stations_df)

# Incidence matrices
Incidence_matrix <- read.csv("data/incidence_matrix.csv", row.names = 1)

sp2del <- sapply(c("Capros.aper","Trachurus","Argentina","Conger.conger","Dicentrarchus","Sardina.pilchardus","Engraulis",
                   "Micromesistius.poutassou","Trisopterus","Lophius","Lepidorhombus","Pagellus.bogaraveo"), 
                 function(x)grep(x,colnames(Incidence_matrix)))
Incidence_matrix <-   Incidence_matrix[,-sp2del]

Incidence_matrix_eDNA <- read.csv("data/incidence_matrix_eDNA.csv", row.names = 1)
sp2del <- sapply(c("Capros.aper","Trachurus","Argentina","Conger.conger","Dicentrarchus","Sardina.pilchardus","Engraulis",
                   "Micromesistius.poutassou","Trisopterus","Lophius","Lepidorhombus","Pagellus.bogaraveo"), 
                 function(x)grep(x,colnames(Incidence_matrix_eDNA)))
Incidence_matrix_eDNA <-   Incidence_matrix_eDNA[,-sp2del]


Incidence_matrix_trawl <- read.csv("data/incidence_matrix_trawl.csv", row.names = 1)
sp2del <- sapply(c("Capros.aper","Trachurus","Argentina","Conger.conger","Dicentrarchus","Sardina.pilchardus","Engraulis",
                   "Micromesistius.poutassou","Trisopterus","Lophius","Lepidorhombus","Pagellus.bogaraveo"), 
                 function(x)grep(x,colnames(Incidence_matrix_trawl)))
Incidence_matrix_trawl <-   Incidence_matrix_trawl[,-sp2del]

Incidence_matrix_merged <- Incidence_matrix_eDNA
rownames(Incidence_matrix_merged) <- gsub("e", "", rownames(Incidence_matrix_merged) )

for (station in rownames(Incidence_matrix_merged)) {
  for (taxa in colnames(Incidence_matrix_merged)) {
    if (Incidence_matrix_trawl[paste0(station, 't'), taxa] > 0) {
      Incidence_matrix_merged[station, taxa] <- 1
    }
  }
}
rm(station, taxa)
# Taxonomic references
Reftax <- read.csv("data/reftax.csv", row.names = 1)
Reftax_species <- read.csv("data/reftax_species.csv", row.names = 1)

# Functional traits
Traits <- read.csv("data/traits.csv", row.names = 1)
Traits$Depth.Range <- as.numeric(Traits$Depth.Range)
Traits$Env_2 <- as.factor(Traits$Env_2)
Traits$Body_shape <- as.factor(Traits$Body_shape)
Traits$mode <- as.factor(Traits$mode)
Traits$fert <- as.factor(Traits$fert)
Traits$care <- as.factor(Traits$care)
Traits_cat <- read.csv("data/traits_cat.csv", row.names = 1)

# Phylogenetic trees
Trees <- readRDS("data/trees.rds")

# Abundance data
EDNA_2019_reads <- read.csv("data/eDNA_2019_reads.csv", row.names = 1)
Trawl_2019_abundance <- read.csv("data/trawl_2019_abundance.csv", row.names = 1)

# Geographical data (for maps)
GISdata <- readRDS("data/GISdata.rds")
