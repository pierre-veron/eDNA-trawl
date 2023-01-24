### ----------------------------------------------------------------------------
#   This script calculates the main indices of taxonomic, functional and
#    phylogenetic diversity and generate null models (+ SES).
#
#   The results are calculated once and stored into results/ directory.
#   If the files already exist, they are only loaded. To run the analysis, remove 
#    the associated results from the folder. 
#
#   Author : Pierre Veron, from Romane Rozanski's pipelines
#   pierre.veron.2017@polytechnique.org
#   Created 2022/04/01
### ----------------------------------------------------------------------------

# Parameters
nb_rd <- 4 # number of random sampling for functional
n_rand_ses <- 4

# PART I : taxonomic ---------------------------------------------------------------------

## STEP 1 : eDNA ---------------------------------------------------------------
# Formating dataframes 
# adding site to colnames
Tab_edna <- EDNA_2019_reads
colnames(Tab_edna) <- as.vector(sapply(colnames(Tab_edna), function(cl) {
  paste0(cl, "_", EDNA_metadata[which(EDNA_metadata$filter == cl), "station"])
}))



Tab_edna$Taxon <- rownames(Tab_edna)

# adding columns 
Tab_edna <- Tab_edna %>% 
  left_join(Reftax %>% select(taxname, class, order, family)
            , by = c("Taxon" = "taxname")) %>%
  rename(Class = class, Order = order, Family = family) %>%
  select(Class, Order, Family, Taxon, everything())

# Calculate taxonomic diversity indices 
path_result <- "output/RESULTS_tax_eDNA.rds"
if (!file.exists(path_result)) {
  RESULTS_tax_eDNA <- taxonomic_diversity(Tab_edna, abundance = F, has.replicates = T)
  saveRDS(RESULTS_tax_eDNA, path_result)
} else {
  RESULTS_tax_eDNA <- readRDS(path_result)
}

## STEP 2 : trawling -----------------------------------------------------------
# Formating dataframes 

Tab_trawl <- Trawl_2019_abundance


Tab_trawl$Taxon <- rownames(Tab_trawl)
# adding columns 
Tab_trawl <- Tab_trawl %>% 
  left_join(Reftax %>% select(taxname, class, order, family)
            , by = c("Taxon" = "taxname")) %>%
  rename(Class = class, Order = order, Family = family) %>%
  select(Class, Order, Family, Taxon, everything())

path_result <- "output/RESULTS_tax_trawl.rds"
if (!file.exists(path_result)) {
  RESULTS_tax_trawl <- taxonomic_diversity(Tab_trawl, abundance = F, has.replicates = F)
  saveRDS(RESULTS_tax_trawl, path_result)
} else { 
  RESULTS_tax_trawl <- readRDS(path_result)
}

## STEP 3 : merged ---------------------------------------------------------------
Tab_merged <- as.data.frame(t(Incidence_matrix_merged))

Tab_merged$Taxon <- rownames(Tab_merged)

# adding columns 
Tab_merged <- Tab_merged %>% 
  left_join(Reftax %>% select(taxname, class, order, family)
            , by = c("Taxon" = "taxname")) %>%
  rename(Class = class, Order = order, Family = family) %>%
  select(Class, Order, Family, Taxon, everything()) 

path_result <- "output/RESULTS_tax_merged.rds"
if (!file.exists(path_result)) {
  RESULTS_tax_merged <- taxonomic_diversity(Tab_merged, abundance = F, has.replicates = F)
  
  saveRDS(RESULTS_tax_merged, path_result)
} else {
  RESULTS_tax_merged <- readRDS(path_result)
}

# PART II : functional --------------------------------------------------------- 
## Normalize traits table ------------------------------------------------------
Traits_norm <- Traits
for (cl in colnames(Traits_norm)) {
  if (class(Traits_norm[, cl]) == "numeric") {
    Traits_norm[, cl] <- (Traits[, cl] - mean(Traits[, cl], na.rm = T)) / sd(Traits[, cl], na.rm = T)
  }
}

# same order as Reftax_species
Traits_norm <- Traits_norm[Reftax_species$species,]

## Classification --------------------------------------------------------------
Classif <- Reftax_species[, c('class', 'order', 'family', 'genus', 'species')]
rownames(Classif) <- Reftax_species$species
colnames(Classif) <- c("Class", "Order", "Family", "Genus", "Species")

## STEP 1 : eDNA and Trawling, -------------------------------------------------
Tab_both <- as.data.frame(t(rbind(Incidence_matrix_eDNA, Incidence_matrix_trawl)))

Tab_both$Taxon <- as.vector(sapply(row.names(Tab_both), function(taxname) {
  gsub('\\.', " ", taxname)
}))

# adding columns 
Tab_both <- Tab_both %>% 
  left_join(Reftax %>% select(taxname, class, order, family)
            , by = c("Taxon" = "taxname")) %>%
  rename(Class = class, Order = order, Family = family) %>%
  select(Class, Order, Family, Taxon, everything())

path_result <- "output/RESULTS_fct.rds"
if (!file.exists(path_result)) {
  RESULTS_fct_both <- functional_diversity(tab_eDNA = Tab_both, 
                                           sp_tr_all_sp = Traits_norm,
                                           sp_tr_cat = Traits_cat, 
                                           Classif = Classif, 
                                           abundance = F,
                                           axes4beta = c(1,2),
                                           nb_rd = nb_rd,
                                           save_fct_space = "output/fct_space.rds",
                                           has.replicates = F,
                                           stop_if_NA = F,
                                           display_progress = F, 
                                           do.hill.chao = F,
                                           do.beta = F)
  
  
  saveRDS(RESULTS_fct_both, path_result)
} else {
  RESULTS_fct_both <- read_rds(path_result)
}


## STEP 2 : creation of null model for SES -------------------------------------
path_result <- "output/RESULTS_null_fct_both.rds"
if (!file.exists(path_result)) {
  
  
  RESULTS_null_fct_both <- lapply(1:n_rand_ses, function(i) {
    print(paste0("Randomisation : ", i, " / ", n_rand_ses))
    Traits_random <- randomize(Traits_norm)
    functional_diversity(tab_eDNA = Tab_both, 
                         sp_tr_all_sp = Traits_random,
                         sp_tr_cat = Traits_cat, 
                         Classif = Classif, 
                         abundance = F,
                         axes4beta = c(1,2),
                         nb_rd = nb_rd,
                         save_fct_space = NULL,
                         has.replicates = F,
                         stop_if_NA = F,
                         display_progress = F, 
                         do.hill.chao = F,
                         do.beta = F)
  })
  
  saveRDS(RESULTS_null_fct_both, path_result)
} else {
  RESULTS_null_fct_both <- readRDS(path_result)
}

# Mean alpha
sites <- rownames((RESULTS_null_fct_both[[1]])$Alpha)

RESULTS_null_mean_fct_both <- as.data.frame(t(sapply(sites, function(st) {
  rowMeans(sapply(RESULTS_null_fct_both, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))

RESULTS_rdm_sd_fct_both <- as.data.frame(t(sapply(sites, function(st) {
  rowSds(sapply(RESULTS_null_fct_both, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))
colnames(RESULTS_rdm_sd_fct_both) <- colnames(RESULTS_null_mean_fct_both)

RESULTS_SES_fct_both <- (RESULTS_fct_both$Alpha - RESULTS_null_mean_fct_both) / RESULTS_rdm_sd_fct_both
saveRDS(RESULTS_SES_fct_both, "output/RESULTS_SES_fct_both.RDS")


## STEP 3 : merged -------------------------------------------------------------
path_result <- "output/RESULTS_fct_merged.rds"

Tab_merged <- as.data.frame(t(Incidence_matrix_merged))

Tab_merged$Taxon <-  as.vector(sapply(row.names(Tab_merged), function(taxname) {
  gsub('\\.', " ", taxname)
}))
# adding columns 
Tab_merged <- Tab_merged %>% 
  left_join(Reftax %>% select(taxname, class, order, family)
            , by = c("Taxon" = "taxname")) %>%
  rename(Class = class, Order = order, Family = family) %>%
  select(Class, Order, Family, Taxon, everything())

if (!file.exists(path_result)) {
  RESULTS_fct_merged <- functional_diversity(tab_eDNA = Tab_merged, 
                                           sp_tr_all_sp = Traits_norm,
                                           sp_tr_cat = Traits_cat, 
                                           Classif = Classif, 
                                           abundance = F,
                                           axes4beta = c(1,2),
                                           nb_rd = nb_rd,
                                           save_fct_space = "output/fct_space_merged.rds",
                                           has.replicates = F,
                                           stop_if_NA = F,
                                           display_progress = F, 
                                           do.hill.chao = F,
                                           do.beta = F)
  
  
  saveRDS(RESULTS_fct_merged, path_result)
} else {
  RESULTS_fct_merged <- read_rds(path_result)
}


## STEP 4 : null model merged --------------------------------------------------
path_result <- "output/RESULTS_null_fct_merged.rds"
if (!file.exists(path_result)) {
  
  
  RESULTS_null_fct_merged <- lapply(1:n_rand_ses, function(i) {
    print(paste0("Randomisation : ", i, " / ", n_rand_ses))
    Traits_random <- randomize(Traits_norm)
    functional_diversity(tab_eDNA = Tab_merged, 
                         sp_tr_all_sp = Traits_random,
                         sp_tr_cat = Traits_cat, 
                         Classif = Classif, 
                         abundance = F,
                         axes4beta = c(1,2),
                         nb_rd = nb_rd,
                         save_fct_space = NULL,
                         has.replicates = F,
                         stop_if_NA = F,
                         display_progress = F, 
                         do.hill.chao = F,
                         do.beta = F)
  })
  
  saveRDS(RESULTS_null_fct_merged, path_result)
} else {
  RESULTS_null_fct_merged <- readRDS(path_result)
}

# Mean alpha
sites <- rownames((RESULTS_null_fct_merged[[1]])$Alpha)

RESULTS_null_mean_fct_merged <- as.data.frame(t(sapply(sites, function(st) {
  rowMeans(sapply(RESULTS_null_fct_merged, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))

RESULTS_rdm_sd_fct_merged <- as.data.frame(t(sapply(sites, function(st) {
  rowSds(sapply(RESULTS_null_fct_merged, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))
colnames(RESULTS_rdm_sd_fct_merged) <- colnames(RESULTS_null_mean_fct_merged)

RESULTS_SES_fct_merged <- (RESULTS_fct_merged$Alpha - RESULTS_null_mean_fct_merged) / RESULTS_rdm_sd_fct_merged
saveRDS(RESULTS_SES_fct_merged, "output/RESULTS_SES_fct_merged.rds")

# PART III : phylo -------------------------------------------------------------
# remove one species that is not on the phylo tree
Classif <- Classif[which(rownames(Classif) != "Chelidonichthys lastoviza"), ]

## STEP 1 :  Phylogenetic indices ----------------------------------------------
path_result <- "output/RESULTS_phy.rds"
if (!file.exists( path_result)) {
  RESULTS_phyl_both <- phylo_diversity(trees = Trees,
                                       tab_eDNA = Tab_both,
                                       Classif = Classif, 
                                       savepathindic ="output/trees/",
                                       abundance = F, 
                                       has.replicates = F, 
                                       do.beta = F)
  saveRDS(RESULTS_phyl_both, path_result)
} else {
  RESULTS_phyl_both <- readRDS(path_result)
}

## STEP 2 : null model ---------------------------------------------------------
path_result <- "output/RESULTS_null_phy.rds"
if (!file.exists(path_result)) {
  
  
  n_sp_tree <- length(Trees[[1]]$tip.label)
  RESULTS_null_phyl_both <- lapply(1:n_rand_ses, function(i) {
    print(paste0("Randomisation : ", i, " / ", n_rand_ses))
    spl <- sample(1:n_sp_tree, size = n_sp_tree, replace = F)
    Trees_random <- lapply(Trees, function(tree) {
      random_tree <- tree
      random_tree$tip.label <- tree$tip.label[spl]
      random_tree
    })
    
    phylo_diversity(trees = Trees_random,
                    tab_eDNA = Tab_both,
                    Classif = Classif, 
                    savepathindic = "output/trees_null/",
                    abundance = F, 
                    has.replicates = F, 
                    do.beta = F)
  })
  saveRDS(RESULTS_null_phyl_both, path_result)
} else {
  RESULTS_null_phyl_both <- readRDS(path_result)
}

# Calculate means and SD on null model
sites <- rownames((RESULTS_null_phyl_both[[1]])$Alpha)

RESULTS_null_mean_phyl_both <- as.data.frame(t(sapply(sites, function(st) {
  rowMeans(sapply(RESULTS_null_phyl_both, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))

RESULTS_rdm_sd_phyl_both <- as.data.frame(t(sapply(sites, function(st) {
  rowSds(sapply(RESULTS_null_phyl_both, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))
colnames(RESULTS_rdm_sd_phyl_both) <- colnames(RESULTS_null_mean_phyl_both)

RESULTS_SES_phyl_both <- (RESULTS_phyl_both$Alpha - RESULTS_null_mean_phyl_both) / RESULTS_rdm_sd_phyl_both
saveRDS(RESULTS_SES_phyl_both, "output/RESULTS_SES_phy_both.RDS")

# SES on gamma 
df_null_models <- do.call(rbind, 
                          lapply(RESULTS_null_phyl_both, function(RESULTS) {
                            a <- RESULTS$Gamma
                            a
                          }))
RESULTS_null_mean_phyl_gamma <- colMeans(df_null_models)
RESULTS_null_sd_phyl_gamma <- (colSds(as.matrix(df_null_models)))
names(RESULTS_null_sd_phyl_gamma) <- names(RESULTS_null_mean_phyl_gamma)
RESULTS_SES_phyl_gamma <- (RESULTS_phyl_both$Gamma - RESULTS_null_mean_phyl_gamma) / RESULTS_null_sd_phyl_gamma 

saveRDS(RESULTS_SES_phyl_gamma, "output/RESULTS_SES_phy_gamma.rds")

## STEP 3 : merged -------------------------------------------------------------
path_result <- "output/RESULTS_phy_merged.rds"
if (!file.exists( path_result)) {
  RESULTS_phyl_merged <- phylo_diversity(trees = Trees,
                                       tab_eDNA = Tab_merged,
                                       Classif = Classif, 
                                       savepathindic = "output/trees_merged/",
                                       abundance = F, 
                                       has.replicates = F, 
                                       do.beta = F)
  saveRDS(RESULTS_phyl_merged, path_result)
} else {
  RESULTS_phyl_merged <- readRDS(path_result)
}

## STEP 4 : null model merged --------------------------------------------------
path_result <- "output/RESULTS_null_phy_merged.rds"
if (!file.exists(path_result)) {
  
  n_sp_tree <- length(Trees[[1]]$tip.label)
  RESULTS_null_phyl_merged <- lapply(1:n_rand_ses, function(i) {
    print(paste0("Randomisation : ", i, " / ", n_rand_ses))
    spl <- sample(1:n_sp_tree, size = n_sp_tree, replace = F)
    Trees_random <- lapply(Trees, function(tree) {
      random_tree <- tree
      random_tree$tip.label <- tree$tip.label[spl]
      random_tree
    })
    
    phylo_diversity(trees = Trees_random,
                    tab_eDNA = Tab_merged,
                    Classif = Classif, 
                    savepathindic = "output/trees_null_merged/",
                    abundance = F, 
                    has.replicates = F, 
                    do.beta = F)
  })
  saveRDS(RESULTS_null_phyl_merged, path_result)
} else {
  RESULTS_null_phyl_merged <- readRDS(path_result)
}

# Calculate means and SD on null model
sites <- rownames((RESULTS_null_phyl_merged[[1]])$Alpha)

RESULTS_null_mean_phyl_merged <- as.data.frame(t(sapply(sites, function(st) {
  rowMeans(sapply(RESULTS_null_phyl_merged, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))

RESULTS_rdm_sd_phyl_merged <- as.data.frame(t(sapply(sites, function(st) {
  rowSds(sapply(RESULTS_null_phyl_merged, function(RESULTS) {
    Alpha_rdm <- RESULTS$Alpha
    as.matrix(Alpha_rdm[st,])[1,]
  }))
})))
colnames(RESULTS_rdm_sd_phyl_merged) <- colnames(RESULTS_null_mean_phyl_merged)

RESULTS_SES_phyl_merged <- (RESULTS_phyl_merged$Alpha - RESULTS_null_mean_phyl_merged) / RESULTS_rdm_sd_phyl_merged
saveRDS(RESULTS_SES_phyl_merged,"output/RESULTS_SES_phy_merged.RDS")

# SES on gamma 
df_null_models <- do.call(rbind, 
                          lapply(RESULTS_null_phyl_merged, function(RESULTS) {
                            a <- RESULTS$Gamma
                            a
                          }))
RESULTS_null_mean_phyl_gamma_merged <- colMeans(df_null_models)
RESULTS_null_sd_phyl_gamma_merged <- (colSds(as.matrix(df_null_models)))
names(RESULTS_null_sd_phyl_gamma_merged) <- names(RESULTS_null_mean_phyl_gamma_merged)
RESULTS_SES_phyl_gamma_merged <- (RESULTS_phyl_merged$Gamma - RESULTS_null_mean_phyl_gamma_merged) / RESULTS_null_sd_phyl_gamma_merged 
saveRDS(RESULTS_SES_phyl_gamma_merged, "output/RESULTS_SES_phy_gamma_merged")


# END - remove some variables --------------------------------------------------
rm(Tab_edna, Tab_trawl, Tab_both, Tab_merged, Classif, df_null_models, 
   cl, path_result,sites, nb_rd, n_rand_ses)
