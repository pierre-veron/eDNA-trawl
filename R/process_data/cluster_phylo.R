### ----------------------------------------------------------------------------
#   Calculates the phylogenetic signal of taxat detected by eDNA and trawling
#   using the D-value measure (Fritz 2010)
#   It makes 100 random draws of species for taxa at genus level and repeats
#   it for each phylogenetic tree. 
#   
#   The mean and sd of values of D obtained and of associated p-values is then 
#   calculated.

#   Reference:
#   Fritz, S. A., & Purvis, A. (2010). Selectivity in mammalian extinction risk 
#     and threat types: a new measure of phylogenetic signal strength in binary 
#     traits. Conservation Biology, 24(4), 1042-1051.
#      https://doi.org/10.1111/j.1523-1739.2010.01455.x 
#   
#
# Author : Pierre Veron, pierre.veron.2017@polytechnique.org
# Created 2022/05/03
### ----------------------------------------------------------------------------

# PART 1 : gamma ---------------------------------------------------------------

if (file.exists("output/results_phylo_signal.rds")) {
  results_phylo_signal <- readRDS("output/results_phylo_signal.rds")
} else {
  n_random_draw <- 4 #100 #  ~1 minute for 1 
  n_tree <- length(Trees)
  stations_eDNA <- paste0(names(Name_stations), "e")
  stations_trawl <- paste0(names(Name_stations), "t")
  results_phylo_signal <- as.data.frame(t(sapply(1:(n_random_draw*n_tree), function(j) {
    df <- Reftax 
    
    # Random draw of species by genus
    df$rep_species <- as.vector(sapply(1:nrow(df), function(i) {
      
      if (df[i, "rank"] == "species") {
        gsub(" ", "_", df[i, "species"])
      } else {
        gen <- df[i, "genus"]
        spec <- sample(Reftax_species[which(Reftax_species$genus == gen), "species"], size = 1)
        gsub(" ", "_", spec)
      }
    }))
    
    rownames(df) <- df$rep_species
    
    df$det_eDNA <- as.vector(sapply(df$taxname, function(tax) {
      if (sum(Incidence_matrix[stations_eDNA, gsub("_", ".", tax)]) > 0) {
        1
      } else {
        0
      }
    }))
    
    df$det_trawl <- as.vector(sapply(df$taxname, function(tax) {
      if (sum(Incidence_matrix[stations_trawl, gsub("_", ".", tax)]) > 0) {
        1
      } else {
        0
      }
    }))
    df <- df[, c("rep_species", "det_eDNA", "det_trawl")]
    df$det_eDNA <- as.numeric(df$det_eDNA)
    df$det_trawl <- as.numeric(df$det_trawl)
    # Drop missing species
    tree <- Trees[[1+(j-1) %/% n_random_draw]]
    tree <- di2multi(tree)
    # Calculate D value
    edna_phyl.d.results <- caper::phylo.d(data = df, phy=  tree, names.col = rep_species, binvar = det_eDNA)
    trawl_phyl.d.results <- caper::phylo.d(data = df, phy=  tree, names.col = rep_species, binvar = det_trawl)
    phyl.d.results <- list("eDNA" = edna_phyl.d.results, "trawl" = trawl_phyl.d.results)
    c("eDNA.D-value" = phyl.d.results$eDNA$DEstimate, 
      "eDNA.p_random" = phyl.d.results$eDNA$Pval1, 
      "eDNA.p_brown" = phyl.d.results$eDNA$Pval0,
      "trawl.D-value" = phyl.d.results$trawl$DEstimate, 
      "trawl.p_random" = phyl.d.results$trawl$Pval1, 
      "trawl.p_brown" = phyl.d.results$trawl$Pval0)
  })))
  
  
  saveRDS(results_phylo_signal,"output/results_phylo_signal.rds")
}

# calculate mean and sd
results_phylo_signal.mean <- colMeans(results_phylo_signal)
results_phylo_signal.sd <- as.vector(t(colSds(as.matrix(results_phylo_signal))))
names(results_phylo_signal.sd) <- names(results_phylo_signal.mean)
results_phylo_signal <- rbind(results_phylo_signal.mean, results_phylo_signal.sd)
rownames(results_phylo_signal) <- c("mean", "sd")
write.csv2(results_phylo_signal, "diversity_indices/phylo_signal_gamma.csv")

# PART 2 : alpha ---------------------------------------------------------------

if (file.exists("output/results_phylo_signal_alpha.rds")) {
  results_phylo_signal_alpha <- readRDS("output/results_phylo_signal_alpha.rds")
} else {
  n_random_draw <- 4 #100 # ~1 minute for 1 and by site
  n_tree <- length(Trees)
  
  results_phylo_signal_alpha <- as.data.frame(t(sapply(1:length(Name_stations), function(i_site) { # 
    site <- Name_stations[i_site]
    print(paste0("Site : ", as.numeric(site)))
    df <- Reftax 
    station <- names(Name_stations)[i_site]
    print(station)
    df$det_site_eDNA <- as.vector(sapply(df$taxname, function(tax) {
      Incidence_matrix[paste0(station, "e"),  gsub(" ", ".", tax)]
    }))
    
    df$det_site_trawl <- as.vector(sapply(df$taxname, function(tax) {
      Incidence_matrix[paste0(station, "t"),  gsub(" ", ".", tax)]
    }))
    
    results_phylo_signal <- as.data.frame(t(sapply(1:(n_random_draw*n_tree), function(j) {
      
      
      # Random draw of species by genus
      df$rep_species <- as.vector(sapply(1:nrow(df), function(i) {
        
        if (df[i, "rank"] == "species") {
          gsub(" ", "_", df[i, "species"])
        } else {
          gen <- df[i, "genus"]
          spec <- sample(Reftax_species[which(Reftax_species$genus == gen), "species"], size = 1)
          gsub(" ", "_", spec)
        }
      }))
      
      rownames(df) <- df$rep_species
      ddf <- df[, c("rep_species", "det_site_eDNA", "det_site_trawl")]
      # Drop missing species
      tree <- Trees[[1+(j-1) %/% n_random_draw]]
      tree <- di2multi(tree)
      # Calculate D value
      edna_phyl.d.results <- caper::phylo.d(data = ddf, phy=  tree, names.col = rep_species, binvar = det_site_eDNA)
      trawl_phyl.d.results <- caper::phylo.d(data = ddf, phy=  tree, names.col = rep_species, binvar = det_site_trawl)
      phyl.d.results <- list("eDNA" = edna_phyl.d.results, "trawl" = trawl_phyl.d.results)
      c("eDNA.D-value" = phyl.d.results$eDNA$DEstimate, 
        "eDNA.p_random" = phyl.d.results$eDNA$Pval1, 
        "eDNA.p_brown" = phyl.d.results$eDNA$Pval0,
        "trawl.D-value" = phyl.d.results$trawl$DEstimate, 
        "trawl.p_random" = phyl.d.results$trawl$Pval1, 
        "trawl.p_brown" = phyl.d.results$trawl$Pval0)
    })))
    
    # calculate mean and sd
    results_phylo_signal.mean <- colMeans(results_phylo_signal)
    results_phylo_signal.sd <- as.vector(t(colSds(as.matrix(results_phylo_signal))))
    names(results_phylo_signal.sd) <- paste0(names(results_phylo_signal.mean), ".sd")
    output <- c(results_phylo_signal.mean, results_phylo_signal.sd)
    saveRDS(output, here("results", "phylo_signal", paste0("temp_site_", site, ".rds")))
    output
  })))
  
  saveRDS(results_phylo_signal_alpha, "output/results_phylo_signal_alpha.rds")
}

write.csv(results_phylo_signal_alpha, "diversity_indices/phylo_signal_alpha.csv")
