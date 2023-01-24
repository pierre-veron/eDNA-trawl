####################################################################################################################################################
# Functional pipeline to calculate Functional Diversity gamma,alpha and beta indices and functional Hill numbers
#based on a table of functional traits and a eDNA table 

#Authors: Romane Rozanski, Camille Albouy, David Eme


#inputs:
# tab_eDNA: Data frame containing the number of copies/reads for each taxa detected with eDNA (filters in columns, taxa in rows)
# sp_tr_all_sp: Data frame containing the functional traits for species related to each taxa detected with eDNA (traits in columns, species in rows)
#WARNING: traits need to be already sorted and standardized 
#sp_tr_cat: Data frame containing the type (continuous, ordinal, categorical, binary) and the weight of each trait from the sp_tr_all_sp table
#Classif: Taxonomic classification of each species from sp_tr_all_sp (Class, Order,Family, Genus)
#savepath: path to save some results from the different functions
#abundance =TRUE/FALSE: Choice if the indices are computed using species abundance (here the number of copies/reads from the eDNA table) or presence/absence
#axes4beta: number of axes used for computing beta-diversity indices (warning: time-consuming step the more axes you choose, the longer it takes)
#nb_rd: Number of randomisation to perform  to randomly choose a representative species for taxa detected at the genus/family level (here, default: 100)

#returns: 
#a list divided in four lists containing the alpha-diversity indices, beta-diversity indices, gamma-diversity indices and Hill numbers (from Chao et al. 2019) )


#Last update: 23/03/2022 (by Romane Rozanski)
####################################################################################################################################################

#fonctions source needed in the pipeline
source(here("R", "functions", "fonctions_source_f", "Functional_&_Phylogenetic.Beta.MPFD.MNND.R"))
source(here("R", "functions", "fonctions_source_f", "functional_hill_number_chao.R")) # function to compute the functional hill numbers (by Pierre Veron: pierre.veron.2017@polytechnique.org; See Chao et al. 2019)
source(here("R", "functions", "fonctions_source_f", "get_beta_val.R"))
source(here("R", "functions", "fonctions_source_f", "Funct_indices_DavidEme_V2.R"))#function to compute functional MPD,MNTD and their respective variance (by David Eme : David.eme3@gmail.com)
# here needed to compute Var of Functional Pairwise Distance (VFPD) as alpha regularity indice because not computed in the "mFD" package. 




################################### PIPELINE ###########################################################################

#' functional_diversity
#'
#' Functional pipeline to calculate Functional Diversity gamma,alpha and beta 
#' indices and functional Hill numbers based on a table of functional traits and
#' a eDNA table 
#' @param tab_eDNA a dataframe containing abundance of species by filter (site 
#' as rows, species as columns) and taxonomy. 
#' The format must be:
#' Class     Order     Family     Taxon          Filt1_Site1   Filt2_Site1   Filt3_Site2   ....
#' Actino    Clupei    Engrau     Engrau encra     NA             13200         58000
#' ...
#' If has.replicates : the name of the columns shound be formated like this : 
#' namefilt_namesite 
#' and the function will add together the abundances of filters for a same site.
#' If not has.replicates : the name of the columns is free.
#' @param sp_tr_all_sp Data frame containing the functional traits for species 
#' related to each taxa detected with eDNA (traits in columns, species in rows)
#   WARNING: traits need to be already sorted and standardized 
#' @param sp_tr_cat Data frame containing the type (continuous, ordinal, 
#' categorical, binary) and the weight of each trait from the sp_tr_all_sp table
#' Format must be as the tr_cat variable
#' See ?mFD::sp.tr.summary for the explicit formating of the traits category
#' table.
#' @param Classif Taxonomic classification of each species from sp_tr_all_sp
#'  (Class, Order,Family, Genus)
#' @param abundance TRUE/FALSE: Choice if the indices are computed using species
#' abundance (here the number of copies/reads from the eDNA table) or presence/absence
#' Default : F
#' @param axes4beta number of axes used for computing beta-diversity indices 
#' (warning: time-consuming step the more axes you choose, the longer it takes)
#' Default : c(1,2)
#' @param nb_rd Number of randomisation to perform  to randomly choose a 
#' representative species for taxa detected at the genus/family level 
#' (here, default: 100)
#' @param save_fct_space NULL, or a string with path to a .RDS file where the 
#' functional spaces from the randomisation will be stored.
#' @param has.replicates a boolean (default TRUE) : T if there are > 1 replicates
#' by site (in this case, the function will add abundances by site). 
#' If F each column is treated as a independant site. 
#' @param stop_if_NA (default T) : stop analysis if one of traits is NA 
#' @param display_progress (default T) : shows the progress of the calculation
#' of alpha indices in the console.
#' @param  do.beta (default T) : calculates beta indices 
#' @param do.hill.chao (default T) : calculates Functional hill numbers acc. to Chao 2019
#' @return a list divided in four lists containing the alpha-diversity indices, 
#' beta-diversity indices, gamma-diversity indices and Hill numbers 
#' (from Chao et al. 2019))
#' 
#' @author Romane Rozanski, Camille Albouy, David Eme
#' @note Updated : 2022/03/30 by Pierre Veron. Corrections : 
#'       - the condition has.replicates and save_fct_space
#'       - adding the possibility that none of taxa is detected at family level
#'       - add the parameter stop_if_NA (if set to default, should not change anything)
#'       - remove parameter savepath which seems useless
#'       - add the parameter save_fct_space (to save the functional spaces)
#'       - add a display_progress parameter (to avoid filling the console with too 
#'          much information)
#'       - some minor corrections 
#' @note Minor update 2022/06/24, add the parameters do.beta and do.hill.chao
#' @details From mFD publication : list of all computed indices with signification
#' Functional richness 
#'     (Cornwell et al. 2006, Villéger et al. 2008)	
#'     FRic	
#'     The volume 
#'     of the convex hull shaping the species present in the assemblage	Computed 
#'     using the ‘convhulln' function of the ‘geometry' package. Computed only if 
#'     the number of species (i.e. points) is strictly higher than the number of 
#'     functional axes. If points are coplanar, the convex hull can not be computed 
#'     and the function returns NA.
#' Functional identity 
#'     (Garnier et al. 2004, Mouillot et al. 2013a)	
#'     FIde	
#'     The weighted average position of species of the assemblage along each axis
#'     	None
#' Functional dispersion 
#'     (Laliberté and Legendre 2010)	
#'     FDis	
#'     The weighted deviation to center of gravity (i.e. defined by the FIde 
#'     values) of species in the assemblage	FIde is always computed with FDis.
#' Functional divergence 
#'     (Villéger et al. 2008)	
#'     FDiv	
#'     The deviation of biomass-density to the center of gravity of the vertices
#'      shaping the convex hull of the studied assemblage	FDiv requires computing
#'       first vertices of the convex hull so it has the same constraints as FRic
#'        (see above).
#' Functional evenness 
#'     (Villéger et al. 2008)	
#'     FEve	
#'     The regularity of biomass-density distribution along the minimum spanning 
#'     tree (i.e. the tree linking all species of the assemblage with the lowest cumulative branch length) for the studied assemblage	There must be at least three species to compute the FEve index. The minimum spanning tree is computed using the ‘mst' function of the ‘ape' package.
#' Functional originality 
#'     (Mouillot et al. 2013a)	
#'     FOri	
#'     The weighted mean distance to the nearest species from the global species
#'      pool	None
#' Functional specialisation 
#'     (Bellwood et al. 2006, Mouillot et al. 2013a)	
#'     FSpe	
#'     The weighted mean distance to the centroid of the global species pool 
#'     (i.e. center of the functional space)	None
#' Functional mean pairwise distance 
#'     (Weiher et al. 1998)	
#'     FMPD	
#'     The mean weighted distance between all pairs of species	None
#' Functional mean nearest neighbor distance 
#'     (Weiher et al. 1998)	
#'     FNND	
#'     The weighted distance to the nearest neighbor within the assemblage	 
functional_diversity <- function(tab_eDNA,sp_tr_all_sp, sp_tr_cat,
                                 Classif=Classif_end,
                                 abundance=FALSE,axes4beta=c(1:2),
                                 nb_rd=100, save_fct_space = NULL,
                                 has.replicates = T, stop_if_NA = T,
                                 display_progress = T,
                                 do.beta = T,
                                 do.hill.chao = T){
  
  
  
  ############################################################################################################################
  ###PART 1 : Classification of species and 100 Randomisation of species belonging to the same genus/family that have been detetected at these 2 levels by eDNA
  ## (entry : sp_tr_all_sp, tab_eDNA) 
  
  sp_tr_all_classif <- (cbind(Classif,sp_tr_all_sp)) # PV : na.omit(cbind(Classif,sp_tr_all_sp))
  ##searching for fish detected at the FAMILY level by spygen in the traits table and randomisation 100 times
  sp_tr_fam <- sp_tr_all_classif[which(sp_tr_all_classif$Family %in% tab_eDNA$Taxon),]
  fam_list <- unique(sp_tr_fam$Family)

  list_100_fam_random <- list()
  if (length(fam_list) > 0) {
    for (i in (1:nb_rd)){ #randomisation 100 times
      list_100_fam_random[[i]] <- do.call(rbind,lapply(seq(1,length(fam_list)), function(x){
        DFtempt <- sp_tr_fam[which(sp_tr_fam$Family==fam_list[x]),]
        DFtempt <- DFtempt[sample(1:nrow(DFtempt), 1, replace=FALSE), ]
        row.names (DFtempt) <- fam_list[x]; DFtempt 
      }))
    } #end of i
  } else {
    list_100_fam_random <- lapply(1:nb_rd, function(i) {NULL})
  }
  ##searching for fish detected at GENUS level detectetd by Spygen and randomisation 
  sp_tr_genus <- sp_tr_all_classif[which(sp_tr_all_classif$Genus %in%tab_eDNA$Taxon & !sp_tr_all_classif$Species %in% tab_eDNA$Taxon ),] #condition "&" to avoid randomisation of fish already detetected at the species level
  genus_list <- unique(sp_tr_genus$Genus)
  

  list_100_genus_random <- list()
  if (length(genus_list) > 0) {
    for (i in (1:nb_rd)){ #randomisation
      list_100_genus_random[[i]]<- do.call(rbind,lapply(seq(1,length(genus_list)), function(x){
        DFtempt <- sp_tr_genus[which(sp_tr_genus$Genus==genus_list[x]),] 
        DFtempt <- DFtempt[sample(1:nrow(DFtempt), 1, replace=FALSE), ]
        row.names (DFtempt) <- genus_list[x]; DFtempt 
      }))
    }
  } else {
    list_100_genus_random <- lapply(1:nb_rd, function(i) {NULL})
  }
  cat("######## Randomisation ########","\n")
  #Grouping the 100 list of randomisation of family and genus and adding the rest of fish already detected at the species level by spygen
  sp_tr_species_only <- sp_tr_all_classif[which(sp_tr_all_classif$Species %in% tab_eDNA$Taxon),]
  row.names(sp_tr_species_only) <- sp_tr_species_only$Species #put species as rownames (to match with tab_eDNA in function 3) 
  
  list_100_sp_tr_final <- lapply(seq(1,nb_rd), function(x){
    Test <- rbind(list_100_fam_random[[x]], list_100_genus_random[[x]], sp_tr_species_only)
    Test <- Test[order(rownames(Test)),]
    subset(Test, select= c(-Class, -Order, -Family, -Genus, -Species))
  })
  
  
  #(exit : list_100_sp_tr_final)
  #############################################################################################################################
  
  
  #PART 2 : Functional alpha-diversity indices calculation (with mFD package)
  ##(entry : list_100_sp_tr_final, sp_tr_cat, tab_eDNA)
  
  
  ###step 1 : CREATION MATRICE ASB_SP_W
  cat("######## Data preparation for alpha indices calculation ########","\n")
  
  rownames(tab_eDNA)<-tab_eDNA$Taxon
  tab_nb_1 <- tab_eDNA[,5:ncol(tab_eDNA)] #delete all the useless columns of classification
  tab_nb_1 <- tab_nb_1[order(rownames(tab_nb_1)), ] #alphabetic order
  
  if (has.replicates) {
    #grouping the columns (filters) with the same site
    split_col <- strsplit(colnames(tab_nb_1),"_")
    
    Vect_concat <- sapply(seq(1:length(split_col)), function(x){
      split_col <- split_col[[x]][-1]
      str_c(split_col, collapse="_")
    })
    colnames(tab_nb_1) <- Vect_concat
    
    #pooling the same site together and adding the number of copies of each species
    site_list <- (unique(colnames(tab_nb_1))) 
    
    tab_nb_sites <- do.call(cbind,lapply(seq(1,length(site_list)), function(x){
      DFtempt <- tab_nb_1[,which(colnames(tab_nb_1)==site_list[x])]
      if(length(which(colnames(tab_nb_1)==site_list[x]))>1){
        rowSums(DFtempt, na.rm=TRUE) # total sum of copies for species in the same site 
      } else {
        DFtempt
      } 
    })) 
    colnames(tab_nb_sites) <- site_list
  } else {
    tab_nb_sites <- tab_nb_1
  }
  
  asb_sp_w <- t(tab_nb_sites) #transpose the matrice to get species in col and sites in rows
  asb_sp_w[is.na(asb_sp_w)] <- 0 #transformation "NA" in "0"
  
  #calculating the total number of copies by species for the "total site" to do the gamma diversity
  total_site <- colSums(asb_sp_w)
  asb_sp_w <- rbind(asb_sp_w,total_site)
  
  
  ###step 2 : DISTANCE MATRICES AND PCOA
  cat("######## Distance matrices and PCOA  ########","\n")
  
  #computation of the distance matrices
  list_100_sp_dist <-lapply(seq(1, length(list_100_sp_tr_final)), function(x){
    mFD::funct.dist(list_100_sp_tr_final[[x]], sp_tr_cat, metric = "gower", scale_euclid = 'range',
                    stop_if_NA = stop_if_NA)
    })
  
  #computation of PCOA and the functional space quality associated to each cumulated axes
  list_100_fspaces_quality <- lapply(seq(1, length(list_100_sp_dist)), function(x){
    mFD::quality.fspaces(list_100_sp_dist[[x]], maxdim_pcoa = 10, deviation_weighting = "absolute",
                         fdist_scaling = F, fdendro = "average")})
  
  #pcoa coordinates
  list_100_sp_faxes_coord <-lapply(seq(1, length(list_100_fspaces_quality)), function(x){
    list_100_fspaces_quality[[x]]$details_fspaces$sp_pc_coord})
  
  if (!is.null(save_fct_space)) {
    fct_space_with_eigen <- lapply(1:length(list_100_fspaces_quality), function(i) {
      fct_space <- list_100_fspaces_quality[[i]]$details_fspaces$sp_pc_coord
      eigen <- list_100_fspaces_quality[[i]][["details_fspaces"]][["pc_eigenvalues"]][["Relative_eig"]]
      colnames(fct_space) <- as.vector(sapply(1:ncol(fct_space), function(j) {
        paste0(colnames(fct_space)[j], " (", 
               round(100* eigen[j], digits = 1),
               " %)")
      }))
      fct_space
    })
    saveRDS(fct_space_with_eigen, save_fct_space)
  }
  
  ###step 3 : ALPHA INDICES CALCULATION
  cat("######## alpha diversity indices calculation ########","\n")
  
  ## in presence-absence (PA) or abundance (PA means "abundance= FALSE" in entry of the pipeline) 
  if (abundance ==FALSE){
    asb_sp_w <- replace(asb_sp_w,asb_sp_w[,]>0,1) #pres_abs modif
  }
  
  
  #FD ALPHA & GAMMA DIV INDICES (gamma corresponding to the last row of alpha indices: "total_site")
  list_100_alpha_fd_indices <- lapply(seq(1, length(list_100_sp_faxes_coord)), function(x){
    mFD::alpha.fd.multidim(list_100_sp_faxes_coord[[x]][, c("PC1", "PC2","PC3","PC4","PC5")], asb_sp_w,
                           ind_vect = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", "fspe"),scaling=F,
                           check_input=T,details_returned=T, verbose = display_progress)
    })
  
  
  list_100_fd_ind_values <-lapply(seq(1, length(list_100_alpha_fd_indices)), function(x){
    list_100_alpha_fd_indices[[x]]$functional_diversity_indices})
  
  #indices mean in the list of 100 for each site(=alpha) and in total(=gamma)
  sites_list<-rownames(list_100_fd_ind_values[[1]])
  
  average_alpha_fd_indices<- do.call(rbind, lapply(seq(1, length(sites_list)), function(x){
    DFtemp = do.call(rbind, lapply(seq(1, length(list_100_fd_ind_values)), function(i){
      list_100_fd_ind_values[[i]][rownames(list_100_fd_ind_values[[i]])== sites_list[x],]
    }))
    c(Mean=apply(DFtemp, 2, mean, na.rm=T),Sd=apply(DFtemp, 2, sd, na.rm=T))
  }))
  
  ###Var of functional Pairwise Distance calculation (David Eme; see source function) : regularity indice
  list_100_VFPD<-list()
  for (x in 1:length(list_100_sp_dist)){
    list_100_VFPD[[x]]<-data.frame(VFPD=Var.mPwFD(list_100_sp_dist[[x]], asb_sp_w))
  }
  
  
  #VFPD mean
  mean_VFPD<-do.call(rbind, lapply(seq(1, length(sites_list)), function(x){
    DFtemp = do.call(rbind, lapply(seq(1, length(list_100_VFPD)), function(i){
      list_100_VFPD[[i]][rownames(list_100_VFPD[[i]])== sites_list[x],]
    }))
    c(Mean=apply(DFtemp, 2, mean, na.rm=T),Sd=apply(DFtemp, 2, sd, na.rm=T))
  }))
  rownames(mean_VFPD)<-sites_list
  colnames(mean_VFPD)<-c("Mean.VFPD","Sd.VFPD")
  
  
  ###Step 4: Alpha and gamma Hill number computation (with mFD package)
  cat("######## Hill number alpha_gamma ########","\n")
  
  if (abundance == FALSE){
    q_hill=c(0) #only for q = 0: functional richness
  } else{
    q_hill=c(0,1,2) #abundance accounted for (1: functional evenness / 2: functional divergence. See Chao et al. 2019)
  }
  
  list_alpha_hill <- lapply(1:nb_rd, function(x){
    alpha.fd.hill(as.matrix(asb_sp_w),list_100_sp_dist[[x]], q=q_hill,tau="mean")$asb_FD_Hill})
  
  alpha.fd.hill <- array(data=NA,c(nrow(asb_sp_w),length(q_hill),nb_rd)); i=1
  for (i in 1:nb_rd) alpha.fd.hill[,,i] <- as.matrix(list_alpha_hill[[i]])
  
  a_hill_mean <- apply(alpha.fd.hill, c(1,2),mean) 
  a.hill_sd <- apply(alpha.fd.hill, c(1,2),sd)
  rownames(a_hill_mean) <- rownames(a.hill_sd) <- rownames(list_alpha_hill[[1]])
  colnames(a_hill_mean) <- colnames(a.hill_sd) <- colnames(list_alpha_hill[[1]])
  
  alpha.hill_tot<- data.frame(Mean.Hill=a_hill_mean, Sd.Hill=a.hill_sd)
  
  
  
  #reunion of all alpha indices
  average_all_alpha_fd_indices <- cbind(average_alpha_fd_indices,mean_VFPD, alpha.hill_tot)
  
  
  #(exit: average_all_alpha_fd_indices)
  #############################################################################################################################  
  
  
  #PART 3 : Functional beta-diversity indices calculation 
  if (do.beta) {
  cat("######## Beta-diversity indices only in presence/absence ########","\n")
  
  asb_sp_w_beta<-subset(asb_sp_w,row.names(asb_sp_w) != "total_site") #delete this row because only used for gamma div
  #asb_sp_summ_beta <- asb.sp.summary(asb_sp_w_beta) ##needed to compute the next step 
  
  
  if (abundance == FALSE){
    
    #Step 1 : indices computed with mFD package (bjac, bjac_turnover, bjac_nestedness) : WARNING only with species pres/abs
    list_100_beta_fd_indices <- lapply(1:length(list_100_sp_faxes_coord), function(i){
      cat("x:",i,"\n")
      mFD::beta.fd.multidim(list_100_sp_faxes_coord[[i]][, axes4beta],
                            asb_sp_w_beta,betapart_para=T,check_input=F,beta_family=c("Jaccard"),details_returned=T)
    })
    
    list_100_pairasb <- lapply(seq(1, length(list_100_beta_fd_indices)), function(x){
      list_100_beta_fd_indices[[x]]$pairasb_fbd_indices})
    
    
    #beta indices mean for the 100 data frames
    Beta_jac <- get_beta_val( Data= lapply(list_100_pairasb, function (x) x$jac_diss),nb_rd=nb_rd,sites= sites_list,names="jac_diss")
    Beta_turn <- get_beta_val( Data= lapply(list_100_pairasb, function (x) x$jac_turn),nb_rd=nb_rd,sites= sites_list,names="jac_turn")
    Beta_nes <- get_beta_val( Data= lapply(list_100_pairasb, function (x) x$jac_nest),nb_rd=nb_rd,sites= sites_list,names="jac_nest")
    
    average_beta_fd_indices <- cbind(Beta_jac,Beta_turn[,c(3:4)],Beta_nes[,c(3:4)])
    rownames(average_beta_fd_indices)<- paste( average_beta_fd_indices[,1],average_beta_fd_indices[,2],sep="_")
    
    
    #Step 2: others beta- diversity indices ONLY IN PRES_ABS
    
    #Mean Pairwise Functional Distance between communities (MPFD beta)
    cat("######## others beta_indices: beta-MPFD, beta-VPFD (in pres/abs) ########","\n")
    
    list_100_beta_MPFD <- mclapply(seq(1, length(list_100_sp_dist)),mc.cores=1, function(x){
      comdist.Swenson(asb_sp_w_beta, list_100_sp_dist[[x]])})
    
    
    res <- array(data=NA,c(nrow(asb_sp_w_beta),nrow(asb_sp_w_beta),nb_rd)); i=1
    for ( i in 1:nb_rd) res[,,i] <- as.matrix(list_100_beta_MPFD[[i]])
    
    Mean_MPFD <-apply(res, c(1,2),mean) 
    rownames(Mean_MPFD) <- colnames(Mean_MPFD) <- sites_list[-length(sites_list)]
    Mean_MPFD[upper.tri(Mean_MPFD,diag=T)] <- NA
    Mean_MPFD <- na.omit(melt(Mean_MPFD))$value
    
    SD_MPFD <- apply(res, c(1,2),sd)
    rownames(SD_MPFD) <- colnames(SD_MPFD) <- sites_list[-length(sites_list)]
    SD_MPFD[upper.tri(SD_MPFD,diag=T)] <- NA
    SD_MPFD <- na.omit(melt(SD_MPFD))$value
    
    
    #Mean Nearest Functional Distance among nearest neigbors between communities (MNFD beta)
    list_100_beta_MNND<-mclapply(seq(1, length(list_100_sp_dist)),mc.cores=1,function(x){
      comdistnt.Swenson(asb_sp_w_beta, list_100_sp_dist[[x]])})
    
    res <- array(data=NA,c(nrow(asb_sp_w_beta),nrow(asb_sp_w_beta),nb_rd)); i=1
    for ( i in 1:nb_rd) res[,,i] <- as.matrix(list_100_beta_MNND[[i]])
    MNND <- list(Mean =apply(res, c(1,2),mean), SD=apply(res, c(1,2),sd))
    
    Mean_MNND <-apply(res, c(1,2),mean) 
    rownames(Mean_MNND) <- colnames(Mean_MNND) <- sites_list[-length(sites_list)]
    Mean_MNND[upper.tri(Mean_MNND,diag=T)] <- NA
    Mean_MNND <- na.omit(melt(Mean_MNND))$value
    
    SD_MNND <- apply(res, c(1,2),sd)
    rownames(SD_MNND) <- colnames(SD_MNND) <- sites_list[-length(sites_list)]
    SD_MNND[upper.tri(SD_MNND,diag=T)] <- NA
    SD_MNND <- na.omit(melt(SD_MNND))$value
    
  }   #end of "if abundance==false"
  
  
  ### Step 3: Beta-Hill number computation (with mFD package)
  cat("######## Beta hill numbers ########","\n")
  
  if (abundance == FALSE){
    q_hill=c(0) #only for q = 0: functional richness
  } else{
    q_hill=c(0,1,2) #abundance accounted for (1: functional evenness / 2: functional divergence. See Chao et al. 2019)
  }
  
  #Hill calculation for 100 dist objects
  list_beta_hill<- lapply(1:nb_rd, function(x){
    beta.fd.hill(as.matrix(asb_sp_w_beta),list_100_sp_dist[[x]], q=q_hill, tau="mean", beta_type="Jaccard")})
  
  
  if (abundance== FALSE){
    beta.hill <- get_beta_val( Data= lapply(list_beta_hill, function(x) x$beta_fd_q$q0),nb_rd=nb_rd,sites= sites_list,names="Hill_q0")
    
    #reunion of all beta indices
    beta_fd_indice_tot<- data.frame(average_beta_fd_indices, Mean_MPFD=Mean_MPFD,Sd_MPFD=SD_MPFD,
                                    Mean_MNND=Mean_MNND,Sd_MNND=SD_MNND,beta.hill[,-c(1,2)])
    
    
  } else{
    Beta_q0 <- get_beta_val( Data= lapply( list_beta_hill, function (x) x$beta_fd_q$q0),nb_rd=nb_rd,sites= sites_list,names="Hill_q0")
    Beta_q1 <- get_beta_val( Data= lapply(list_beta_hill, function (x) x$beta_fd_q$q1),nb_rd=nb_rd,sites= sites_list,names="Hill_q1")
    Beta_q2 <- get_beta_val( Data= lapply(list_beta_hill, function (x) x$beta_fd_q$q2),nb_rd=nb_rd,sites= sites_list,names="Hill_q2")
    
    beta.hill <- cbind(Beta_q0,Beta_q1[,c(3,4)],Beta_q2[,c(3,4)])
    
    #reunion of all beta indices
    beta_fd_indice_tot<-beta.hill[,-c(1,2)]#only Hill number computed in abundance
    
  }
  
  } else { # in case we do not calculate beta indices of functional diversity
    beta_fd_indice_tot <- NULL
  }
  #(exit: beta_fd_indice_tot)
  #############################################################################################################################
  
  
  if (do.hill.chao) {
  #PART 4: Hill number computation (alpha, beta, gamma) according to Chao et al. 2019
  cat("######## Hill numbers Chao ########","\n")
  
  List_100_Hill_func<-lapply(seq(1, length(q_hill)), function(i){
    DFtemp= do.call(rbind,lapply(seq(1, length(list_100_sp_dist)), function(x){
      FD_Decomp(data= t(asb_sp_w_beta), dij=as.matrix(list_100_sp_dist[[x]]), #data need to be the transposed matrix to get sites as col and species as rows + need to take the asb_sp_w_beta (without "total_site")
                tau=mean(list_100_sp_dist[[x]]), #tau suggested to be chosen as the "mean" of the "dist" object (Chao et al. 2019)
                q=q_hill[i], type="C")
    }))
  })
  names(List_100_Hill_func)<-paste0("q_hill_",q_hill)    
  
  #mean calculation
  Hill_func_Chao<-do.call(rbind, lapply(seq(1, length(List_100_Hill_func)), function(x){
    c(Mean=apply(List_100_Hill_func[[x]], 2, mean, na.rm=T),Sd=apply(List_100_Hill_func[[x]], 2, sd, na.rm=T))
  }))
  rownames(Hill_func_Chao)<-names(List_100_Hill_func)
  
  #(exit: Hill_func_Chao)
  ############################################################################################################################# 
  } else {
    Hill_func_Chao <- NULL
  }
  
  # FINAL RETURN
  return (list(Alpha=average_all_alpha_fd_indices,
               Beta=beta_fd_indice_tot,
               Gamma=average_all_alpha_fd_indices[nrow(average_all_alpha_fd_indices),],
               Hill_Func_Chao = Hill_func_Chao))
  
  
}# end of functional diversity pipeline  
##################################################################################################################################################

