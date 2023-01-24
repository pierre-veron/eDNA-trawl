####################################################################################################################################################
# Phylogenetic pipeline to calculate phylogenetic Diversity gamma,alpha and beta indices and phylogenetic Hill numbers
#based on a phylogenetic tree and a eDNA table 

#Authors: Romane Rozanski, Camille Albouy, David Eme


#inputs:
#tab_eDNA: Data frame containing the number of copies/reads for each taxa detected with eDNA (filters in columns, taxa in rows)
#trees: list of phylogenetic trees (class multiphylo)
#Classif_sp: a data frame containg the classification of a set of species belonging to a genus/family detected by eDNA (at least Family, genus and species ranks)
#savepathindic= name of the path storing the .rds files of the alpgha indices function (e.g. "PHYLO/Indices_100trees/")
#abundance =TRUE/FALSE: Choice if the indices are computed using species abundance (here the number of copies/reads from the eDNA table) or presence/absence
#nb_rd:number of randomisation depending on the number of trees contained in trees object


#returns: 
#a list divided in four lists containing the alpha-diversity indices, beta-diversity indices, gamma-diversity indices and partitioned Hill numbers (from Chao et al. 2014) )


#Last update: 28/03/2022 (by Romane Rozanski)
####################################################################################################################################################




#source functions
source(here("R", "functions", "fonctions_source_p", "PD.NRI.NTI.MPD.MNTD.VPD.VNTD.multi.trees.R"))
source(here("R", "functions", "fonctions_source_p","Summary.Stat.R"))
source(here("R", "functions", "fonctions_source_p","get_q_&_get_resbet.R"))
source(here("R", "functions", "fonctions_source_p","Functional_&_Phylogenetic.Beta.MPFD.MNND.R"))




################################### PIPELINE ###########################################################################

#' phylo_diversity
#'
#' Phylogenetic pipeline to calculate phylogenetic Diversity gamma,alpha and beta indices and phylogenetic Hill numbers
#' based on a phylogenetic tree and a eDNA table 
#' @param trees list of phylogenetic trees (class multiphylo)
#' @param tab_eDNA  dataframe containing abundance of species by filter (site 
#' as rows, species as columns) and taxonomy. 
#' The format must be:
#' Class     Order     Family     Taxon          Filt1_Site1   Filt2_Site1   Filt3_Site2   ....
#' Actino    Clupei    Engrau     Engrau encra     NA             13200         58000
#' ...
#' If has.replicates : the name of the columns shound be formated like this : 
#' namefilt_namesite 
#' and the function will add together the abundances of filters for a same site.
#' If not has.replicates : the name of the columns is free.
#' @param Classif a data frame containg the classification of a set of species 
#' belonging to a genus/family detected by eDNA (at least Family, genus and species ranks)
#' @param savepathindic name of the path storing the .rds files of the alpha 
#' indices function (e.g. "PHYLO/Indices_100trees/"). 
#' WARNING : must be a directory, ended by a "/" (at least on Windows). 
#' @param abundance TRUE/FALSE: Choice if the indices are computed using species
#' abundance (here the number of copies/reads from the eDNA table) or presence/absence
#' @param has.replicates a boolean (default TRUE) : T if there are > 1 replicates
#' @param do.beta a boolean (default T) : calculates beta indices 
#' by site (in this case, the function will add abundances by site). 
#' If F each column is treated as a independant site. 
#' @return #a list divided in four lists containing the alpha-diversity indices,
#' beta-diversity indices, gamma-diversity indices and partitioned Hill numbers 
#' (from Chao et al. 2014) )
#' @export
#'
#' @author Romane Rozanski, Camille Albouy, David Eme
#' @note Updated : 2022/03/30 by Pierre Veron. Corrections : 
#'       - the condition has.replicates and save_fct_space
#'       - add the possibility that none of taxa is detected at genus/family level
#'       - minor bug : Classif_sp was an external variable, set it to Classif
#'       - minor bug : the replacement Genus species to Genus_species did not take 
#'          into account genus or family random sampling 
#'       - remove the default of "savepathindic" because it must be specified
#'       - remove the parameter nb_rd which is always equal to length(trees)
#'  @note Minor update 2022/06/24, add the parameter do.beta 
phylo_diversity<-function(trees = trees_acti_chond, tab_eDNA,Classif, savepathindic,  abundance=FALSE,
                          has.replicates = T, do.beta = T){
  nb_rd <- length(trees)
  Classif_sp <- Classif # PV
  ############################################################################################################################################ 
  
  ###PART 1 : transformation of tab_eDNA in community matrice whose species are matching with the species in the tree
  cat("########  Data preparation for indicators calculation  ########", "\n")
  #(entry: tab_eDNA, trees, Classif_sp)
  
  # Step 1: Choosing a random representative species for taxa detected at the genus/family level in tab_eDNA 
  #because phylo tree only contains species so to match with the taxa from tab_eDNA, all taxa need to be at the species level
  
  #to  match with the species name in Classif_sp and phylo trees
  
  
  ##searching for fish detected at GENUS or FAMILY level by eDNA
  genus_vect<-c()
  fam_vect<-c()
  
  for (i in 1:nrow(tab_eDNA)){
    if (tab_eDNA$Taxon[[i]] %in% Classif_sp$Genus){ #genus
      genus_vect<-c(genus_vect,tab_eDNA$Taxon[[i]])
    }
    if (tab_eDNA$Taxon[[i]] %in% Classif_sp$Family){  #family
      fam_vect<-c(fam_vect,tab_eDNA$Taxon[[i]])
    }
  } #end of "for"
  
  
  #searching all the species representing the genus/family in the Classif and replacing the taxa by the species randomly chosen
  #based on the hypothesis that a genus or a family represents a monophyletic group so the choice of a species over another doesn't change anything
  
  #genus
  if (length(genus_vect) > 0) {
    for (x in (1:length(genus_vect))){
      vec_temp=Classif_sp[which(genus_vect[x]==Classif_sp$Genus),"Species"]#species belonging to the same genus
      vec_temp<-subset(vec_temp, ! (vec_temp %in% tab_eDNA$Taxon)) # WARNING: some taxa identified at both genus level and species level: need to have 2 different species
      sample_sp=sample(vec_temp,size=1)  
      tab_eDNA$Taxon[which(tab_eDNA$Taxon==genus_vect[x])] <-sample_sp
      
    }
  }
  
  #family (no possibility of taxa identified at both family level and species)
  if (length(fam_vect) > 0) {
    for (x in (1:length(fam_vect))){
      sample_sp=sample(Classif_sp[which(fam_vect[x]==Classif_sp$Family),"Species"], size=1)#random choice of the species among the ones belonging to the same genus
      tab_eDNA$Taxon[which(tab_eDNA$Taxon==fam_vect[x])] <-sample_sp
    }
  }
  tab_eDNA$Taxon<-gsub(" ","_",tab_eDNA$Taxon) 
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
  
  # Step 3: Comparison species of fish data and the tree containing actino and chondri
  SpToRemove <- setdiff(trees[[1]]$tip.label, rownames(tab_nb_sites)) # species list to remove from the trees because not in the eDNA table.
  
  data_acti_chond_trees <- lapply(seq(1,length(trees)), function(i){
    drop.tip(trees[[i]], SpToRemove) #removing of the species of the trees that are not in the datas tab_eDNA
  })
  
  
  #Step 4: transformation into matrice
  asb_sp_w<-t(tab_nb_sites) #transpose the matrice to get species in col and sites in rows
  asb_sp_w[is.na(asb_sp_w)]<-0 #transformation "NA" in "0"
  
  #calculating the total number of copies by species for the "global site" to do the gamma diversity
  total_site<-colSums(asb_sp_w)
  asb_sp_w<-rbind(asb_sp_w,total_site)
  
  
  #asb_sp_w_pa<-subset(asb_sp_w_pa,rowSums(asb_sp_w_pa)>=3) #3 species/site minimum to calculate the indices
  
  
  
  #(exit: asb_sp_w_pa, data_acti_chond_trees)
  ############################################################################################################################################ 
  
  
  ### PART 2 : calculation phylogenetic alpha indices
  #(entry:asb_sp_w_pa, data_acti_chond_trees)
  cat("######## phylogenetic alpha-diversity indices ########", "\n")
  
  
  ### Step 1: alpha, gamma, phylo indices: WARNING only in presence/absence 
  if (abundance==FALSE){
    asb_sp_w<-replace(asb_sp_w,asb_sp_w[,]>0,1) #pres_abs modif # PV : corrected > 0
    
    indices_alpha_gamma_phylo_acti_chond <- PD.NRI.NTI.MPD.MNTD.VPD.VNTD.multi.trees(trees = data_acti_chond_trees, Com.Mat=as.data.frame(asb_sp_w),
                                                                                     save.path= savepathindic,nthreads=3,
                                                                                     null.model="uniform",reps = 1)
    
    #calculation mean & sd alpha, gamma
    Mean_alpha_phylo <- Summary.Stat(path=savepathindic,prefix=NULL)
    
  } 
  ###Step 2: Alpha Hill number computation
  cat("########  alpha Hill number  ########", "\n")
  
  if (abundance== TRUE){
    q=c(0,1,2) 
    
    #partition of hill number gamma, alpha, beta
    list_alpha_hill_part <- lapply(seq(1:length(data_acti_chond_trees)), function(i){
      q0 <- hillR::hill_phylo_parti(asb_sp_w, data_acti_chond_trees[[i]], q=q[1], phy_abund=NULL,rel_then_pool=FALSE)
      q1 <- hillR::hill_phylo_parti(asb_sp_w, data_acti_chond_trees[[i]], q=q[2], phy_abund=NULL,rel_then_pool=FALSE)
      q2 <- hillR::hill_phylo_parti(asb_sp_w, data_acti_chond_trees[[i]], q=q[3], phy_abund=NULL,rel_then_pool=FALSE)
      rbind(q0,q1,q2)
    })
    
    #alpha hill on each site
    alpha_hill <- lapply(seq(1:length(data_acti_chond_trees)), function(i){
      q0 <- hillR::hill_phylo(asb_sp_w, data_acti_chond_trees[[i]], q=q[1],rel_then_pool=FALSE)
      q1 <- hillR::hill_phylo(asb_sp_w, data_acti_chond_trees[[i]], q=q[2],rel_then_pool=FALSE)
      q2 <- hillR::hill_phylo(asb_sp_w, data_acti_chond_trees[[i]], q=q[3],rel_then_pool=FALSE)
      cbind(q0,q1,q2)
    })
    
  } else{
    q=c(0) #q=0 phylogenetic richness
    
    
    #partition of hill number gamma, alpha, beta
    list_alpha_hill_part <- lapply(seq(1:length(data_acti_chond_trees)), function(i){
      q0 <- hillR::hill_phylo_parti(asb_sp_w, data_acti_chond_trees[[i]], q=q[1], phy_abund=NULL,rel_then_pool=FALSE)
    })
    
    #alpha hill on each site
    alpha_hill <- lapply(seq(1:length(data_acti_chond_trees)), function(i){
      q0 <- data.frame(q0=hillR::hill_phylo(asb_sp_w, data_acti_chond_trees[[i]], q=q,rel_then_pool=FALSE))
    })
  }
  
  
  #mean  and sd hill alpha partition
  a_hill_meanpart <- do.call(rbind, lapply(q, function(x){
    DFtemp = do.call(rbind, lapply(seq(1, length(list_alpha_hill_part)), function(i){
      list_alpha_hill_part[[i]][list_alpha_hill_part[[i]][,1]== x,]
    }))
    c(Mean=apply(DFtemp, 2, mean, na.rm=T),Sd=apply(DFtemp, 2, sd, na.rm=T))
  }))
  
  #mean and sd alpha Hill on each site 
  res <- array(data=NA,c(nrow(alpha_hill[[1]]),length(q),nb_rd)); i=1
  for ( i in 1:nb_rd) res[,,i] <-  as.matrix(alpha_hill[[i]])
  res_alpha_hill <- cbind(Mean = apply(res,c(1,2),mean),SD = apply(res,c(1,2),sd))
  colnames(res_alpha_hill) <- c(paste0("q",q,"_mean"),paste0("q",q,"_sd"))
  rownames(res_alpha_hill) <- rownames(alpha_hill[[1]])    
  
  
  #reunion of all alpha indices
  if (abundance==FALSE){
    ind_alpha_phylo_tot_final<-cbind(Mean_alpha_phylo[-nrow(Mean_alpha_phylo),],res_alpha_hill[-nrow(res_alpha_hill),])
    ind_gamma_phylo_tot_final<-Mean_alpha_phylo[nrow(Mean_alpha_phylo),]
    
  }else{
    ind_alpha_phylo_tot_final<-res_alpha_hill
    ind_gamma_phylo_tot_final<-res_alpha_hill[nrow(res_alpha_hill),]
  }
  
  #(exit: ind_alpha_phylo_tot_final, ind_gamma_phylo_tot_final)
  ############################################################################################################################################ 
  
  
  ### PART 3 : calculation phylogenetic beta diversity indices
  if (do.beta) {
    cat("########  beta phylo indices  ########", "\n")
    
    #modif asb_sp_w 
    asb_sp_w_beta<-subset(asb_sp_w,row.names(asb_sp_w) != "total_site") #delete this row because only used for gamma div
    
    #Step 1 : dissimilarity indices : only in presence/absence
    if (abundance==FALSE){
      
      list_beta_indices<-mclapply(seq(1,length(data_acti_chond_trees)), mc.cores=1,function(i){
        betapart::phylo.beta.pair(as.data.frame(asb_sp_w_beta), data_acti_chond_trees[[i]], index.family="jaccard")
      })
      
      
      #split the big list of 100 lists in 3 different list for each indice:
      list_bjac<-list() ; list_bjtu<-list() ; list_bjne<-list()
      for (i in (1: length(list_beta_indices))){
        list_bjac[[i]]<-list_beta_indices[[i]][["phylo.beta.jac"]]   #jaccard
        list_bjtu[[i]]<-list_beta_indices[[i]][["phylo.beta.jtu"]]  # jaccard turnover
        list_bjne[[i]]<-list_beta_indices[[i]][["phylo.beta.jne"]]  # jaccard nestedness
      }
      
      #mean and sd on each list 
      #bjac : jaccard
      Mean_bjac <- get_resbet(x=asb_sp_w_beta,data=list_bjac,nb_rd=nb_rd,type="mean")
      SD_bjac <- get_resbet(x=asb_sp_w_beta,data=list_bjac,nb_rd=nb_rd,type="sd")
      
      #bjtu : turnover
      Mean_bjtu <- get_resbet(x=asb_sp_w_beta,data=list_bjtu,nb_rd=nb_rd,type="mean")
      SD_bjtu <- get_resbet(x=asb_sp_w_beta,data=list_bjtu,nb_rd=nb_rd,type="sd")
      
      #bjne : nestedness
      Mean_bjne <- get_resbet(x=asb_sp_w_beta,data=list_bjne,nb_rd=nb_rd,type="mean")
      SD_bjne <- get_resbet(x=asb_sp_w_beta,data=list_bjne,nb_rd=nb_rd,type="sd")
      
      
      
      #Step 2: others beta  diversity indices in presence/absence
      cat("########  others beta indices MPD, MNND  ########", "\n")
      
      #extraction of a distance matrix between species from the trees
      list_sp_dist <-lapply(seq(1, length(data_acti_chond_trees)), function(x){
        as.dist(cophenetic.phylo(data_acti_chond_trees[[x]]))})
      
      #Mean Pairwise Distance (MPD)
      list_beta_MPFD<- mclapply(seq(1, length(list_sp_dist)),mc.cores=1, function(x){
        comdist.Swenson(asb_sp_w_beta, list_sp_dist[[x]])})
      
      res <- array(data=NA,c(nrow(asb_sp_w_beta),nrow(asb_sp_w_beta),nb_rd)); i=1
      for ( i in 1:nb_rd) res[,,i] <-  as.matrix(list_beta_MPFD[[i]])
      
      Mean_MPFD <-apply(res, c(1,2),mean) 
      rownames(Mean_MPFD) <- colnames(Mean_MPFD) <- rownames(asb_sp_w_beta)
      Mean_MPFD[upper.tri(Mean_MPFD,diag=T)] <- NA
      Mean_MPFD <- na.omit(melt(Mean_MPFD))
      
      SD_MPFD <- apply(res, c(1,2),sd)
      rownames(SD_MPFD) <- colnames(SD_MPFD) <- rownames(asb_sp_w_beta)
      SD_MPFD[upper.tri(SD_MPFD,diag=T)] <- NA
      SD_MPFD <- na.omit(melt(SD_MPFD))
      
      
      #Mean Nearest Neighbour Distance (MNND)
      list_beta_MNND<-mclapply(seq(1, length(list_sp_dist)),mc.cores=1,function(x){
        comdistnt.Swenson(asb_sp_w_beta, list_sp_dist[[x]])})
      
      res <- array(data=NA,c(nrow(asb_sp_w_beta),nrow(asb_sp_w_beta),nb_rd)); i=1
      for ( i in 1:nb_rd) res[,,i] <-  as.matrix(list_beta_MNND[[i]])
      MNND <- list(Mean =apply(res, c(1,2),mean), SD=apply(res, c(1,2),sd))
      
      Mean_MNND <-apply(res, c(1,2),mean) 
      rownames(Mean_MNND) <- colnames(Mean_MNND) <- rownames(asb_sp_w_beta)
      Mean_MNND[upper.tri(Mean_MNND,diag=T)] <- NA
      Mean_MNND <- na.omit(melt(Mean_MNND))
      
      SD_MNND <- apply(res, c(1,2),sd)
      rownames(SD_MNND) <- colnames(SD_MNND) <- rownames(asb_sp_w_beta)
      SD_MNND[upper.tri(SD_MNND,diag=T)] <- NA
      SD_MNND <- na.omit(melt(SD_MNND))
      
      #reunion of all the computed indices
      beta_indice_phylo<- data.frame(Mean_bjac=Mean_bjac$value,Sd_bjac=SD_bjac$value,
                                     Mean_bjtu=Mean_bjtu$value,Sd_bjtu=SD_bjtu$value,
                                     Mean_bjne=Mean_bjne$value,SD_bjne=SD_bjne$value, 
                                     Mean_MPFD=Mean_MPFD$value,Sd_MPFD=SD_MPFD$value,
                                     Mean_MNND=Mean_MNND$value,Sd_MNND=SD_MNND$value)
      
      rownames(beta_indice_phylo)<-str_c(Mean_bjac$Var1,paste="-",Mean_bjac$Var2) 
      
      
    }#end of "if abundance==FALSE"
    
    
    
    #Step 3: beta Hill number computation 
    cat("########  beta Hill number computation  ########", "\n")
    
    
    list_q<-lapply(q, function(x){  # q value depending on "abundance=T/F" in the entry of the function
      get_q(q=x,nb_rd,trees=data_acti_chond_trees,comm=asb_sp_w_beta)
    })
    
    
    if (abundance==TRUE){
      Hill_beta <- cbind(list_q[[1]], list_q[[2]][,-c(1,2)], list_q[[3]][,-c(1,2)]) #q=0,1,2
    } else{
      Hill_beta<-list_q[[1]] #q=0 only
    }
    
    
    rownames(Hill_beta) <- str_c(Hill_beta$site2, paste="-", Hill_beta$site1)  #to get the same order as beta_indices_phylo
    Hill_beta[,1:2]<-NULL
    
    if (abundance == FALSE){
      beta_indice_phylo_tot_final<-merge(beta_indice_phylo, Hill_beta, by="row.names", all=T)
    }else{
      beta_indice_phylo_tot_final<- Hill_beta
    }
  } else { # does not calculate beta indices
    beta_indice_phylo_tot_final <- NULL
  }
  
  #exit: beta_indice_phylo_tot_final
  ############################################################################################################################################ 
  
  #final return:
  
  return(list(Alpha=ind_alpha_phylo_tot_final,
              Beta=beta_indice_phylo_tot_final,
              Gamma=ind_gamma_phylo_tot_final,
              Hill_Partition=a_hill_meanpart))
  
  
  
}# end of pipeline phylo
