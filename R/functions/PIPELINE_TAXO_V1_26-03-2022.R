####################################################################################################################################################
# Taxonomic pipeline to calculate phylogenetic Diversity gamma,alpha and beta indices and Hill numbers
#based on an eDNA table 

#Author: Romane Rozanski


#inputs:
#tab_eDNA: Data frame containing the number of copies/reads for each taxa detected with eDNA (filters in columns, taxa in rows)
#abundance =TRUE/FALSE: Choice if the indices are computed using species abundance (here the number of copies/reads from the eDNA table) or presence/absence

#returns: 
#a list divided in four lists containing the alpha-diversity indices(with hill), beta-diversity indices (with hill), 
#gamma-diversity indices and partitioned Hill numbers (from Chao et al. 2014) )



####################################################################################################################################################




################################### PIPELINE ###########################################################################

#' taxonomic_diversity
#' Calculates diversity indexes based on taxonomy
#' @param tab_eDNA a dataframe containing abundance of species by filter (site 
#' as rows, species as columns) and taxonomy. 
#' The format must be:
#' Class     Order     Family     Taxon          Filt1_Site1   Filt2_Site1   Filt3_Site2   ....
#' Actino    Clupei    Engrau     Engrau encra     NA             13200         58000
#' ...
#' If has.replicates : the name of the columns should be formatted like this : 
#' namefilt_namesite 
#' and the function will add together the abundances of filters for a same site.
#' If not has.replicates : the name of the columns is free.
#' @param abundance a boolean (default FALSE) : if T the analysis is based on abundance
#' data. 
#' If F we only calculate incidence (presence-absence) indices
#' @param has.replicates a boolean (default TRUE) : T if there are > 1 replicates
#' by site (in this case, the function will add abundances by site). 
#' If F each column is treated as a independant site. 
#'
#' @return a list divided in four lists containing the alpha-diversity indices 
#' (with hill), beta-diversity indices (with hill), #gamma-diversity indices and
#'  partitioned Hill numbers (from Chao et al. 2014))
#'
#' @author Romane Rozanski
#' @note Last update 2022/04/08 (Pierre Veron) : correction of a minor bug for the 
#' case abundance == T.
#' @note Update 2022/03/30 (Pierre Veron) : added the condition has.replicates
taxonomic_diversity<-function(tab_eDNA, abundance=FALSE, has.replicates = T){
  
  ### Part 1 : transformation of tab_EdDNAin community matrice 
  cat("########  Data preparation for indicators calculation  ########", "\n")
  #(entry: tab_eDNA)
  
  #delete the useless columns of classification
  rownames(tab_eDNA)<-tab_eDNA$Taxon
  tab_eDNA<-tab_eDNA[,-c(1:4)] #classification columns
  
  if (has.replicates) {
    #Step 1: grouping the columns (filters) with the same site
    split_col<-strsplit(colnames(tab_eDNA),"_")
    
    list_split_col_modif<-lapply(seq(1:length (split_col)), function(x){
      split_col[[x]]<-split_col[[x]][-1]})
    
    list_concat<-lapply(seq(1:length(list_split_col_modif)), function(x){
      str_c(list_split_col_modif[[x]], collapse="_")})
    
    names(tab_eDNA)<-list_concat
    
    #Step 2: pooling the same site together and adding the number of copies of each species
    site_list<- (unique(colnames(tab_eDNA))) 
    
    tab_eDNA_sites<- do.call(cbind,lapply(seq(1,length(site_list)), function(x){
      DFtempt <-  tab_eDNA[,which(colnames(tab_eDNA)==site_list[x])]
      if(length(which(colnames(tab_eDNA)==site_list[x]))>1){
        rowSums(DFtempt, na.rm=TRUE) # total sum of copies for species in the same site 
      } else {
        DFtempt
      } 
    }))  
    colnames(tab_eDNA_sites)<-site_list
  } else {
    tab_eDNA_sites <- tab_eDNA
  }
  #Step 3: transformation into matrice
  asb_sp_w<-t(tab_eDNA_sites) #transpose the matrice to get species in col and sites in rows
  #asb_sp_w[is.na(asb_sp_w)]<-0 #transformation "NA" in "0"
  
  #calculating the total number of copies by species for the "global site" to do the gamma diversity
  total_site<-colSums(asb_sp_w)
  asb_sp_w<-rbind(asb_sp_w,total_site)
  
  if (abundance==FALSE){
    q=c(0) #for hill number computation
    asb_sp_w<-replace(asb_sp_w,asb_sp_w[,]>0,1) #pres_abs modif
  } else {
    q=c(0,1,2)
  }
  
  #asb_sp_w<-subset(asb_sp_w,rowSums(asb_sp_w)>=3) #
  #(exit: asb_sp_w)
  ###########################################################################################################################################
  
  
  ###Part 2: Taxonomic alpha indices calculation
  
  
  #Step 1: alpha indices
  cat("########  alpha taxonomic indices ########", "\n")
  
  #species richness calculation
  sp_rich<-rowSums(asb_sp_w)
  
  #indices in abundance : WARNING 3 species/site minimum to calculate the indices
  if (abundance==TRUE){
    shannon_ind<-vegan::diversity(asb_sp_w, index="shannon", base=exp(1)) #shannon
    simpson_ind<-vegan::diversity(asb_sp_w, index="simpson") #simpson indice
    pielou_ind <-shannon_ind/log(specnumber(asb_sp_w))  #Pielou indice
    alpha_ind<-data.frame(sp_richness=sp_rich,Shannon=shannon_ind, Simpson=simpson_ind, Pielou=pielou_ind)
  } else{
    alpha_ind<-data.frame(sp_richness=sp_rich)
  }
  
  
  
  #Step 2: Alpha Hill number 
  cat("########  alpha Hill number  ########", "\n")
  
  #hill alpha partition
  Hill_q0_a_part <- hillR::hill_taxa_parti(asb_sp_w, q=q[1],rel_then_pool=FALSE)
  
  #hill alpha sites
  Hill_a_q0 <- hillR::hill_taxa(asb_sp_w, q=q[1],MARGIN = 1, base = exp(1))
  
  if (abundance ==TRUE){
    #partition
    Hill_q1_a_part <- hillR::hill_taxa_parti(asb_sp_w,q=q[2],rel_then_pool=FALSE)
    Hill_q2_a_part<- hillR::hill_taxa_parti(asb_sp_w, q=q[3],rel_then_pool=FALSE)
    alpha_hill_part<-rbind(Hill_q0_a_part,Hill_q1_a_part,Hill_q2_a_part)
    
    #site
    Hill_a_q1 <- hillR::hill_taxa(asb_sp_w, q=q[2],MARGIN = 1, base = exp(1))
    Hill_a_q2 <- hillR::hill_taxa(asb_sp_w, q=q[3],MARGIN = 1, base = exp(1))
    alpha_hill<-cbind(Hill_a_q0,Hill_a_q1,Hill_a_q2)
    
  } else {
    alpha_hill_part<- Hill_q0_a_part #partition
    alpha_hill<-Hill_a_q0            #site : equal to species richness
  }
  
  #merging every alpha indices
  alpha_ind_tot<-cbind(alpha_ind,alpha_hill)
  
  #(exit: alpha_ind_tot)
  ###########################################################################################################################################
  
  
  ### Part 3: beta diversity indices
  cat("########  taxonomic beta diversity indices  ########", "\n")
  
  #modif asb_sp_w 
  asb_sp_w_beta<-subset(asb_sp_w,row.names(asb_sp_w) != "total_site") #delete this row because only used for gamma div
  
  
  #Step1: Classic beta indices : jaccard dissimilarities only in presence/absence (bjac=bjtu+bjne, see Baselga 2012)
  if (abundance==FALSE){
    #conversion asb_sp_w_beta in pres/abs
    #asb_sp_w_beta_pa<-replace(asb_sp_w_beta,asb_sp_w_beta[,]>2,1)
    
    #function to get bjac, bjtu, bjne dissimilarity indices
    beta<-betapart::beta.pair(as.data.frame(asb_sp_w_beta), index.family = "jaccard") #only pres/abs for these indice
    
    #convert the dist object into  a data frame
    list_beta_mat<-lapply(seq(1:length(beta)), function(x){
      beta_temp<-as.matrix(beta[[x]])
      beta_temp[upper.tri(beta_temp,diag=T)] <- NA
      #beta_temp<-beta_temp[order(rownames(beta_temp)),]
      beta_temp2<-na.omit(melt(beta_temp))
    })
    
    names(list_beta_mat)<-names(beta)
    
    
    #merging the 3 indices together
    beta_ind<-data.frame(site1=list_beta_mat[[1]][,1], site2=list_beta_mat[[1]][,2],
                         bjtu=list_beta_mat[[1]][,3], bjne=list_beta_mat[[2]][,3], bjac=list_beta_mat[[3]][,3])
    
    
    
    rownames(beta_ind)<-str_c(as.character(beta_ind$site2),paste="-",as.character(beta_ind$site1))
    beta_ind[,1:2]<-NULL
    
  } # end of "if abundance==false"
  
  
  
  #Step 2: Beta Hill number 
  cat("########  taxonomic beta Hill number computation  ########", "\n")
  
  #function  to compute beta hill number
  get_q_taxo<-function(q, comm=asb_sp_w_beta){
    b_hill <- as.data.frame(hillR::hill_taxa_parti_pairwise(comm,q =q,pairs = "unique"))
  }
  
  #application of the previous function: "get_q_taxo"
  list_q<-lapply(q, function(x){
    get_q_taxo(q=x,comm=asb_sp_w_beta)
  })
  
  for (x in (1:length(list_q))){
    colnames(list_q[[x]])[4:ncol(list_q[[x]])]<-paste0(colnames(list_q[[x]])[4:ncol(list_q[[x]])], "_q", list_q[[x]][1,1])
  }
  
  if (abundance==TRUE){
    Hill_beta <- cbind(list_q[[1]][,-1], list_q[[2]][,-c(1,2,3)], list_q[[3]][,-c(1,2,3)])
  } else{
    Hill_beta<-list_q[[1]][,-1]
  }
  
  
  rownames(Hill_beta) <- str_c(Hill_beta$site1, paste="-", Hill_beta$site2) #warning: need to be the same order as beta_ind to be merged
  Hill_beta[,1:2]<-NULL
  
  #merging all the beta indices
  if (abundance) {
    beta_ind_tot <- Hill_beta
  } else {
    beta_ind_tot<-merge(beta_ind, Hill_beta, by=0, all=T)
  }
  #(exit: beta_ind_tot)
  ###########################################################################################################################################
  
  
  #final return
  return(list(Alpha=alpha_ind_tot[-nrow(alpha_ind_tot),],
              Beta=beta_ind_tot,
              Gamma=alpha_ind_tot[nrow(alpha_ind_tot),],
              Hill_partition=alpha_hill_part))
  
}##end of taxonomic diversity function

#########################################################################################



