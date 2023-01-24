###############
## 2 functions used in the phylo PIPELINE 
##############


# function 1 : beta Hill number computation: calculate beta hill number for the 3 orders of q :
#q=0 species richness/ q=1 shannon entropy / q=2 inverse Simpson
#give the mean and sd of beta hill number on each pair of sites

get_q<-function(q,nb_rd,trees=data_acti_chond_trees, comm=asb_sp_w_beta){
  b_hill <- lapply(seq(1:length(trees)), function(i){
    as.data.frame(hillR::hill_phylo_parti_pairwise(comm, trees[[i]], q=q, pairs="unique",rel_then_pool=FALSE, .progress = F))
  }) 
  
  n <- ((nrow(comm)*nrow(comm))-nrow(comm))/2
  res <- array(data=NA,c(n,5,nb_rd)); i=1
  for ( i in 1:nb_rd) res[,,i] <- as.matrix(b_hill[[i]][,4:8])
  Mean_q<-apply(res, c(1,2),mean) #mean calculation
  data<-cbind(b_hill[[1]][,2:3],Mean_q)
  colnames(data)<-colnames(b_hill[[1]])[-1]
  colnames(data)[3:7]<-paste0("Mean_",colnames(data)[3:7], "_q", q)
  
  Sd_q<-apply(res, c(1,2),sd) #sd calculation
  data_sd<-cbind(b_hill[[1]][,2:3],Sd_q)
  colnames(data_sd)<-colnames(b_hill[[1]])[-1]
  colnames(data_sd)[3:7]<-paste0("sd_",colnames(data_sd)[3:7], "_q", q)
  
  data_tot<-cbind(data,data_sd[,3:7])
  data_tot
}



#function 2 : Calculation of mean and sd of beta dissimilarity indices (Jac, jtu, jne) for each pair of sites
get_resbet <- function(x=asb_sp_w_beta_pa,data=list_bjac,nb_rd=nb_rd,type="mean") {
  res <- array(data=NA,c(nrow(x),nrow(x),nb_rd)); i=1
  for ( i in 1:nb_rd) res[,,i] <-  as.matrix(data[[i]])
  
  if (type=="mean"){
    Mean<-apply(res, c(1,2),mean) 
    rownames(Mean) <- colnames(Mean) <- rownames(x)
    Mean[upper.tri(Mean,diag=T)] <- NA
    return(na.omit(melt(Mean)))
  }
  
  if(type=="sd") {
    SD <- apply(res, c(1,2),sd)
    rownames(SD) <- colnames(SD) <- rownames(x)
    SD[upper.tri(SD,diag=T)] <- NA
    return(na.omit(melt(SD)))
  }
  
}


