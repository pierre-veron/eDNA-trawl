
#function to use in  the functional pipeline to ge the beta hill number 


get_beta_val <- function(Data=Beta1,nb_rd=nb_rd,sites= sites_list,names="jac_diss"){
  
  DFtemp <- array(data=NA,c(nrow(as.matrix(Data[[1]])),ncol(as.matrix(Data[[1]])),nb_rd)); i=1
  for ( i in 1:nb_rd) DFtemp[,,i] <- as.matrix(Data[[i]])
  
  res=apply(DFtemp,c(1,2),mean) 
  res_SD=apply(DFtemp,c(1,2),sd) 
  
  rownames(res) <- colnames(res) <- sites[-length(sites)]
  res[upper.tri(res,diag=T)] <- NA
  res_SD[upper.tri(res_SD,diag=T)] <- NA
  res_tot <- na.omit(cbind(expand.grid(dimnames(res)), Mean = as.vector(res),SD=as.vector(res_SD)))
  colnames(res_tot)[c(3,4)] <- paste(names,colnames(res_tot)[c(3,4)],sep="_")
  res_tot
}