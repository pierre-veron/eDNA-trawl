## Script description ----------------------------------------------------------
# This script implements Anne Chao's function for calculation of functional Hill 
#  number for alpha, beta and gamma diversity.
# It tests this function on an example taken from the mFD package. 
#
# Author : Pierre Veron (pierre.veron.2017@polytechnique.org), M2 student under
#  Camille Albouy's supervision at IFREMER.
#
# Created : 2022-03-09
#
# Remark : the function FD_Beta was originally designed to calculate 
#  dissimilarity indices as defined in Chao et. al 2019, (DOI:10.1002/ECM.1343)
#  but I changed it to obtain the desired Hill numbers that were already
#  calculated by the function. Only a few line were changed and all changes are 
#  labelled with a comment #PV


## Load Anne Chao's function ---------------------------------------------------

# Source : https://github.com/AnneChao/FunD/blob/master/FunD_Rcode.txt

#' FD_Decomp(data, dij, tau, q, type) is a function of obtaining functional dissimilarity measure of order q.
#' @param data a SxN matrix of species sample frequencies with S species, N sites.
#' @param dij a matrix of species pairwise distances.
#' @param tau a numeric of a specified level of threshold distinctiveness.
#' @param q a numeric of a specified diversity order.
#' @param type a character to choose index: "C" for 1-CqN ; "U" for 1-UqN.
#' @return PV: a vector containing gamma, alpha and beta Hill numbers

FD_Decomp <- function(data, dij, tau, q, type){
  N <- ncol(data)
  if (class(dij)=="matrix") {
    dij <- lapply(seq_len(N), function(i) dij)
  }
  dij <- lapply(seq_along(dij), function(i) {
    d <- as.matrix(dij[[i]])
    d[which(d>tau,arr.ind = T)] <- tau
    d
  })
  aik <- sapply(seq_along(dij), function(i) {
    (1-dij[[i]]/tau) %*% data[, i]  #matrices multiplication
  })
  aiplus <- apply(aik, 1, sum)
  vi <- apply(data, 1, sum)/aiplus  
  alpha_v <- rep(vi, N)
  nplus <- sum(data)
  aik <- as.vector(aik)
  alpha_v <- alpha_v[aik!=0]
  aik <- aik[aik!=0]
  if(q==1){
    gamma=exp(sum(-vi*aiplus/nplus*log(aiplus/nplus)))
    alpha=1/N*exp(sum(-alpha_v*aik/nplus*log(aik/nplus)))
    beta = max(1, gamma/alpha)
    #PV : I do not need the following line. 
    #out <- log(beta)/log(N)
  }else{
    gamma=(sum(vi*(aiplus/nplus)^q))^(1 / (1-q))
    alpha=1/N*(sum(alpha_v*(aik/nplus)^q))^(1 / (1-q))
    beta = max(1, gamma/alpha)
    #PV : I do not need the following two lines 
    #if(type == "C") out <- (1-(beta)^(1-q))/(1-(N)^(1-q))
    #if(type == "U") out <- (1-(beta)^(q-1))/(1-(N)^(q-1))
  }
  out <- c("gamma" = gamma, "alpha" = alpha, "beta" = beta) #PV : this line added to get Hill numbers
  out
}
