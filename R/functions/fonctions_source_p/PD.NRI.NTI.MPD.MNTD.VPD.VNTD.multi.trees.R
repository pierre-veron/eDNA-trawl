# Compute SR, MPD, MNTD, VPD and  VNTD metrics (See Tucker et al. 2017 for the definition of the metrics)
# Author David Eme, 29/10/2019.

# tree =  a phylogenetic tree, a phylo object,
# matrix =  a data.frame with only one row (one community),
# return a data frame with SR (species richness), MPD (Mean pairwise distance), MNTD (Mean nearest taxonomic distance), 
# VPD (Variance of the pairwise distance) and VNTD (Variance of the nearest taxonomic distance).

require(ape)
require(PhyloMeasures)

PD.NRI.NTI.MPD.MNTD.VPD.VNTD = function(tree = NULL, matrix = NULL, null.model = NULL, reps = 999){

  SR = sum(matrix)
  NTI<-mntd.query(tree=tree, matrix=matrix , standardize = TRUE, null.model=null.model, reps=reps)
  NRI<-mpd.query(tree=tree, matrix=matrix, standardize = TRUE, null.model=null.model, reps=reps)
  PD<-pd.query(tree=tree, matrix= matrix, null.model=null.model, reps=reps)
  
  sampleTaxa = colnames(matrix)[which(matrix > 0)]
  
  SpToRemove = setdiff(tree$tip.label, sampleTaxa) # species list to remove from the tree.
  tmp.tree=drop.tip(tree, SpToRemove) # prune tree
  
  if(SR > 1) {
    com.distD = cophenetic(tmp.tree)
    com.dist=as.dist(com.distD)
    MPD=mean(com.dist) ### MPD=AvTD=Delta+. ="Dispersion index"
    VPD=((sum(com.dist^2))/((SR*(SR-1))/2))-(mean(com.dist)^2) ### VPD=Var.MPD (Tucker et al. 2017) give the same results that the Lambda (VarTD: Clarke & Warwick 2001) computed by taxondive function from vegan. ="Regularity index"
    diag(com.distD)=NA
    MNTD=mean(apply(com.distD, MARGIN=1, min, na.rm=T)) ### Mean nearest taxon distance (mntd) See Swenson 2014. ="Dispersion index"
    VNTD=sum((apply(com.distD, MARGIN=1, min, na.rm=T)-mean(apply(com.distD, MARGIN=1, min, na.rm=T)))^2)/SR ### Variance of MNTD (Tucker et al. 2017). ="Regularity index"
  } else {
    MPD = NA
    VPD = NA
    MNTD = NA
    VNTD = NA
  }
  data.frame(SR = SR, PD=PD, NTI=NTI,NRI=NRI, MPD = MPD, MNTD = MNTD, VPD = VPD, VNTD = VNTD)
}

# Compute the function on multiple trees and large community data matrix, it runs in parallel.

# trees = a "multiPhylo" object,
# Com.Mat =  the community data matrix with the species names matching perfectly
# the tip label of the trees.
# save.path =  name of the path storing the .rds files (e.g. "Results/Indices_100trees/")
# nthreads  = number of threads that shoukd be used to run the analyses in parallel.

require(parallel)

PD.NRI.NTI.MPD.MNTD.VPD.VNTD.multi.trees = function(trees = NULL, Com.Mat = NULL, save.path = NULL, nthreads = NULL, null.model = NULL, reps = 999){ 
  
  if(class(trees) == "phylo" ){
    trees = list(trees)
  }

  lapply(1:length(trees), function(y){
    
    Index = mclapply(1:nrow(Com.Mat), mc.cores = 1, function(x){
      #Index = lapply(1:nrow(Com.Mat),  function(x){
      #cat("x=",x,"\n")
      PD.NRI.NTI.MPD.MNTD.VPD.VNTD (tree = trees[[y]], matrix = Com.Mat[x,], null.model = null.model, reps = reps)
    })
    
    Index_results <- do.call(rbind,Index)
    row.names(Index_results) = row.names(Com.Mat)
    
    saveRDS(Index_results, paste0(save.path, "Index_results_Tree_",y,".rds"))
  })
}

