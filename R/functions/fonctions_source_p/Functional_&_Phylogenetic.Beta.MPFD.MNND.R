

### Function name=comdist.Swenson (= Dpw, MPFD).
### This function computes the mean pairwise distance among species between two communities. This metric is also called Dpw or comdist. This function use a fast code based on apply function developed by Nathan Swenson in Swenson 2014 "Functional and Phylogenetic Ecology in R" Springer p99-100
### This function can also be used to compute the Mean pairwise Functional Distance (MPFD) based on any functional distance matrix (e.g. Gower distance matrix).
### Input:
# comm is a community data matrix (all the communities have to contain at least 1 species).
# dis should be a matrix containing the distances among species (e.g. coming from cophenetic() function for instance). If the dis object is a "dist" object
# the function atomatically convert it to a "matrix" object

## The conspecific distances among communities are included in the computation (as for the comdist function in the R package picante).


comdist.Swenson=function(comm, dis){
list.of.names=apply(comm, 1, function(x) names(x[x>0])) ### list of names species present in each community.

if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
    }
diss=as.matrix(dis)

return(as.dist(do.call(cbind,lapply(list.of.names, function(x) lapply(list.of.names, function(z) mean(diss[x,z], na.rm=T))) )))
}


### Function name=comdistnt.Swenson (=Dnn, MNND, Gamma+)
### This function computes the mean nearest distance among species between two communities. This metric is also called Dnn in Swenson which is equivalent of the Gamma + metric/100 defined by Clarke and Warwick 2006 or the MNND (mean nearest neighbour distance).
### this is the phylo beta version of the MNTD, comparing the closest relative of each species between communities.
### This function can also be used to estimate functional beta diversity using the mean nearest neighbour distance (MNND) approach, using any type of functional distance matrix (e.g. Gower distance matrix).
### This metric is influcence by the terminal branches of the phylogenetic tree (or close functional distance).
### Input:
# comm is a community data matrix (all the communities have to contain at least 1 species).
# dis should be a matrix containing the distances among species (e.g. coming from cophenetic() function for instance). If the dis object is a "dist" object
# the function atomatically convert it to a "matrix" object.

## The conspecific distances among communities are included in the computation (as for the comdist function in the R package picante).

comdistnt.Swenson=function(comm, dis){
list.of.names=apply(comm, 1, function(x) names(x[x>0])) ### list of names species present in each community.

if(class(list.of.names) == "matrix"){
list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
}
diss=as.matrix(dis)
matres=as.dist(do.call(cbind,lapply(list.of.names, function(x) lapply(list.of.names, function(z) {
if(length(x)==1 | length(z)==1){
mean(c(min(diss[x,z], na.rm=T), diss[x,z]), na.rm=T)
} else {
mean(c(apply(diss[x,z], MARGIN=1, min, na.rm=T), apply(diss[x,z], MARGIN=2, min, na.rm=T)), na.rm=T)
}
}
))))
return(matres)
}


