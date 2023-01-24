#### Compute four Functional indices:

### With thses four indidces the idea is to use the functional space intact by using the full functional distance matrix, 
# no need for reducing the functional space using the x best functional axes.

## 2 indices documenting the divergence (=dispersion)
# mPwFD = Mean pairwise functional distance, = average functional distance among two species randomly selected in the community.
# mNFD = Mean nearest functional distance, this index provides an idea about how unique a species or an assemblage is in comparison of the other.

## 2 indices documenting the regularity
# Var.mPwFD =  variance of the mean pariwise functional distances. If the variance is high it means you have the presence of species closely located in the functional space and species that are far appart.
# Var.mNFD =  variance of the nearest functional distances. 


### input for the 4 functions.

# dis= is  "dist" object, so a distance object typically a distance matrix based on euclidian distance or Gower distances.
# CommData =  is a data.frame with the information about the community data, the sampling sites are in rows and the species are in column. Abundance or presence/absence (0/1) can be provided.


### Calculating FDis manually which is also MPD.
# Functional dispersion = Mean pairwise distance = functional dispersion (Laliberte &  Legendre, 2010).
mPwFD = function(dis=NULL, CommData=NULL){
  list.of.names=apply(CommData, 1, function(x) names(x[x>0]))
  if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
  }
  diss=as.matrix(dis)
  unlist(lapply(list.of.names, function(x) mean(as.dist(diss[x,x]), na.rm=T)))
}

### Calculating the variation of the MPD (Var.mPwFD).
Var.mPwFD = function(dis=NULL, CommData=NULL){
  diss=as.matrix(dis)
  list.of.names=apply(CommData, 1, function(x) names(x[x>0]))
  if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
  }
  unlist(lapply(list.of.names, function(x) {
    com.dist=as.dist(diss[x,x])
    ((sum(com.dist^2))/((length(x)*(length(x)-1))/2))-(mean(com.dist)^2)
  }))
}

### Mean nearest taxon distance (mntd) See Swenson 2014. ="Dispersion index"
mNFD = function(dis=NULL, CommData=NULL){
  diss=as.matrix(dis)
  diag(diss) = NA
  list.of.names=apply(CommData, 1, function(x) names(x[x>0]))
  if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
  }
  cnt = unlist(lapply(list.of.names, function(i) length(i)))
  ID = seq(1, length(list.of.names))
  a = which(cnt < 2)
  if(length(a) > 0){
    list.of.names = list.of.names[-a]
    res = unlist(lapply(list.of.names, function(x){
      mean(apply(diss[x,x], MARGIN=1, min, na.rm=T))
    }))
    IDres = c(ID[-a], a)
    res1 = c(res, rep(NA, length(a)))
    aa = cbind(IDres, res1)
    aa[order(aa[,1]),2]
  } else {
    unlist(lapply(list.of.names, function(x){
      mean(apply(diss[x,x], MARGIN=1, min, na.rm=T))
    }))
  }
}

### Variance of MNTD (Tucker et al. 2016). ="Regularity index"
Var.mNFD = function(dis=NULL, CommData=NULL){
  diss=as.matrix(dis)
  diag(diss) = NA
  list.of.names=apply(CommData, 1, function(x) names(x[x>0]))
  if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
  }
  cnt = unlist(lapply(list.of.names, function(i) length(i)))
  ID = seq(1, length(list.of.names))
  a = which(cnt < 2)
  if(length(a) > 0){
    list.of.names = list.of.names[-a]
    res = unlist(lapply(list.of.names, function(x) {
      sum((apply(diss[x,x], MARGIN=1, min, na.rm=T)-mean(apply(diss[x,x], MARGIN=1, min, na.rm=T)))^2)/length(x)
    }))  
    IDres = c(ID[-a], a)
    res1 = c(res, rep(NA, length(a)))
    aa = cbind(IDres, res1)
    aa[order(aa[,1]),2]
  } else {
    unlist(lapply(list.of.names, function(x) {
      sum((apply(diss[x,x], MARGIN=1, min, na.rm=T)-mean(apply(diss[x,x], MARGIN=1, min, na.rm=T)))^2)/length(x)
    }))  
  }
}

#END FUNCTIONS
#############################################################################################





## Fof debug
#sp1 = c(1,1,0,0)
#sp2 = c(1,1,0,0)
#sp3 = c(1,1,0,0)
#CommData = data.frame(sp1, sp2, sp3)

#tr1 = rnorm(4, 20, 10)
#tr2 = rnorm(4, 12, 4)
#tr3 = rnorm(4, 15, 4)
#tr4 = rnorm(4, 12, 8)

#dataTr = data.frame(tr1, tr2, tr3, tr4)
#row.names(dataTr) = c("sp1", "sp2", "sp3", "sp4")

#dis = dist(dataTr)
#diss=as.matrix(dis)

#mPwFD(dis=dis, CommData=CommData)
#Var.mPwFD(dis=dis, CommData=CommData)
#mNFD(dis=dis, CommData=CommData)
#Var.mNFD(dis=dis, CommData=CommData)

