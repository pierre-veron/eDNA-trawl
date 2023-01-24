# function to extract the Average and sd and 0.025 and 97.5 quantiles of different statistics present in a tables across multiple rds tables 
# stored in a folder.

# path = "Results/PhyloAlpha/100trees/"

Summary.Stat = function(path = NULL, prefix = NULL){

if(is.null(prefix)){
NewInd_100tr = list.files(path)
} else {
NewInd_100tr = list.files(path)
NewInd_100tr = NewInd_100tr[grep(prefix, NewInd_100tr, fixed = TRUE)]
}

# Convert all the rds tables into a single list.
NewInd_list = list()
i = 1
for(i in 1:length(NewInd_100tr)){
#cat("i = ", i, "\n")
NewInd_list[[i]] = readRDS(paste(path, NewInd_100tr[i], sep = ""))
}

nb.met = dim(NewInd_list[[1]])[2]

# Compute the mean
Met.Mean = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
apply(a.100, 1, function(x) mean(x, na.rm = TRUE))
})))
names(Met.Mean) = paste(names(NewInd_list[[1]]), ".mean", sep = "")
row.names(Met.Mean) = row.names(NewInd_list[[1]])

# compute the sd.
Met.Sd = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
apply(a.100, 1, function(x) sd(x, na.rm = TRUE))
})))
names(Met.Sd) = paste(names(NewInd_list[[1]]), ".sd", sep = "")
row.names(Met.Sd) = row.names(NewInd_list[[1]])

# 0.025 quantiles
Met.Quant0.025 = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
apply(a.100, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
})))
names(Met.Quant0.025) = paste(names(NewInd_list[[1]]), ".Quant0.025", sep = "")
row.names(Met.Quant0.025) = row.names(NewInd_list[[1]])

# 0.975 quantiles
Met.Quant0.975 = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
apply(a.100, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
})))
names(Met.Quant0.975) = paste(names(NewInd_list[[1]]), ".Quant0.975", sep = "")
row.names(Met.Quant0.975) = row.names(NewInd_list[[1]])

return(data.frame(Met.Mean, Met.Sd, Met.Quant0.025, Met.Quant0.975))
}
