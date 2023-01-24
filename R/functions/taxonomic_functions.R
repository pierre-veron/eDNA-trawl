### ----------------------------------------------------------------------------
#    Defines some functions that are useful to manipulate taxonomic information.
#
# Author : Pierre Veron, pierre.veron.2017@polytechnique.org
# Created 2022/02/21
### ----------------------------------------------------------------------------


#' class_to_species
#' Retrieves the name of a species from a classification object (from the library
#' taxize). If the classification object correspond to a rank higher than species
#' it will return NA. If the classification correspond to a Not found taxa, it 
#' returns "not_found".
#' 
#' @param class a list object from a classification query 
#'
#' @return string : name of the species or NA or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "bold")
#' class_to_species(class)
#' [1] "Sardina pilchardus"
#' 
#' class <- classification("Clupeidae", db = "bold")
#' class_to_species(class)
#' [1] NA
#' 
#' class <- classification("foo", db = "bold")
#' class_to_species(class)
#' [1] "not_found"
#' 
class_to_species <- function(class) {
  if (typeof(class[[1]]) == "logical") { # means empty classification
    sp <- "not_found"
  } else if (nrow(class[[1]] %>% filter(rank == "species")) == 0) {
    sp <- NA
  } else {
    sp <- (class[[1]] %>% filter(rank == "species"))$name
  }
  sp
}

#' class_to_genus
#' Retrieves the name of a genus from a classification object (from the library
#' taxize). If the classification object correspond to a rank higher than genus
#' it will return NA. If the classification correspond to a Not found taxa, it 
#' returns "not_found".
#' 
#' @param class a list object from a classification query 
#'
#' @return string : name of the genus or NA or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "bold")
#' class_to_genus(class)
#' [1] "Sardina"
#' 
#' class <- classification("Clupeidae", db = "bold")
#' class_to_genus(class)
#' [1] NA
#' 
#' class <- classification("foo", db = "bold")
#' class_to_genus(class)
#' [1] "not_found"
#' 
class_to_genus <- function(class) { 
  if (typeof(class[[1]]) == "logical") { # means empty classification
    ge <- "not_found"
  } else if (nrow(class[[1]] %>% filter(rank == "genus")) == 0) {
    ge <- NA
  } else {
    ge <- (class[[1]] %>% filter(rank == "genus"))$name
  }
  ge
}

#' class_to_family
#' Retrieves the name of a family from a classification object (from the library
#' taxize). If the classification object correspond to a rank higher than family
#' it will return NA. If the classification correspond to a Not found taxa, it 
#' returns "not_found".
#' 
#' @param class a list object from a classification query 
#'
#' @return string : name of the species or NA or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "bold")
#' class_to_family(class)
#' [1] "Clupeidae"
#' 
#' class <- classification("Actinopterygii", db = "bold")
#' class_to_family(class)
#' [1] NA
#' 
#' class <- classification("foo", db = "bold")
#' class_to_family(class)
#' [1] "not_found"
#' 
class_to_family <- function(class) { 
  if (typeof(class[[1]]) == "logical") { # means empty classification
    fa <- "not_found"
  } else if (nrow(class[[1]] %>% filter(rank == "family")) == 0) {
    fa <- NA
  } else {
    fa <- (class[[1]] %>% filter(rank == "family"))$name
  }
  fa
}

#' class_to_order
#' Retrieves the name of an order from a classification object (from the library
#' taxize). If the classification object correspond to a rank higher than order
#' it will return NA. If the classification correspond to a Not found taxa, it 
#' returns "not_found".
#' 
#' @param class a list object from a classification query 
#'
#' @return string : name of the species or NA or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "bold")
#' class_to_order(class)
#' [1] "Clupeiformes"
#' 
#' class <- classification("Actinopterygii", db = "bold")
#' class_to_order(class)
#' [1] NA
#' 
#' class <- classification("foo", db = "bold")
#' class_to_order(class)
#' [1] "not_found"
#' 
class_to_order <- function(class) { 
  if (typeof(class[[1]]) == "logical") { # means empty classification
    or <- "not_found"
  } else if (nrow(class[[1]] %>% filter(rank == "order")) == 0) {
    or <- NA
  } else {
    or <- (class[[1]] %>% filter(rank == "order"))$name
  }
  or
}

#' class_to_class
#' Retrieves the name of a class from a classification object (from the library
#' taxize). If the classification object correspond to a rank higher than class
#' it will return NA. If the classification correspond to a Not found taxa, it 
#' returns "not_found".
#' 
#' @param class a list object from a classification query 
#'
#' @return string : name of the species or NA or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "bold")
#' class_to_class(class)
#' [1] "Actinopterygii"
#' 
#' class <- classification("Chordata", db = "bold")
#' class_to_class(class)
#' [1] NA
#' 
#' class <- classification("foo", db = "bold")
#' class_to_class(class)
#' [1] "not_found"
#' 
class_to_class <- function(class) {
  if (typeof(class[[1]]) == "logical") { # means empty classification
    cl <- "not_found"
  } else if (nrow(class[[1]] %>% filter(rank == "class")) == 0) {
    cl <- NA
  } else {
    cl <- (class[[1]] %>% filter(rank == "class"))$name
  }
  cl
}

#' lower_rank
#' Retrieves the lower taxonomic level (rank) of a "taxize" classification.
#' If the classification comes from an unsuccessful query, it returns "not_found"
#' 
#' @param class a list object from a classification query 
#'
#' @return string : name of the tax. rank (e.g. "species", "genus", "family"...)
#' or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "bold")
#' lower_rank(class)
#' [1] "species"
#' 
#' class <- classification("Sardina", db = "bold")
#' lower_rank(class)
#' [1] "genus"
#' 
#' class <- classification("Clupeidae", db = "bold")
#' lower_rank(class)
#' [1] "family"
#' 
#' class <- classification("foo", db = "bold")
#' lower_rank(class)
#' [1] "not_found"
#' 
lower_rank <- function(class) { 
  rks = c("class", "order","family","genus","species")
  if (typeof(class[[1]]) == "logical") { # means empty classification
    rk <- "not_found"
  } else {
    rk <- NA
    for (test_rank in rks) {
      if (nrow(class[[1]] %>% filter(rank == test_rank)) > 0){
        rk <- test_rank
      } 
    }
  }
  rk
}


#' get_specific_rank
#' Convert a taxize classification into a string readable as a tree. 
#' Output will be of format : 
#' superkingdom/kingdom/phylum/class/order/family/species
#' @param class a list object from a classification query 
#' @param rk a string with the name of the rank (e.g. species or phylum)
#'
#' @return string : path of the species/taxa or NA or "not_found"
#' @export
#'
#' @examples
#' class <- classification("Sardina pilchardus", db = "ncbi")
#' get_specific_rank(class, "genus")
#' [1] "Sardina"
#' 
#' get_specific_rank(class, "order")
#' [1] "Clupeiformes"
#' 
get_specific_rank <- function(class,rk) { 
  if (typeof(class[[1]]) == "logical") { # means empty classification
    out <- "not_found"
  } else if (nrow(class[[1]] %>% filter(rank == rk)) == 0) {
    out <- NA
  } else {
    out <- (class[[1]] %>% filter(rank == rk))$name
  }
  out
}