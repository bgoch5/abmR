#' A class for the input species
#'
#' @slot x Species origin longitude value (degrees)
#' @slot y Species origin latitude value (degrees)
#' @slot p1 Observed value for morphological parameter 1 (numeric)
#' @slot p1mean Population mean for morphological parameter 1 (numeric)
#' @slot p1sd Population standard deviation for morphological parameter 1 (numeric)
#' @slot p1sign Do higher values of p1 lead to faster or slower movement?
#' Specify "Pos" if faster and "Neg" if slower.
#' @slot p2 Observed value for morphological parameter 2 (numeric)
#' @slot p2mean Population mean for morphological parameter 2 (numeric)
#' @slot p2sd Population standard deviation for morphological parameter 2 (numeric)
#' @slot p2sign Do higher values of p2 lead to faster or slower movement?
#' Specify "Pos" if faster and "Neg" if slower.
#' @examples
#' my_species=species(x=-100,y=26,p1=15,p1mean=10,p1sd=2,p1sign="Pos",
#' p2=7,p2mean=6,p2sd=1,p2sign="Pos")

species <- setClass("species", slots=c(x="numeric", y="numeric",
                                       p1="numeric", p1mean="numeric", p1sd="numeric",p1sign="character",
                                       p2="numeric", p2mean="numeric", p2sd="numeric",p2sign="character"))

