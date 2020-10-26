#' A class for the input species
#'
#' @slot x Species origin longitude value (degrees)
#' @slot y Species origin latitude value (degrees)
#' @slot mass Species mass
#' @slot wing Species wing
#'
species <- setClass("species", slots=c(x="numeric", y="numeric",mass="numeric",
                                       wing="numeric"))
