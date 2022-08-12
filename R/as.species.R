#' Creates object of "species" class for input into moveSIM() and energySIM()
#' 
#' Here we define the geographical origin of the agents whose movement we will be modeling.
#' 
#' @param x Species origin longitude value (degrees). Required.
#' @param y Species origin latitude value (degrees). Required.
#' @return
#' Object of class "species" for input into moveSIM() or energySIM()
#' @examples
#' myspecies <- as.species(x = -98.7, y = 34.7)
#' @export

as.species <- function(x=NA, y=NA){
  sp <- setClass("species", slots=c(x="numeric", y="numeric"),
                 where=topenv(parent.frame()))
      res <- sp(x=x, y=y)
  return(res)
}
