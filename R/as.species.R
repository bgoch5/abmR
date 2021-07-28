#' Creates object of "species" class for input into moveSIM() and energySIM()
#' 
#' Here we define the agents whose movement we will be modeling. The user must 
#' indicate the geographical origin of the agents but can also optionally
#' specify values of morphological parameters, which will be standardized using
#' (observed-mean/sd). The standardized values will affect the movement
#' simulations by serving as a multiplier on "mot_x" and "mot_y" in the
#' direction specified: "Pos" or "Neg".
#' 
#' In example 1 below, we have a bird with origin (-100,26)
#' with observed wing chord length of 15 mm, while the population
#' mean for this measure is 10 mm with a SD of  2 mm (these birds are bigger
#' than average). We declare morphpar1sign = "Pos" because we assume longer
#' wingchord length leads to longer flight. Here, "morphpar2" represents mass,
#' and we want to model heavier than average birds. We assume that heavier birds
#' will fly longer distance, so specify `morphpar2sign` = `Pos` . If we assumed
#' that heavier birds will fly shorter distance, we would  specify 
#' `morphpar2sign`=`"Neg"` to indicate the inverse relationship.
#' 
#' @param x Species origin longitude value (degrees). Required.
#' @param y Species origin latitude value (degrees). Required.
#' @param morphpar1 Observed value for morphological parameter 1 (numeric)
#' @param morphpar1mean Population mean for morphological parameter 1 (numeric)
#' @param morphpar1sd Population standard deviation for morphological parameter 1 (numeric)
#' @param morphpar1sign Do higher values of morphpar1 lead to longer or shorter distances traveled each day? Specify "Pos" if longer and "Neg" if shorter.
#' @param morphpar2 Observed value for morphological parameter 2 (numeric)
#' @param morphpar2mean Population mean for morphological parameter 2 (numeric)
#' @param morphpar2sd Population standard deviation for morphological parameter 2 (numeric)
#' @param morphpar2sign Do higher values of morphpar2 lead to longer or shorter distances traveled each day? Specify "Pos" if longer and "Neg" if shorter.
#' 
#' @examples
#' 
#' # Example 1 -- Birds
#' # -Origin (-100,26)
#' #-Morphpar1=Wing Chord: Observed=15 mm, Pop. Mean=10, Pop. SD=2; Pos effect
#' #on movement speed
#' #-Morphpar2=Mass; Observed=15 g, Pop. Mean=6, Pop. SD=2; Pos effect on
#' # movement speed
#' 
#' my_species=as.species(x=-100,y=26,morphpar1=15,morphpar1mean=10,
#' morphpar1sd=2, morphpar1sign="Pos", morphpar2=7,morphpar2mean=6,
#' morphpar2sd=1,morphpar2sign="Pos")
#' 
#' # Example 2 -- Terrestrial Mammals
#' #-Origin (-90,40)
#' #-Morphpar1=Leg length: Observed=1.1 m, Pop. Mean=1, Pop. SD=.2; Pos effect
#' #on movement speed
#' #-Morphpar2=Mass; Observed=50 kg, Pop. Mean=55 kg, Pop. SD=10 kg; Neg effect
#' #on movement speed
#' 
#' my_species2=as.species(x=-90,y=40,morphpar1=1.1,morphpar1mean=1,
#' morphpar1sd=0.2, morphpar1sign="Pos", morphpar2=50,morphpar2mean=55,
#' morphpar2sd=10 ,morphpar2sign="Neg")
#' 
#' # Example 3 -- Unspecified agents
#' # -Origin (-90,40)
#' # -Not interested in modeling effect of morphology
#' 
#' my_species3=as.species(x=-90,y=40)
#' 
#' @export

as.species <- function(x=NA, y=NA, morphpar1=NA, morphpar1mean=NA, morphpar1sd=NA,morphpar1sign=NA,
                       morphpar2=NA,morphpar2mean=NA,morphpar2sd=NA,morphpar2sign=NA){
  sp <- setClass("species", slots=c(x="numeric", y="numeric", morphpar1="numeric", morphpar1mean="numeric",
                                    morphpar1sd="numeric",morphpar1sign="character",
                                    morphpar2="numeric", 
                                    morphpar2mean="numeric", 
                                    morphpar2sd="numeric",morphpar2sign="character"),
                 where=topenv(parent.frame()))
  if(is.na(morphpar1) & is.na(morphpar2)){
    res <- sp(x=x, y=y)} else{
      if (is.na(morphpar2)){
        cat("Error: Specified morphpar1 without specifying morphpar2")
        stop()
      }
      if (is.na(morphpar1)){
        cat("Error: Specified morphpar2 without specifying morphpar1")
        stop()
      }
      if (is.na(morphpar1mean)|is.na(morphpar1sd)){
        cat("Error: Specified morphpar1 without specifying morphpar1mean and/or morphpar1sd")
        stop()
      }
      if (is.na(morphpar2mean)|is.na(morphpar2sd)){
        cat("Error: Specified morphpar2 without specifying morphpar2mean and/or morphpar2sd")
        stop()
      }
      res <- sp(x=x, y=y, morphpar1=morphpar1, morphpar1mean=morphpar1mean, morphpar1sd=morphpar1sd,morphpar1sign=morphpar1sign,
                morphpar2=morphpar2,morphpar2mean=morphpar2mean,morphpar2sd=morphpar2sd,morphpar2sign=morphpar2sign)
    } 
  return(res)
}
