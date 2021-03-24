#' A class for the input species. 
#' 
#' The units for the morphological parameters are irrelevant as long as
#' they are consistent within `morphpar1` and within `morphpar2`. Information for each
#' parameter will be standardized (observed-mean/sd), and will affect
#' the speed of migration in the direction specified. Slots x and y are only required
#' arguments, morphpar1 through morphpar2 sign are optional--leave blank if you'd
#' like to run a simple model that doesn't consider morphology.
#' 
#' In example 1 below, we have a bird with destination (-100,26)
#' with observed wing chord length of 15 mm, while the population
#' mean for this measure is 10 mm with a SD of  2 mm (these birds are bigger than
#' average). We declare `morphpar1sign`=`"Pos"` because we believe longer wingchord
#' length leads to faster flight. Here, `morphpar2` represents mass, and we want to model
#' heavier than average birds. We believe that heavier birds will fly faster so
#' specify `morphpar2sign`=`"Pos"`. If we believed that heavier birds will fly slower,
#' we would  specify `morphpar2sign`=`"Neg"` to indicate the inverse relationship.
#' 
#' @slot x Species origin longitude value (degrees)
#' @slot y Species origin latitude value (degrees)
#' @slot morphpar1 Observed value for morphological parameter 1 (numeric)
#' @slot morphpar1mean Population mean for morphological parameter 1 (numeric)
#' @slot morphpar1sd Population standard deviation for morphological parameter 1 (numeric)
#' @slot morphpar1sign Do higher values of morphpar1 lead to faster or slower movement? Specify "Pos" if faster and "Neg" if slower.
#' @slot morphpar2 Observed value for morphological parameter 2 (numeric)
#' @slot morphpar2mean Population mean for morphological parameter 2 (numeric)
#' @slot morphpar2sd Population standard deviation for morphological parameter 2 (numeric)
#' @slot morphpar2sign Do higher values of morphpar2 lead to faster or slower movement? Specify "Pos" if faster and "Neg" if slower.
#' 
#' @examples
#' 
#' Example 1
#' -Destination (-100,26)
#' -Morphpar1=Wing Chord: Observed=15 mm, Pop. Mean=10, Pop. SD=2; Pos effect
#' on movement speed
#' -Morphpar2=Mass; Observed=15 g, Pop. Mean=6, Pop. SD=2; Pos effect on
#' movement speed
#' 
#' my_species=as.species(x=-100,y=26,morphpar1=15,morphpar1mean=10,morphpar1sd=2,
#' morphpar1sign="Pos", morphpar2=7,morphpar2mean=6,morphpar2sd=1,morphpar2sign="Pos")
#' 
#' Example 2
#' -Destination (-90,40)
#' -Morphpar1=Leg length: Observed=1.1 m, Pop. Mean=1, Pop. SD=.2; Pos effect
#' on movement speed
#' -Morphpar2=Mass; Observed=50 kg, Pop. Mean=55 kg, Pop. SD=10 kg; Neg effect on
#' movement speed
#' 
#' my_species2=as.species(x=-90,y=40,morphpar1=1.1,morphpar1mean=1,morphpar1sd=0.2,
#' morphpar1sign="Pos", morphpar2=50,morphpar2mean=55,morphpar2sd=10
#' ,morphpar2sign="Neg")
#' 
#' Example 3
#' -Destination (-90,40)
#' -Not interested in modeling effect of morphology
#' 
#' my_species3=as.species(x=-90,y=40)
#' 
#' @export

as.species =  function(x=NA, y=NA, morphpar1=NA, morphpar1mean=NA, morphpar1sd=NA,morphpar1sign=NA,
                       morphpar2=NA,morphpar2mean=NA,morphpar2sd=NA,morphpar2sign=NA)
  {
  sp <- setClass("species", slots=c(x="numeric", y="numeric", morphpar1="numeric", morphpar1mean="numeric",
                                    morphpar1sd="numeric",morphpar1sign="character",
                                    morphpar2="numeric", 
                                    morphpar2mean="numeric", 
                                    morphpar2sd="numeric",morphpar2sign="character"))
  res <- sp(x=x, y=y,
            morphpar1=morphpar1, morphpar1mean=morphpar1mean, morphpar1sd=morphpar1sd,morphpar1sign=morphpar1sign,
            morphpar2=morphpar2,morphpar2mean=morphpar2mean,morphpar2sd=morphpar2sd,morphpar2sign=morphpar2sign)
  return(res)
}