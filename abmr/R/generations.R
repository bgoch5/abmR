#' Runs model for multiple replicates over multiple generations
#'
#' @import raster
#' @import sp
#' @import rgdal
#'
#' @param n_generations Integer, desired number of generations to run
#' @param replicates Integer, desired number of replicates per generation
#' @param days Integer, How many days (timesteps), would you like to model?
#' @param ndvi_rast COME BACK
#' @param my_sigma Numeric, randomness parameter
#' @param modeled_species Object of class "species"
#' @param my_NOAM COME BACK
#' @param mot_y Numeric, movement motivation in y direction
#' @param sp_poly Come back to this
#' @param current_generation Fed into function by function `generations`
#' @param search_radius Radius of semicircle to South of current location to search for next timestep (in km)
#'
#' @return A matrix
#' @examples
#' Come back to this
#' @export

generations=function(n_generations=3, replicates=200,days=27,ndvi_rast=ndvi_raster, my_sigma,
                     modeled_species,my_NOAM=NOAM)

{
  my_env=ndvi_rast-modeled_species@opt
  my_list=list()
  my_shapefiles=list()

  results=array(rep(0,replicates*days*2),c(replicates,days,2))

  for (j in 1:n_generations)
  {
    print(paste0("Starting Generation ", j))

    if(j==1){
      sp=my_NOAM
    }

    for(i in 1:replicates){

      Species=data.matrix(moveSIM(modeled_species,my_env,days,my_sigma, -99.11,19.15,1,1,sp_poly,j),375)

      if (length(Species)==days*2){
        results[i,,]=Species
      }
      else{
        my_matrix=matrix(rep(NA,days*2),nrow=days,ncol=2,byrow=TRUE)
        results[i,,]=my_matrix
      }

      if (i%%5 == 0) {
        print(paste0("Number of birds processed this gen.: ",i))
      }
    }

    my_list[[j]]=results

    aggregate_species=apply(results, c(2,3), mean,na.rm=T)
    SD_species=apply(results, c(2,3), sd,na.rm=T)
    lower_bound=aggregate_species-1.64*SD_species  #95 was too wide try 90%
    upper_bound=aggregate_species+1.64*SD_species

    my_df=rbind(lower_bound,upper_bound)
    my_df=na.omit(my_df)
    ch <- chull(my_df)
    coords <- my_df[c(ch, ch[1]), ]

    sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))

    my_shapefiles[[j]]=sp_poly

    prelim_results=list(this_gen=results, species_full=my_list,shapefiles=my_shapefiles)
    save(prelim_results,file=paste0("PrelimResults_Gen",j,"_",Sys.Date(),".RData"))
  }
  results_to_return=list(last_gen=results,
               species_full=my_list,shapefiles=my_shapefiles)
  return(results_to_return)
}
