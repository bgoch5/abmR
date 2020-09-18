#' Runs model for multiple replicates over multiple generations
#'
#' @import raster
#' @import sp
#' @import rgdal
#' @import swfscMisc
#'
#' @param replicates Integer, desired number of replicates per generation
#' @param days Integer, How many days (timesteps), would you like to model?
#' @param env_rast Rasterstack or Rasterbrick with number of layers equal to n_days
#' @param search_radius Radius of semicircle to South of current location to search for next timestep (in km)
#' @param my_sigma Numeric, randomness parameter
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction
#' @param mot_y Numeric, movement motivation in y direction
#' @param modeled_species Object of class "species"
#' @param my_shapefile COME BACK

#'
#' @return
#' 1. species_full: A (days x replicates X n_generations) X 5 dataframe containing the columns
#' lon, lat, day, gen, and bird_id (e.g. 1_2 means second bird of first gen.)
#'
#' 2. last_gen: species_full subsetted to last gen.
#'
#' 3. shapefiles: A list of shapefiles produced, one for each generation
#' @examples
#' my_results=genSIM(n_generations = 2,replicates = 3,my_sigma = 0.7,env_rast = ndvi_raster,
#' modeled_species = N_pop, dest_x = -100,dest_y = 25,my_shapefile = NOAM,
#' mot_x = 0.9,mot_y = 0.8,days = 27)
#' my_results$last_gen
#' my_results$species_full
#' my_results$shapefiles
#' @export

moveSIM=function(replicates=200,days=27,env_rast=ndvi_raster, search_radius=375,
                 my_sigma, dest_x, dest_y, mot_x, mot_y, modeled_species, my_shapefile=NOAM,optimum)

{
  my_env=env_rast-optimum

  long=data.frame(lon=numeric(),lat=numeric(),day=numeric(),
                  agent_id=character())
  sp_poly=my_shapefile

  for(i in 1:replicates){
    print(i)
    Species=moveSIM_helper(sp=modeled_species,env=my_env,days=days,sigma=my_sigma,
                    dest_x=dest_x,dest_y=dest_y,mot_x=mot_x,mot_y=mot_y,
                    sp_poly=sp_poly,search_radius=search_radius,optimum=optimum)
    names(Species)=c("lon","lat")
    Species$day=1:nrow(Species)
    Species$agent_id=paste("Agent",as.character(i),sep="_")

    if (nrow(Species)==days){
      long=rbind(long,Species)
    }
    if (i%%5 == 0) {
      print(paste0("Number of agents processed: ",i))
    }
  }

  long[,"distance"]=NA
  for (i in 2:nrow(long)){
    if(i%%days!=0){
      long$distance[i]<-distHaversine(long[(i-1),1:2], long[i,1:2])/1000
    }
  }
  return(long)
}

