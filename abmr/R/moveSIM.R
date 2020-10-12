#'
#' Runs basic Brownian/Ornstein–Uhlenbeck agent-based model for multiple replicates.
#'
#' Here, agent mortality occurs when agent fails to achieve suitable raster values
#' at least 5 days (timesteps) in a row. Agent energy stores are not dynamic, so movement
#' speed isn't directly affected by quality of raster cells achieved. Results may be analyzed
#' with `moveVIZ()`. Relies on underlying function `moveSIM_helper`, which is not to be used
#' alone.
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
#' @param sigma Numeric, randomness parameter
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction
#' @param mot_y Numeric, movement motivation in y direction
#' @param modeled_species Object of class "species"
#' @param my_shapefile COME BACK
#' @param optimum Numeric, optimal environmental value
#' @param direction Character, mig direction, one of "N","S","E","W"
#' @param write_results Logical, save results to csv?
#' @param single_rast Logical, are you using a one-layer raster for all timesteps?
#'
#' @return
#' A (days X replicates) X 5 dataframe containing data on latitude, longitude,
#' day, agent ID, and distance traveled between each timestep (in km).
#' @examples
#' testing=moveSIM(replicates=1,days=27,env_rast=ndvi_raster, search_radius=375,
#' sigma=.4, dest_x=-100, dest_y=25, mot_x=1, mot_y=1,modeled_species=N_pop,
#' my_shapefile=NOAM,optimum=.5,direction="S",write_results=TRUE,single_rast=FALSE)
#' @export

moveSIM=function(replicates=200,days=27,env_rast=ndvi_raster, search_radius=375,
                 sigma, dest_x, dest_y, mot_x, mot_y, modeled_species, my_shapefile=NOAM,optimum,
                 direction="S",write_results=FALSE,single_rast=FALSE)

{
  if(nlayers(env_rast)==1 & single_rast==FALSE)
  {print("Single layer environmental raster with single_rast=FALSE specified.
         Please check your raster or change tosingle_rast=TRUE. Exiting
         function")
  return()}
  if(nlayers(env_rast)!=1 & single_rast==TRUE)
  {print("Multiple layer environmental raster with single_rast=TRUE specified.
    Using only first layer of raster")
    }

  my_env=env_rast-optimum

  long=data.frame(lon=numeric(),lat=numeric(),day=numeric(),
                  agent_id=character())
  sp_poly=my_shapefile

  for(i in 1:replicates){
    print("Agent #")
    Species=moveSIM_helper(sp=modeled_species,env=my_env,days=days,sigma=sigma,
                    dest_x=dest_x,dest_y=dest_y,mot_x=mot_x,mot_y=mot_y,
                    sp_poly=sp_poly,search_radius=search_radius,optimum=optimum,
                    direction=direction,single_rast=single_rast)
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
  if (write_results){
  currentDate=format(Sys.time(), "%d-%b-%Y %H.%M")
  file_name <- paste("moveSIM_results_",currentDate,".csv",sep="")
  write.csv(long,file_name)
  }

  missing_pct=sum(is.na(long$lon))/nrow(long)*100

  params=data.frame(replicates=replicates,days=days,
                    env_raster=deparse(substitute(env_raster)),
                    search_radius=search_radius,sigma=sigma,dest_x=dest_x,
                    dest_y=dest_y,mot_x=mot_x,mot_y=mot_y,
                    modeled_species=deparse(substitute(modeled_species)),
                    optimum=optimum,direction=direction,write_results=write_results,
                    single_rast=single_rast,missing_pct=missing_pct)
  return(list(results=long,run_params=params))
}
