#' Runs more advanced Brownian/Ornstein–Uhlenbeck agent-based model for
#' multiple replicates.
#'
#' Here, agent mortality occurs when agent reaches energy = 0 (out of 100). Agent energy
#' stores are not dynamic, and effect search area as a multiplier, so movement speed
#' is directly affected by the quality of raster cells achieved. Results may be analyzed
#' with `energyVIZ()`. Relies on underlying function `energySIM_helper`, which is
#' not to be used alone.
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
#' @param optimum_lo Numeric, optimal environmental value (low)
#' @param optimum_hi Numeric, optimal environmental value (high)
#' @param init_energy Numeric, initial energy in interval (0,100]
#' @param direction Character, mig direction, one of "N","S","E","W"
#' @param write_results Logical, save results to csv?
#' @param single_rast Logical, are you using a one-layer raster for all timesteps?
#'
#' @return
#' #' A (days X replicates) X 7 dataframe containing data on latitude, longitude, energy,
#' day, agent ID, distance traveled between each timestep (in km), and change in
#' energy from last timestep.
#' @examples
#' testing=energySIM(replicates=1,days=27,env_rast=ndvi_raster, search_radius=375,
#' sigma=.4, dest_x=-100, dest_y=25, mot_x=1, mot_y=1,
#' modeled_species=N_pop, my_shapefile=NOAM, optimum_lo=.4,optimum_hi=.6,
#' init_energy=100,direction="S",write_results=TRUE,single_rast=FALSE)
#' @export

energySIM=function(replicates=200,days=27,env_rast=ndvi_raster, search_radius=375,
                sigma, dest_x, dest_y, mot_x, mot_y, modeled_species, my_shapefile=NOAM,
                optimum_lo,optimum_hi,init_energy,direction="S",
                write_results=FALSE,single_rast=FALSE)

{
  if(optimum_hi<optimum_lo){
  print("my_opt_hi smaller than my_opt_lo--please check your work")
  }

  if(init_energy<0 | init_energy>100)
  {print("Initial Energy should be such that 0<init_energy<=100")}

  if(nlayers(env_rast)==1 & single_rast==FALSE)
  {print("Single layer environmental raster with single_rast=FALSE specified.
         Please check your raster or change tosingle_rast=TRUE. Exiting
         function")
    return()}

  if(nlayers(env_rast)!=1 & single_rast==TRUE)
  {print("Multiple layer environmental raster with single_rast=TRUE specified.
    Using only first layer of raster")
  }
  optimum=(optimum_hi+optimum_lo)/2
  my_env=env_rast-optimum
  long=data.frame(lon=numeric(),lat=numeric(),energy=numeric(),
                  day=numeric(),agent_id=character())
  for(i in 1:replicates){
      Species=energySIM_helper(sp = modeled_species,
                                           env_orig = env_rast,
                                           env_subtract=my_env,
                                           days = days,
                                           sigma = sigma,
                                           dest_x = dest_x,
                                           dest_y = dest_y,
                                           mot_x = mot_x,
                                           mot_y = mot_y,
                                           sp_poly = my_shapefile,
                                           search_radius = search_radius,
                                           optimum_lo = optimum_lo,
                                           optimum_hi = optimum_hi,
                                           init_energy=init_energy,
                                           direction=direction,
                                           single_rast=single_rast)
      names(Species)=c("lon","lat","energy")
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
    if(long$day[i]!=1){
    long$distance[i]<-distHaversine(long[(i-1),1:2], long[i,1:2])/1000
    }
  }
  long[,"delta_energy"]=NA
  for (i in 2:nrow(long)){
    if(long$day[i]!=1){
      long$delta_energy[i]<-long[i,3]-long[(i-1),3]
    }
  }
  if (write_results){
    currentDate=format(Sys.time(), "%d-%b-%Y %H.%M")
    file_name <- paste("energySIM_results_",currentDate,".csv",sep="")
    write.csv(long,file_name)
  }
  return(long)
}
