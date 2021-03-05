#' Runs more advanced Brownian/Ornstein–Uhlenbeck agent-based model for
#' multiple replicates.
#'
#' Here, agent mortality occurs when agent reaches energy = 0 (out of 100). Agent energy
#' stores are dynamic, and affect search area as a multiplier, so movement speed
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
#' @param optimum_lo Numeric, optimal environmental value (low)
#' @param optimum_hi Numeric, optimal environmental value (high)
#' @param init_energy Numeric, initial energy in interval (0,100]
#' @param direction Character, mig direction, one of "N","S","E","W"
#' @param mortality Logical, should low energy levels result in death?
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

energySIM2=function(replicates=200,days=27,env_rast=ndvi_raster, search_radius=375,
                sigma, dest_x, dest_y, mot_x, mot_y, modeled_species,
                optimum_lo,optimum_hi,init_energy,direction="S", mortality=TRUE,
                energy_adj=c(15,10,5,0,-5,-10,-15,-20),write_results=FALSE,
                single_rast=FALSE)

{
  sp=modeled_species
  if(optimum_hi<optimum_lo){
  print("opt_hi smaller than opt_lo--please check your work")
  stop()
  }

  if(init_energy<0 | init_energy>100)
  {print("Initial Energy should be such that 0<init_energy<=100")
    stop()}

  if(nlayers(env_rast)==1 & single_rast==FALSE)
  {print("Single layer environmental raster with single_rast=FALSE specified.
         Please check your raster or change to single_rast=TRUE. Exiting
         function")
    stop()}

  if(nlayers(env_rast)!=1 & single_rast==TRUE)
  {print("Multiple layer environmental raster with single_rast=TRUE specified.
    Using only first layer of raster")
  }
  
  if(length(sp@morphpar1)==1 & length(sp@morphpar2)==1){
  
  if (sp@morphpar1>sp@morphpar1mean+3.5*sp@morphpar1sd | sp@morphpar1<sp@morphpar1mean-3.5*sp@morphpar1sd)
  {cat("Error: Specified value for morphpar1 is greater than 3.5 SDs away from the
         specified mean, which is extremely unlikely. Consider adjusting morphpar1, morphpar1mean,
         or morphpar1SD. Function terminated.")
    stop()}
  
  if (sp@morphpar2>sp@morphpar2mean+3.5*sp@morphpar2sd | sp@morphpar2<sp@morphpar2mean-3.5*sp@morphpar2sd)
  {cat("Error: Specified value for morphpar2 is greater than 3.5 SDs away from the
         specified mean, which is extremely unlikely. Consider adjusting morphpar2, morphpar2mean,
         or morphpar2SD. Function terminated.")
    stop()}
  }
  optimum=(optimum_hi+optimum_lo)/2
  
  if(direction != "R"){
  my_env=env_rast-optimum
  }
  else{
    my_env=env_rast
    print("Direction=R specified--Raster will be ignored")
  }
  long=data.frame(lon=numeric(),lat=numeric(),energy=numeric(),
                  day=numeric(),agent_id=character())
  for(i in 1:replicates){
      Species=energySIM_helper2(sp = modeled_species,
                                           env_orig = env_rast,
                                           env_subtract=my_env,
                                           days = days,
                                           sigma = sigma,
                                           dest_x = dest_x,
                                           dest_y = dest_y,
                                           mot_x = mot_x,
                                           mot_y = mot_y,
                                           search_radius = search_radius,
                                           optimum_lo = optimum_lo,
                                           optimum_hi = optimum_hi,
                                           init_energy=init_energy,
                                           direction=direction,
                                           mortality=mortality,
                                           energy_adj=energy_adj,
                                           single_rast=single_rast)
      names(Species)=c("lon","lat","energy","Ending")
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
    currentDate=format(Sys.time(), "%d-%b-%Y %H.%M.%S")
    file_name <- paste("energySIM_results_",currentDate,".csv",sep="")
    write.csv(long,file_name)
  }

  missing_pct=sum(is.na(long$lon))/nrow(long)*100
  mortality_pct=length(which(long$energy==0))/replicates*100

  params=data.frame(replicates=replicates,days=days,
                    env_raster=deparse(substitute(env_raster)),
                    search_radius=search_radius,sigma=sigma,dest_x=dest_x,
                    dest_y=dest_y,mot_x=mot_x,mot_y=mot_y,
                    modeled_species=deparse(substitute(modeled_species)),
                    optimum_lo=optimum_lo,optimum_hi=optimum_hi,
                    direction=direction, mortality=mortality, write_results=write_results,
                    single_rast=single_rast,missing_pct=missing_pct,
                    mortality_pct=mortality_pct)
  return(list(results=long,run_params=params))
}
