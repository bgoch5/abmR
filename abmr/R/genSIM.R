#' Runs model for multiple replicates over multiple generations
#'
#' @import raster
#' @import sp
#' @import rgdal
#' @import swfscMisc
#'
#' @param n_generations Integer, desired number of generations to run
#' @param replicates Integer, desired number of replicates per generation
#' @param days Integer, How many days (timesteps), would you like to model?
#' @param env_rast COME BACK
#' @param search_radius Radius of semicircle to South of current location to search for next timestep (in km)
#' @param my_sigma Numeric, randomness parameter
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction
#' @param mot_y Numeric, movement motivation in y direction
#' @param modeled_species Object of class "species"
#' @param my_shapefile COME BACK
#' @param progress_save logical, Whether or not to save results after each gen;
#' default=FALSE; makes sense to use if you have a huge run and you want
#' to be sure that intermediate results are saved
#'
#' @return
#' species_full: A days*replicates*n_generations X 5 dataframe containing the columns
#' lon, lat, day, gen, and bird_id (e.g. 1_2 means second bird of first gen.)
#' last_gen: species_full subsetted to last gen.
#' shapefiles: A list of shapefiles produced, one for each generation
#' @examples
#' my_results=genSIM(n_generations=2,replicates=3,my_sigma=.7,env_rast=ndvi_raster,
#' modeled_species=N_pop,dest_x=-90,dest_y=34,my_shapefile=NOAM,
#' mot_x=.9,mot_y=.8,days=27)
#' my_results$last_gen
#' my_results$species_full
#' my_results$shapefiles
#' @export

genSIM=function(n_generations=3, replicates=200,days=27,env_rast=ndvi_raster, search_radius=375,
                my_sigma, dest_x, dest_y, mot_x, mot_y, modeled_species, my_shapefile=NOAM,
                progress_save=FALSE)

{
  my_env=env_rast-modeled_species@opt
  my_list=list()
  my_shapefiles=list()

  results=array(rep(0,replicates*days*2),c(replicates,days,2))
  long=data.frame(lon=numeric(),lat=numeric(),day=numeric(),gen=numeric(),
                     bird_id=character())

  for (j in 1:n_generations)
  {
    print(paste0("Starting Generation ", j))

    if(j==1){
      sp_poly=my_shapefile
    }

    for(i in 1:replicates){

      Species=moveSIM(sp=modeled_species,env=my_env,n=days,sigma=my_sigma,
      dest_x=dest_x,dest_y=dest_y,mot_x=mot_x,mot_y=mot_y,
      sp_poly=sp_poly,current_gen=j,search_radius=search_radius)
      names(Species)=c("lon","lat")
      Species$day=1:nrow(Species)
      Species$gen=j
      Species$bird_id=paste(as.character(j),as.character(i),sep="_")

      if (nrow(Species)==days){
        long=rbind(long,Species)
        results[i,,]=data.matrix(Species[,c("lon","lat")])
      }
      else{
        my_matrix=matrix(rep(NA,days*2),nrow=days,ncol=2,byrow=TRUE)
        results[i,,]=my_matrix
      }

      if (i%%5 == 0) {
        print(paste0("Number of birds processed this gen.: ",i))
      }
    }

    #my_list[[j]]=results

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

    my_last_gen=long[long$gen==n_generations,]
    rownames(my_last_gen)=NULL

    if(progress_save){
    prelim_results=list(species_full=long,shapefiles=my_shapefiles)
    save(prelim_results,file=paste0("PrelimResults_Gen",j,"_",Sys.Date(),".RData"))
  }}
  results_to_return=list(last_gen=my_last_gen,
               species_full=long,shapefiles=my_shapefiles)
  return(results_to_return)
}
