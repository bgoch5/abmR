
#' Creates a plot of energySIM() results
#'
#' Compares results with straight line path (null model)
#'
#' @import raster
#' @import sp
#' @import rgdal
#' @import swfscMisc
#' @import rnaturalearth
#' @import rnaturalearthdata
#' @import ggplot2
#'
#'
#' @param data Data to be plotted, this object should be the output from
#' moveSIM().
#' @param title Title for the plot that is output.
#'
#' @return
#' A plot showing model output compared to a dashed line that represents straight line
#' movement from the starting point to the final destination.
#' @examples
#' 
#' 1. Run moveSIM()
#' 
#' my_result=moveSIM(replicates=1,days=27,env_rast=ndvi_raster, search_radius=375,
#' sigma=.4, dest_x=-100, dest_y=25, mot_x=1, mot_y=1,modeled_species=N_pop,
#' my_shapefile=NOAM,optimum=.5,direction="S",write_results=TRUE,single_rast=FALSE)
#' 
#' 2. Run energySIM() on your result
#' moveVIZ(my_result,title="Visualizing MoveSIM results",type="plot",aspect_ratio=4/3,
#' label=TRUE)
#'
#' @export

moveVIZ=function(data,title="Visualizing MoveSIM results")
{
    world <- ne_countries(scale = "medium", returnclass = "sf")
    start.p <- cbind(data$results[1,1], data$results[1,2])
    # Generalize this soon
    start.p.df <- as.data.frame(start.p)
    colnames(start.p.df)[1:2] = c("Lon", "Lat")
    run = "ideal"
    start.p.df <- cbind(start.p.df, run)
    end.p <- cbind(data$run_params["dest_x"],data$run_params["dest_y"])
    end.p.df <- as.data.frame(end.p)
    colnames(end.p.df)[1:2] = c("Lon", "Lat")
    end.p.df <- cbind(end.p.df, run)
    ideal.df <- rbind(start.p.df, end.p.df)
    t.energy.res <- data$results
    my_xlim = c((min(t.energy.res$lon,na.rm=T)-6), (max(t.energy.res$lon,na.rm=T)+6))
    my_ylim = c((min(t.energy.res$lat,na.rm=T)-6), (max(t.energy.res$lat,na.rm=T)+6))
    myplot=ggplot(data = world) +
      geom_sf() +
      coord_sf(xlim = my_xlim, ylim = my_ylim, expand = FALSE) +
      geom_path(data = t.energy.res,
                aes(x=lon, y=lat),
                color = "blue", size = 0.6, alpha = 0.4, lineend = "round") +
      geom_path(data = ideal.df,
                aes(x=Lon, y=Lat),
                color = "black", size = 1.2, alpha = 1, linetype = 2) +
      ggtitle(title)
    return(myplot)
}
