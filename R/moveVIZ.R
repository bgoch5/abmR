
#' Creates a plot of energySIM() results
#'
#' Compares results with straight line path between origin and destination.
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
#' @param type "plot" or "summary_table", default "plot".
#' @param title Title for the plot that is output.
#'
#' @return
#' A plot showing model output compared to a dashed line that represents straight line
#' movement from the starting point to the final destination.
#' @examples
#' 
#' 1. Run moveSIM()
#' 
#' EX2=moveSIM(replicates=5,days=27,env_rast=ndvi_raster, search_radius=550,
#' sigma=.1, dest_x=-108.6, dest_y=26.2, mot_x=.8, mot_y=.8,modeled_species=pabu.pop,optimum=.6, n_failures=5, fail_thresh=.40,
#'  direction="S",write_results=TRUE,single_rast=FALSE,mortality = T)
#' 
#' 2. Run energySIM() on your result
#' moveVIZ(EX2,title="Visualizing MoveSIM results",type="plot",aspect_ratio=4/3,
#' label=TRUE)
#'
#' @export

moveVIZ=function(data, type="plot", title="MoveSIM results")
{
  if(type=="plot"){
    dest_x=data$run_params$dest_x
    dest_y=data$run_params$dest_y
    world <- ne_countries(scale = "medium", returnclass = "sf")
    start.p <- cbind(data$results[1,"lon"], data$results[1,"lat"])
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
    if(dest_x!=999 & dest_y!=999){
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
    }
    else{
      myplot=ggplot(data = world) +
        geom_sf() +
        coord_sf(xlim = my_xlim, ylim = my_ylim, expand = FALSE) +
        geom_path(data = t.energy.res,
                  aes(x=lon, y=lat),
                  color = "blue", size = 0.6, alpha = 0.4, lineend = "round") +
        ggtitle(title)  
    }
    return(myplot)
  }
  else if(type=="summary_table")
  {
  test=tbl_summary(data$results[,-c(1,2,7)])
  return(test)
  }
}
