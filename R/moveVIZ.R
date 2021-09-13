#' Creates a plot/table of moveSIM() results
#'
#' When type="plot", function plots the movement tracks versus the the straight
#' line track between the origin and destination (unless the destination was
#' unspecified in the call to moveSIM(), then straight line track is omitted).
#' When type="summary_table", a summary table is output.
#'
#' @import raster sp rgdal rnaturalearth rnaturalearthdata ggplot2
#' @importFrom gtsummary tbl_summary
#' @param data Data to be plotted, this object should be the output from
#' moveSIM().
#' @param type "plot" or "summary_table", default "plot".
#' @param title Title for the plot that is output.
#' @param aspect_ratio Aspect ratio, defaults to 1.
#' @param xlim Optionally specify desired x limits as a numeric vector: c(low,hi)
#' @param ylim Optionally specify desired y limits as a numeric vector: c(low,hi)
#'
#' @examples
#' 
#' # 1. Define Population and Run moveSIM()
#' 
#' pop1 <- as.species(x=-100, y=55,
#' morphpar1=15, morphpar1mean=16, morphpar1sd=2,morphpar1sign="Pos",
#' morphpar2=19,morphpar2mean=18,morphpar2sd=1,morphpar2sign="Pos")
#' 
#' EX2=moveSIM(replicates=5,days=27,env_rast=ex_raster, search_radius=550,
#' sigma=.1, dest_x=-108.6, dest_y=26.2, mot_x=.8, mot_y=.8,
#' modeled_species=pop1,optimum=.6, n_failures=5, fail_thresh=.40,
#' direction="S",write_results=FALSE,single_rast=TRUE,mortality = TRUE)
#' 
#' # 2. Run moveVIZ() on your result
#' moveVIZ(EX2,title="Visualizing MoveSIM results",type="plot")
#' 
#' moveVIZ(EX2, type="summary_table")
#'
#' @export

moveVIZ=function(data, type="plot", title="moveSIM results", aspect_ratio=1, xlim=NULL, ylim=NULL)
{
  if(type=="plot"){
    dest_x=data$run_params$dest_x
    dest_y=data$run_params$dest_y
    world <- ne_countries(scale = "medium", returnclass = "sf")
    start.p <- cbind(data$results[1,3], data$results[1,4])
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
    t.move.res <- data$results
    my_xlim = c((min(t.move.res$lon,na.rm=T)-6), (max(t.move.res$lon,na.rm=T)+6))
    my_ylim = c((min(t.move.res$lat,na.rm=T)-6), (max(t.move.res$lat,na.rm=T)+6))
    
    if(!is.null(xlim)){
      my_xlim=xlim
    }
    
    if(!is.null(ylim)){
      my_ylim=ylim
    }
    
    if(dest_x!=999 & dest_y!=999){
    myplot=ggplot(data = world) +
      geom_sf() +
      coord_sf(xlim = my_xlim, ylim = my_ylim, expand = FALSE) +
      geom_path(data = t.move.res,
                aes(x=t.move.res$lon, y=t.move.res$lat,group=t.move.res$agent_id),
                color = "blue", size = 0.6, alpha = 0.4, lineend = "round") +
      geom_path(data = ideal.df,
                aes(x=ideal.df$Lon, y=ideal.df$Lat),
                color = "black", size = 1.2, alpha = 1, linetype = 2) + theme(aspect.ratio=aspect_ratio) + 
      ggtitle(title)
    }
    else{
      myplot=ggplot(data = world) +
        geom_sf() +
        coord_sf(xlim = my_xlim, ylim = my_ylim, expand = FALSE) +
        geom_path(data = t.move.res,
                  aes(x=t.move.res$lon, y=t.move.res$lat,group=t.move.res$agent_id),
                  color = "blue", size = 0.6, alpha = 0.4, lineend = "round") + theme(aspect.ratio=aspect_ratio) + 
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
