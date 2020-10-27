
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
    myplot=ggplot(data = world) +
      geom_sf() +
      coord_sf(xlim = c(-127, -78), ylim = c(10, 55), expand = FALSE) +
      geom_path(data = t.energy.res,
                aes(x=lon, y=lat),
                color = "blue", size = 0.6, alpha = 0.4, lineend = "round") +
      geom_path(data = ideal.df,
                aes(x=Lon, y=Lat),
                color = "black", size = 1.2, alpha = 1, linetype = 2) +
      ggtitle(title)
    return(myplot)
}
