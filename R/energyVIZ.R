
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
#' @import table1
#'
#'
#' @param data Data to be plotted, this object should be the output from
#' energySIM().
#' @param title Title for the plot that is output.
#' @param type String from "plot", "gradient", or "table"?
#' @param aspect_ratio Aspect ratio, defaults to 4/3
#' @param label Logical, label the origin and specified final destination?
#' @param xlim Optionally specify desired x limits as a numeric vector: c(low,hi)
#' @param ylim Optionally specify desired y limits as a numeric vector: c(low,hi)
#'
#' @return
#' A plot showing model output compared to a dashed line that represents straight line
#' movement from the starting point to the final destination.
#' @examples
#' testing=energySIM(replicates=2,days=27,env_rast=ndvi_raster, search_radius=375,
#' sigma=.4, dest_x=-100, dest_y=20, mot_x=1, mot_y=1, 
#' modeled_species=my_species, my_shapefile=NOAM, optimum_lo=.4,optimum_hi=.6,init_energy=100,
#' direction="S",write_results=TRUE,single_rast=FALSE)
#' energyVIZ(testing,title="Visualizing EnergySIM results",type="plot", aspect_ratio=5/3,
#' label=TRUE)
#' @export

energyVIZ=function(data,title="Visualizing EnergySIM results",type="plot",
                   aspect_ratio=4/3, label=FALSE,
                   xlim=NULL,ylim=NULL)
{
if(type=="plot"){
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

#Latitude
if(my_ylim[1]>ideal.df[2,2]){
  my_ylim[1]=ideal.df[2,2]-4
}

if(my_ylim[2]<ideal.df[2,2]){
  my_ylim[2]=ideal.df[2,2]+4
}

# Longitude
if(my_xlim[1]>ideal.df[2,1]){
  my_xlim[1]=ideal.df[2,1]-4
}

if(my_xlim[2]<ideal.df[2,1]){
  my_xlim[2]=ideal.df[2,1]+4
}

y_diff=abs(my_ylim[2]-my_ylim[1])
x_diff=abs(my_xlim[2]-my_xlim[1])

if(y_diff>2*x_diff){
  my_xlim[1]=my_xlim[1]-.5*x_diff
  my_xlim[2]=my_xlim[2]+.5*x_diff
}

if(x_diff>2*y_diff){
  my_ylim[1]=my_ylim[1]-.5*y_diff
  my_ylim[2]=my_ylim[2]+.5*y_diff
}

if(!is.null(xlim)){
  my_xlim=xlim
}

if(!is.null(ylim)){
  my_ylim=ylim
}

myplot=ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim =my_xlim,
           ylim = my_ylim, 
           expand = FALSE) +
  geom_path(data = t.energy.res,
            aes(x=lon, y=lat),
            color = "red", size = 0.6, alpha = 0.4, lineend = "round") +
  geom_path(data = ideal.df,
            aes(x=Lon, y=Lat),
            color = "black", size = 1.2, alpha = 1, linetype = "dashed") + theme(aspect.ratio=aspect_ratio)
  ggtitle(title)
if(label){
  ideal.df[,"type"]=NA
  ideal.df[1,4]="Origin"
  ideal.df[2,4]="Ideal Final"
  myplot=myplot+geom_point(data=ideal.df,aes(x=Lon,y=Lat,color=type))
}
return(myplot)
}
if(type=="gradient")
{}
if(type=="table")
{
t.energy.res <- data$results
table1(~energy + day + distance | delta_energy, data = t.energy.res)}
}
