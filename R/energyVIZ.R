
#' Creates a plot/table of energySIM() results
#'
#' When type="plot", function plots the movement tracks versus the the straight
#' line track between the origin and destination (unless the destination was
#' unspecified in the call to energySIM(), then straight line track is omitted).
#' When type="gradient", creates a gradient plot showing what regions cause
#' agents to gain/lose energy. Two table  options are also available using
#' type="summary_table" or type="strat_table" (table with results stratified
#' by energy gain or loss). Please see Vignette for examples of this output.
#'
#' @import raster sp rgdal rnaturalearth rnaturalearthdata ggplot2 table1 gstat sf tmap
#' @importFrom gtsummary tbl_summary
#' @param data Data to be plotted, this object should be the output from
#' energySIM().
#' @param type String from "plot", "gradient", "summary_table", or "strat_table"?
#' @param title Title for the plot that is output.
#' @param aspect_ratio Aspect ratio, defaults to 1. 
#' @param label Logical, label the origin and specified final destination?
#' @param xlim Optionally specify desired x limits as a numeric vector: c(low,hi)
#' @param ylim Optionally specify desired y limits as a numeric vector: c(low,hi)
#' @examples
#' 
#' # 1. Define Population and Run energySIM()
#' 
#' pop1 <- as.species(x=-98.7, y=34.7,
#' morphpar1=15, morphpar1mean=16, morphpar1sd=2,morphpar1sign="Pos",
#' morphpar2=19,morphpar2mean=18,morphpar2sd=1,morphpar2sign="Pos")
#' 
#' EX1=energySIM(replicates=2,days=27,env_rast=ex_raster, search_radius=200,
#' sigma=.1, dest_x=-108.6, dest_y=26.2, mot_x=.9, mot_y=.9,
#' modeled_species=pop1,
#' optimum_lo=.8,optimum_hi=.9,init_energy=100,
#' direction="S",write_results=FALSE,single_rast=TRUE,mortality = TRUE)
#' 
#' # 2. Run energyVIZ() on your result
#' 
#' energyVIZ(EX1,title="Visualizing EnergySIM results",type="plot", aspect_ratio=5/3,
#' label=TRUE)
#' 
#' energyVIZ(EX1,type="summary_table")
#' 
#' energyVIZ(EX1,type="strat_table")
#' 
#' \dontrun{energyVIZ(EX1,type="gradient")}
#' @export

energyVIZ=function(data, type="plot", title="energySIM results",
                   aspect_ratio=1, label=FALSE,
                   xlim=NULL,ylim=NULL)
{
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
t.energy.res <- data$results

my_xlim = c((min(t.energy.res$lon,na.rm=T)-6), (max(t.energy.res$lon,na.rm=T)+6))
my_ylim = c((min(t.energy.res$lat,na.rm=T)-6), (max(t.energy.res$lat,na.rm=T)+6))

if (dest_x!=999){
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
}
if(!is.null(xlim)){
  my_xlim=xlim
}

if(!is.null(ylim)){
  my_ylim=ylim
}

if(type=="plot"){
  
if(dest_x!=999 & dest_y!=999){
myplot=ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim =my_xlim,
           ylim = my_ylim, 
           expand = FALSE) +
  geom_path(data = t.energy.res,
            aes(x=t.energy.res$lon, y=t.energy.res$lat,group=t.energy.res$agent_id),
            color = "red", size = 0.6, alpha = 0.4, lineend = "round") +
  geom_path(data = ideal.df,
            aes(x=ideal.df$Lon, y=ideal.df$Lat),
            color = "black", size = 1.2, alpha = 1, linetype = "dashed") + theme(aspect.ratio=aspect_ratio) + 
  ggtitle(title)
}
else{
  myplot=ggplot(data = world) +
    geom_sf() +
   coord_sf(xlim =my_xlim,
             ylim = my_ylim, 
             expand = FALSE) +
    geom_path(data = t.energy.res,
              aes(x=t.energy.res$lon, y=t.energy.res$lat,group=t.energy.res$agent_id),
              color = "red", size = 0.6, alpha = 0.4, lineend = "round") + theme(aspect.ratio=aspect_ratio) +
    ggtitle(title) 
  label=FALSE
}
if(label){
  ideal.df[,"type"]=NA
  ideal.df[1,4]="Origin"
  ideal.df[2,4]="Ideal Final"
  myplot=myplot+geom_point(data=ideal.df,aes(x=ideal.df$Lon,y=ideal.df$Lat,color=type))
}
return(myplot)
}
if(type=="gradient")
{my.df = data$results
  my.sf.point = my.df
  my_vector=!is.na(my.df$x)
  my.sf.point <- my.sf.point[my_vector,]
  my.sf.point$energy = NULL
  my.sf.point$day = NULL
  my.sf.point$agent_id = NULL
  my.sf.point$distance = NULL
  my.df$X = NULL
  my.df=na.omit(my.df)
  #------------------------------------------------
  # First Energy interpolation (Basic POINT Plot)
  my.sf.point <- st_as_sf(x = my.df,
                          coords = c("lon", "lat"),
                          crs = "+proj=longlat +datum=WGS84")
  my.sp.point <- as(my.sf.point, "Spatial")
  ###############################################
  # reducing the extent of that huge NOAM shp
  # over a target area (doing this manually for now)
  # but maybe the function will need to use the input raster?
  #----------------------------------------------
  world.redu <- st_crop(world, extent(c(my_xlim[1],my_xlim[2],my_ylim[1],my_ylim[2])), snap="out")
  #=============================================
  # plotting energy +/-
  grd <- as.data.frame(spsample(my.sp.point, "regular", n=50000))
  names(grd) <- c("X", "Y")
  coordinates(grd) <- c("X", "Y")
  gridded(grd) <- TRUE
  fullgrid(grd) <- TRUE
  proj4string(grd) <- proj4string(my.sp.point)
  P.idw <- gstat::idw(delta_energy ~ 1, my.sp.point, newdata=grd, idp=2.0)
  r <- raster(P.idw)
  r.m <- mask(r, world.redu)
  my_plot=tm_shape(r.m) +
    tm_raster(n=10,palette = "RdBu", midpoint = NA,
              title="energy") +
    tm_shape(my.sp.point) + tm_dots(size=0.01, alpha=0.1) +
    tm_legend(legend.outside=T)
  return(my_plot)
  }
if(type=="summary_table"){
  test=tbl_summary(data$results[,-c(1,2,9)])
  return(test)
}
if(type=="strat_table")
{
t.energy.res <- data$results
table1(~energy + day + distance | delta_energy, data = t.energy.res)}
}
