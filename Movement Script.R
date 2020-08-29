qt(.975,8)
# Last Updated August 22, 2020
# Environmental Data
#https://www.ncdc.noaa.gov/cdr/terrestrial/normalized-difference-vegetation-index

# Load required libraries and shapefiles
require(swfscMisc)
require(rsMove)
require(raster)
require(rgdal)
require(ggplot2)
require(knitr)
require(kableExtra)
require(caret)
require(lattice)
require(igraph)
require(randomForest)
require(e1071)
require(randomForest)
require(rsMove)
require(MODISTools)
require(raster)  
require(dismo)  
require(ceramic)
require(jpeg)
setwd("./Input_Data")
NOAM=readOGR(layer="NOAM",dsn=".")
ndvi_raster=stack("NDVI_2013.gri")
e=extent(ndvi_raster)
ndvi_composite=raster("NDVI_2013_composite.gri") #formed by taking mean of all these
ndvi_composite=crop(ndvi_composite,e)
setwd("~/move-model/Numerical_Results")


test=ndvi_raster-.1

species <- setClass("species", slots=c(x="numeric", y="numeric", opt="numeric",mass="numeric",
                                       wing="numeric"))


# The following function takes an an object of class `species`, an environmental raster, a 
# number of replicates n, theta_x and theta_y (destination coordinates), alpha_x and alpha_y
# (motivation in x and y generation), sp_poly, and current_gen (see explanation before
# generations function)

go <- function (sp, env, n, sigma, theta_x, alpha_x, theta_y, alpha_y,sp_poly,current_gen) { 
  my_env=env-sp@opt
  track <- data.frame()  
  track[1,1] <- sp@x  # 1st row 1st col is input x coord
  track[1,2] <- sp@y  # 1st row 2nd col is input y coord
  
  # We recognize that morphological characteristics of a species may affect the speed at which
  # they move. Thus we added these parameters to species class and had them affect motivation
  # (preliminary).
  # Idea: Bigger birds can fly further/faster
  alpha_x_new=alpha_x+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  alpha_y_new=alpha_y+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  failures=0
  in_box=FALSE
  
  for (step in 2:n) { # These are days
    if(current_gen>1){ # If you are in generation 2 or following, restrict movement based on
    # last generations "confidence interval" (see below)
      
    if(!is.null(intersect(env[[step]],sp_poly))){
    my_rast=crop(env[[step]],extent(sp_poly))
    my_rast<-mask(my_rast,sp_poly,inverse=FALSE)
    }
    }
    
    else{
    my_rast=env[[step]]
    }
    
    lon_candidate<--9999  
    lat_candidate<--9999  
    
    # Birds search area is a semicircle of radius 375 km (bird can't move North)
    
    test=circle.polygon(track[step-1,1], track[step-1,2], 375,
                        units = "km")
    test=data.frame(test)
    test=subset(test,test[,2]<=track[step-1,2])
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps),proj4string=crs(NOAM))
    if(!is.null(intersect(my_rast,sps))){
    my_rast=crop(my_rast,extent(sps))
    my_rast<-mask(my_rast,sps,inverse=FALSE)
    }
    pt=SpatialPoints(cbind(-99.11,19.15))
    proj4string(pt)=proj4string(NOAM)
    
    
    # We are simulating birds that were captured at a study site in Mexico (-99.11, 19.15).
    # We didn't want to force birds there from the start, but if this study site falls
    # within the search area, we want birds to head in that direction.
    
    if(!is.na(over(pt,sps,fn=NULL)[1]))
    {best_coordinates=c(-99.11,19.15)
    in_box=TRUE}
    else if(in_box==TRUE){
    best_coordinates=c(-99.11,19.15)
    }
    else{
    # If it doesn't fall within, then just take environmental cell
    # within search area that has minimal distance from optimal value
    cell_num=which.min(abs(my_rast))
    cell_num=sample(cell_num,1) # There may be ties so we need to sample 1

    if (is.na(cell_num)){ #Ignore--edge case error handling
      track[step:n,1]=NA 
      track[step:n,2]=NA
      break
    }
    best_coordinates=xyFromCell(my_rast,cell_num)
    }
    target_x=best_coordinates[1]
    target_y=best_coordinates[2]
    i=1
    while (is.na(extract(env[[step]], matrix(c(lon_candidate,lat_candidate),1,2)))) {  
      lon_candidate <- track[step-1,1]+ (sigma * rnorm(1)) + (alpha_x_new * (target_x - track[step-1,1]))  
      lat_candidate <- track[step-1,2]+ (sigma * rnorm(1)) + (alpha_y_new * (target_y - track[step-1,2])) 
      i=i+1
      # How to select candidate destination, this is as you originally had it.
      if(i>20){ # Avoid infite loop
        track[step:n,1]=NA 
        track[step:n,2]=NA
        break
      }
    }  
    
    pt=SpatialPoints(cbind(lon_candidate,lat_candidate))
    proj4string(pt)=proj4string(NOAM)
    
    if(is.na(over(pt,NOAM,fn=NULL)$OBJECTID)){ #Birds can't stop over ocean (they must be over
      # North America)
      track[step:n,1]=NA 
      track[step:n,2]=NA
      break
    }
    
    # Second searching step: now that we've added a destination for this step (and added some)
    # randomness, we want to simulate bird finding best location in close proximity to where it
    # ended up (small scale searching behavior, whereas earlier is larger scale from evolutinary
    # memory.)
    
    neig <- adjacent(env[[step]],   
                     cellFromXY(env[[step]], matrix(c(lon_candidate, #put step in brackets here 
                                              lat_candidate), 1,2)),   
                     directions=8, pairs=FALSE )  
    # Get cell numbers for adjacent cells
    options <- data.frame() # Create blank dataframe
    for (i in 1:length(neig)){  
      options[i,1]<-neig[i]   # ith row first column is each neighboring cell
      options[i,2]<- sp@opt - env[[step]][neig[i]] # 2nd col is difference of environmental
      # value and optimal
    }
    
    option <- c(options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ],   
                options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ])  
    
    if (is.null(option)){ # Ignore--edge case error handling
      track[step:n,1]=NA 
      track[step:n,2]=NA
      break
    }
    new_cell <- sample(option,1)  
    new_coords <- xyFromCell(env[[step]],new_cell) #put step in brackets here
    track[step,1] <- new_coords[1]
    track[step,2] <- new_coords[2]  
    if(is.na(env[[step]][new_cell])){
      failures=failures
    }
    
    # Want to simulate bird mortality. If a bird's actual environmental value, as determined
    # by its landing point, is more than 50% off from its optimal, that counts as a "failure"
    # for that day. 5 or more consecutive failures leads to death.
    
    else if(abs(env[[step]][new_cell]-sp@opt)>.50*sp@opt){ 
      failures=failures+1
    }
    else{
      failures=0
    }
    
    if(failures>4){
      print('Bird died')
      track[step:n,1]=NA #Bird died, rest of points are N/A
      track[step:n,2]=NA 
      break
    }
  }  
  return(track)  
}  


# Define 3 populations: N_pop and S_pop refer to two different populations of Orange Crowned
# Warbler, and WIWA_pop is Wilson's Warbler
N_pop <- species(x=-112.24, y =52.24, opt= .55,mass=7.5,wing=15.5) #change mass and wing based on species
S_pop <- species(x=-106.61,y=34.23,opt=.55,mass=7.5,wing=15.5) 
WIWA_pop= species(x=-131.74,y=64.87,opt=.55,mass=19,wing=9) 
# Their origin locations were determined through stable isotope analysis


# We want to run multiple generations of birds to see if their migration evolves. Thus, we
# run the generations function, which essentially runs the 'go' function multiple times,
# creating a confidence interval polygon each time that will restrict the movements of future 
# generations. All the data and polygons from each generation are returned by the function.

generations=function(n_generations=3, replicates=200,N_pop_days=18,S_pop_days=18,WIWA_pop_days=18,ndvi_rast=ndvi_raster,
                     my_N_pop=N_pop,my_S_pop=S_pop,my_WIWA_pop=WIWA_pop,my_NOAM=NOAM)
  {
  my_N_list=list()
  my_S_list=list()
  my_WIWA_list=list()
  
  my_N_shapefiles=list()
  my_S_shapefiles=list()
  my_WIWA_shapefiles=list()
  
  N_pop_results=array(rep(0,replicates*N_pop_days*2),c(replicates,N_pop_days,2))
  S_pop_results=array(rep(0,replicates*S_pop_days*2),c(replicates,S_pop_days,2))
  WIWA_pop_results=array(rep(0,replicates*WIWA_pop_days*2),c(replicates,WIWA_pop_days,2))
  
  
  for (j in 1:n_generations)
  {
    print(paste0("Starting Generation ", j))
    
    if(j==1){
      sp_poly_N=my_NOAM
      sp_poly_S=my_NOAM
      sp_poly_WIWA=my_NOAM
    }
    
  for(i in 1:replicates){
    
          
    North=data.matrix(go(my_N_pop,ndvi_rast,N_pop_days,.6, -99.11,1,19.15,1,sp_poly_N,j))
    South=data.matrix(go(my_S_pop,ndvi_rast,S_pop_days,.6,-99.11,1,19.15,1,sp_poly_S,j))
    WIWA=data.matrix(go(my_WIWA_pop,ndvi_rast,WIWA_pop_days,.6,-99.11,1,19.15,1,sp_poly_WIWA,j))
    if (length(North)==N_pop_days*2){
    N_pop_results[i,,]=North
    }
    else{
    my_matrix=matrix(rep(NA,N_pop_days*2),nrow=N_pop_days,ncol=2,byrow=TRUE)
    N_pop_results[i,,]=my_matrix
    }
    if (length(South)==S_pop_days*2){
    S_pop_results[i,,]=South
    }
    else{
    my_matrix=matrix(rep(NA,S_pop_days*2),nrow=S_pop_days,ncol=2,byrow=TRUE)
    S_pop_results[i,,]=my_matrix
    }
     if (length(WIWA)==WIWA_pop_days*2){
       WIWA_pop_results[i,,]=WIWA
     }
    else{
      my_matrix=matrix(rep(NA,WIWA_pop_days*2),nrow=WIWA_pop_days,ncol=2,byrow=TRUE)
      WIWA_pop_results[i,,]=my_matrix
    }
    if (i%%5 == 0) {
      print(paste0("Number of birds processed per population this gen.: ",i))
    }
  }
    
  my_N_list[[j]]=N_pop_results
  my_S_list[[j]]=S_pop_results
  my_WIWA_list[[j]]=WIWA_pop_results
  
  aggregate_N=apply(N_pop_results, c(2,3), mean,na.rm=T)
  SD_N=apply(N_pop_results, c(2,3), sd,na.rm=T)
  lower_bound=aggregate_N-1.64*SD_N  #95 was too wide try 90%
  upper_bound=aggregate_N+1.64*SD_N

  my_df=rbind(lower_bound,upper_bound)
  my_df=na.omit(my_df)
  ch <- chull(my_df)
  coords <- my_df[c(ch, ch[1]), ]

  sp_poly_N <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))

  aggregate_S=apply(S_pop_results, c(2,3), mean,na.rm=T)
  SD_S=apply(S_pop_results, c(2,3), sd,na.rm=T)
  lower_bound=aggregate_S-1.64*SD_S  #95 was too wide try 90%
  upper_bound=aggregate_S+1.64*SD_S

  my_df=rbind(lower_bound,upper_bound)
  my_df=na.omit(my_df)
  my_df
  ch <- chull(my_df)
  coords <- my_df[c(ch, ch[1]),]
  sp_poly_S <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
  
  
  aggregate_WIWA=apply(WIWA_pop_results, c(2,3), mean,na.rm=T)
  SD_WIWA=apply(WIWA_pop_results, c(2,3), sd,na.rm=T)
  lower_bound=aggregate_WIWA-1.64*SD_WIWA  #95 was too wide try 90%
  upper_bound=aggregate_WIWA+1.64*SD_WIWA
  
  my_df=rbind(lower_bound,upper_bound)
  my_df=na.omit(my_df)
  my_df
  ch <- chull(my_df)
  coords <- my_df[c(ch, ch[1]),]
  sp_poly_WIWA <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
  

  my_N_shapefiles[[j]]=sp_poly_N
  my_S_shapefiles[[j]]=sp_poly_S
  my_WIWA_shapefiles[[j]]=sp_poly_WIWA

  prelim_results=list(N_pop_final=N_pop_results,S_pop_final=S_pop_results,
                      N_pop_full=my_N_list,S_pop_full=my_S_list,shapefiles_N=my_N_shapefiles,
                      shapefiles_S=my_S_shapefiles,
                      WIWA_final=WIWA_pop_results,WIWA_full=my_WIWA_list,shapefiles_WIWA=my_WIWA_shapefiles)
  save(prelim_results,file=paste0("PrelimResults_Gen",j,"_",Sys.Date(),".RData"))
  }

  results=list(N_pop_final=N_pop_results,S_pop_final=S_pop_results,
               N_pop_full=my_N_list,S_pop_full=my_S_list,shapefiles_N=my_N_shapefiles,
               shapefiles_S=my_S_shapefiles,
               WIWA_final=WIWA_pop_results,WIWA_full=my_WIWA_list,shapefiles_WIWA=my_WIWA_shapefiles)
  return(results)

}


# Run the function
# Will have to bump this up later, this is just to see if code is working for you
my_results=generations(n_generations=3,replicates=10)
my_results
save(my_results,file="July 22 Results.RData")
proc.time()-ptm

WIWA_results=my_results$WIWA_final
N_pop_results=my_results$N_pop_final
S_pop_results=my_results$S_pop_final


aggregate_WIWA=apply(WIWA_results, c(2,3), mean,na.rm=T)
aggregate_N=apply(N_pop_results, c(2,3), mean,na.rm=T)
aggregate_S=apply(S_pop_results, c(2,3), mean,na.rm=T)

getwd()
setwd("~/move-model/Plots")


NEW_NOAM=crop(NOAM,extent(ndvi_composite))
extent(NEW_NOAM)
par(mar = c(2.5, 2.5, 2.5, 2.5))

jpeg(filename = "Plot of Aggregated Tracks new.jpeg",
     pointsize = 12, width=8, height=8, units="in", res=1000)
plot(ndvi_composite,xlim=c(-150,-60),ylim=c(10,70),
     xlab="Longitude",ylab="Latitude",main="Migration Tracks and NDVI for 3
     Migrant Populations")
legend(-144,24, legend=c("WIWA", "N. OCWA","S. OCWA"),
       col=c("light blue", "black","black"),lty=c(1,1,3),cex=1.25,lwd=c(3,3,2),bg="white",box.lwd=2.5)
lines(aggregate_N,lwd=2.2,col="black")
lines(aggregate_S,lwd=3,col="black",lty="dotted")
lines(aggregate_WIWA,lwd=2.2,col="light blue")
points(-112.24,52.24)
points(-106.61,34.23)
points(-131.74,64.87)
points(-99.11,19.15,pch=19,col="black")
plot(NEW_NOAM,add=TRUE)
dev.off()

