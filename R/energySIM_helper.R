#' Run energy-dynamic based model for one replicate
#'
#' Runs agent based modeling for one replicate of a single species.
#'
#' @import raster
#' @import sp
#' @import rgdal
#' @import geosphere
#'
#' @param sp A species object
#' @param env Raster, should represent NDVI or your environmental variable of interest
#' @param days Integer, how many days (timesteps), would you like to model
#' @param sigma Numeric, amount of random error
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction
#' @param mot_y Numeric, movement motivation in y direction
#' @param sp_poly Come back to this
#' @param search_radius Radius of semicircle to South of current location to search for next timestep (in km)
#' @param optimum_lo come back
#' @param optimum_hi come back
#' @param init_energy come back
#' @param direction come back
#' @param mortality Logical, should low energy levels result in death?
#' @return A nx3 dataset containing longitude and latitude and energy
#' points for all n timesteps
#' @examples
#' my_results = moveSIM(sp = wiwa.pop, env = ndvi_raster, n = 27, sigma = 0.6,
#' dest_x = -100, dest_y = 25, mot_x = 0.9, mot_y = 0.9, search_radius = 200, current_gen = 1)
#' @export

energySIM_helper <- function (sp, env_orig,env_subtract, days, sigma, dest_x, dest_y, mot_x, mot_y,sp_poly,
                     search_radius,optimum_lo,optimum_hi,init_energy,direction,single_rast,mortality)
  {
  track <- data.frame()
  track[1,1] <- sp@x  
  track[1,2] <- sp@y  
  track[1,3]=init_energy

  in_interval=FALSE

  optimum=(optimum_hi+optimum_lo)/2

  # We recognize that morphological characteristics of a species may affect the speed at which
  # they move. Thus we added these parameters to species class and had them affect motivation
  # (preliminary).
  # Idea: Bigger birds can fly further/faster
  
  if(length(sp@morphpar1==1) & length(sp@morphpar2==1)){
  mot_x_new=(mot_x+(sp@morphpar1-sp@morphpar1mean)/sp@morphpar1sd*.1*ifelse(sp@morphpar1sign == "Pos", 1, -1)
             +(sp@morphpar2-sp@morphpar2mean)/sp@morphpar2sd*.1*ifelse(sp@morphpar2sign == "Pos", 1, -1))

  mot_y_new=(mot_y+(sp@morphpar1-sp@morphpar1mean)/sp@morphpar1sd*.1*ifelse(sp@morphpar1sign == "Pos", 1, -1)
             +(sp@morphpar2-sp@morphpar2mean)/sp@morphpar2sd*.1*ifelse(sp@morphpar2sign == "Pos", 1, -1))
  }
  else{
    mot_x_new=mot_x
    mot_y_new=mot_y
  }
  if(mot_x_new<0){
    mot_x_new=.001
    mot_y_new=.001
  }

  in_box=FALSE

  energy=init_energy

  for (step in 2:days) {

    if(single_rast){
    curr_env_subtract=env_subtract[[1]]
    curr_env_orig=env_orig[[1]]
    }
    else{
    curr_env_subtract=env_subtract[[step]]
    curr_env_orig=env_orig[[step]]
    }

    lon_candidate<--9999
    lat_candidate<--9999

    # Birds search area is a semicircle of search_radius in direction specified
    if(mortality){
    search_radius_update=search_radius*(energy/100)
    }
    
    else{
    search_radius_update=search_radius  # Mortality also turns off search radius multiplier
    }
    
    test=circle.polygon(track[step-1,1], track[step-1,2], search_radius_update,
                        units = "km")
    test=data.frame(test)
    if(direction=="S")
    {test=subset(test,test[,2]<=track[step-1,2])}
    else if (direction=="N")
    {test=subset(test,test[,2]>=track[step-1,2])}
    else if (direction=="E")
    {test=subset(test,test[,1]>=track[step-1,1])}
    else if (direction=="W")
    {test=subset(test,test[,1]<=track[step-1,1])}
    else if (direction=="R"){
    test=test}
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps),proj4string=crs(sp_poly))
    my_bool=tryCatch(!is.null(intersect(curr_env_subtract,sps)), error=function(e) return(FALSE))

    if(my_bool){
      curr_env_subtract=crop(curr_env_subtract,extent(sps))
      curr_env_subtract<-mask(curr_env_subtract,sps,inverse=FALSE)
      curr_env_orig=crop(curr_env_orig,extent(sps))
      curr_env_orig<-mask(curr_env_orig,sps,inverse=FALSE)
    }
    
    if(direction=="R"){
      random=sampleRandom(curr_env_orig,1,xy=TRUE)
      track[step,1] <- random[1]
      track[step,2] <- random[2]
    }
    
    else{
      
    if(dest_x!=999 & dest_y!=999){
    pt=SpatialPoints(cbind(dest_x,dest_y))
    proj4string(pt)=proj4string(sp_poly)
    }

    # We are simulating birds that were captured at a study site in Mexico (-99.11, 19.15).
    # We didn't want to force birds there from the start, but if this study site falls
    # within the search area, we want birds to head in that direction.
    if(dest_x==999 & dest_y==999){
      cell_num=which.min(abs(curr_env_subtract))
      if (length(which.min(abs(curr_env_subtract)))==0){ #Ignore--edge case error handling
        print("Edge Case 1")
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      cell_num=sample(cell_num,1) # There may be ties so we need to sample 1
      best_coordinates=xyFromCell(curr_env_subtract,cell_num)
    }
    
    else if(!is.na(over(pt,sps,fn=NULL)[1]))
    {best_coordinates=c(dest_x,dest_y)
    in_box=TRUE}
  
    else if(in_box==TRUE & dest_x!=999 & dest_y!=999){
      best_coordinates=c(dest_x,dest_y)
    }
    else{
      # If it doesn't fall within, then just take environmental cell
      # within search area that has minimal distance from optimal value
      cell_num=which.min(abs(curr_env_subtract)) # had my_rast here, need curr_env_subtract
      if (length(which.min(abs(curr_env_subtract)))==0){ #Ignore--edge case error handling
        print("Edge Case 1")
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      cell_num=sample(cell_num,1) # There may be ties so we need to sample 1
      best_coordinates=xyFromCell(curr_env_subtract,cell_num)
    }
    target_x=best_coordinates[1]
    target_y=best_coordinates[2]
    i=1
    while(is.na(extract(curr_env_subtract, matrix(c(lon_candidate,lat_candidate),1,2)))) {
      lon_candidate <- track[step-1,1]+ (sigma * rnorm(1)) + (mot_x_new * (target_x - track[step-1,1]))
      lat_candidate <- track[step-1,2]+ (sigma * rnorm(1)) + (mot_y_new * (target_y - track[step-1,2]))
      pt=SpatialPoints(cbind(lon_candidate,lat_candidate))
      proj4string(pt)=proj4string(sp_poly)
      i=i+1
      # How to select candidate destination, this is as you originally had it.
      if(i>35){ # Avoid infinite loop
        print("Edge Case 2")
        track[step:days,1]=NA
        track[step:days,2]=NA
        return(track)
      }
    }

    pt=SpatialPoints(cbind(lon_candidate,lat_candidate))
    proj4string(pt)=proj4string(sp_poly)

    if(is.na(over(pt,sp_poly,fn=NULL)$OBJECTID)){ #Birds can't stop over ocean (they must be over
      # North America)
      print("Edge Case 3")
      track[step:days,1]=NA
      track[step:days,2]=NA
      break
    }

    # Second searching step: now that we've added a destination for this step (and added some)
    # randomness, we want to simulate bird finding best location in close proximity to where it
    # ended up (small scale searching behavior, whereas earlier is larger scale from evolutinary
    # memory.)

    neig <- adjacent(curr_env_subtract,
                     cellFromXY(curr_env_subtract, matrix(c(lon_candidate, #put step in brackets here
                    lat_candidate), 1,2)), directions=8, pairs=FALSE )
    # Get cell numbers for adjacent cells
    options <- data.frame() # Create blank dataframe
    for (i in 1:length(neig)){
      options[i,1]<-neig[i]   # ith row first column is each neighboring cell
      options[i,2]<- curr_env_subtract[neig[i]] # 2nd col is difference of environmental
      # value and optimal
    }
    option <- c(options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ],
                options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ])

    if (is.null(option) | length(option)==0){ # Ignore--edge case error handling
      print("New Edge Case")
      track[step:days,1]=NA
      track[step:days,2]=NA
      break
    }
    
    new_cell <- sample(option,1)

    if(length(option==8)){
      new_cell=cellFromXY(curr_env_subtract, matrix(c(lon_candidate, #put step in brackets here
                                                lat_candidate), 1,2))
    }
    # If everything in the neighborhood is NA use the cell itself

    new_coords <- xyFromCell(curr_env_subtract,new_cell) #put step in brackets here
    track[step,1] <- new_coords[1]
    track[step,2] <- new_coords[2]

    dist_from_opt=curr_env_orig[new_cell]-optimum
    dist_from_opt_hi=curr_env_orig[new_cell]-optimum_hi
    dist_from_opt_lo=curr_env_orig[new_cell]-optimum_lo
    opts=list(dist_from_opt,dist_from_opt_hi,dist_from_opt_lo)
    abs_opts=list(abs(dist_from_opt),abs(dist_from_opt_hi),abs(dist_from_opt_lo))
    my_min=which.min(abs_opts)

    if(my_min==1 | (my_min==2 & dist_from_opt_hi<0)|(my_min==3 & dist_from_opt_lo>0))
    {in_interval=TRUE}
    else if (my_min==2){
    diff=abs(dist_from_opt_hi)
    }
    else if (my_min==3){
    diff=abs(dist_from_opt_lo)
    }
    else{
    diff=.21*optimum #this is just a fix for now when all cells in neighborhood have
    # NA values; will have to fix later
    }

    if(in_interval){
        energy=energy+15
    }
    else if (diff<.1*optimum){
      energy=energy+10
    }
    else if (diff<.2*optimum){
      energy=energy+5
    }
    else if (diff<.3*optimum){
      energy=energy
    }
    else if (diff<.4*optimum){
      energy=energy-5
    }
    else if (diff<.5*optimum){
      energy=energy-10
    }
    else if (diff<.6*optimum){
      energy=energy-15
    }
    else{
      energy=energy-20
    }

    if(energy>100)
    {energy=100}
    if(energy<0)
    {energy=0}

    track[step,3]=energy


    in_interval=FALSE
    if(mortality==TRUE){
    if(energy==0 & step<days){
      print('Agent died')
      track[(step+1):days,1]=NA #Bird died, rest of points are N/A
      track[(step+1):days,2]=NA
      break
    }
    }
    }
    }
  return(track)
}
