################################################################################
# This is a txt file with changes on moveSIM_helper. Changes are indicated in comments at margin
# library(swfscMisc) needed when using this function from script

moveSIM_helper_JFL <- function (sp, env, days, sigma, dest_x, dest_y, mot_x, mot_y,sp_poly,
                            search_radius,optimum,direction,single_rast,mortality) {
  track <- data.frame()
  track[1,1] <- sp@x  # 1st row 1st col is input x coord
  track[1,2] <- sp@y  # 1st row 2nd col is input y coord
  
  # We recognize that morphological characteristics of a species may affect the speed at which
  # they move. Thus we added these parameters to species class and had them affect motivation
  # (preliminary).
  # Idea: Bigger birds can fly further/faster
  #mot_x_new=mot_x+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  #mot_y_new=mot_y+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  
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
  failures=0
  in_box=FALSE
  
  
  for (step in 2:days) { # These are days
    if(single_rast){
      my_rast=env[[1]]
    }
    else{
      my_rast=env[[step]]
    }
    
    lon_candidate<--9999
    lat_candidate<--9999
    
    # Birds search area is a semicircle of radius (bird can't move North)
    
    test=circle.polygon(track[step-1,1], track[step-1,2], search_radius,
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
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps),proj4string=crs(NOAM))                                # change crs(NOAM) by crs(sp_poly)
    my_bool=tryCatch(!is.null(intersect(my_rast,sps)), error=function(e) return(FALSE))
    #print(my_bool)
    if(my_bool){
      my_rast=crop(my_rast,extent(sps))
      my_rast<-mask(my_rast,sps,inverse=FALSE)
    }
    if(dest_x!=999 & dest_y!=999){
      pt=SpatialPoints(cbind(dest_x,dest_y))
      proj4string(pt)=proj4string(NOAM)                                                  # change proj4string(NOAM) by proj4string(sp_poly)
    }
    
    # We are simulating birds that were captured at a study site in Mexico (-99.11, 19.15).
    # We didn't want to force birds there from the start, but if this study site falls
    # within the search area, we want birds to head in that direction.
    if(dest_x==999 & dest_y==999){
      cell_num=which.min(abs(my_rast))
      if (length(which.min(abs(my_rast)))==0){ #Ignore--edge case error handling
        print("Edge Case 1")
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      cell_num=sample(cell_num,1) # There may be ties so we need to sample 1
      
      best_coordinates=xyFromCell(my_rast,cell_num)
      
    }
    
    else if(!is.na(over(pt,sps,fn=NULL)[1]))
    {best_coordinates=c(dest_x,dest_y)
    in_box=TRUE}
    
    else if(in_box==TRUE){
      best_coordinates=c(dest_x,dest_y)
    }
    
    else{
      # If it doesn't fall within, then just take environmental cell
      # within search area that has minimal distance from optimal value
      cell_num=which.min(abs(my_rast))
      if (length(which.min(abs(my_rast)))==0){ #Ignore--edge case error handling
        print("Edge Case 1")
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      cell_num=sample(cell_num,1) # There may be ties so we need to sample 1
      
      best_coordinates=xyFromCell(my_rast,cell_num)
    }
    
    target_x=best_coordinates[1]
    target_y=best_coordinates[2]
    i=1
    while (is.na(extract(my_rast, matrix(c(lon_candidate,lat_candidate),1,2)))) {
      lon_candidate <- track[step-1,1]+ (sigma * rnorm(1)) + (mot_x_new * (target_x - track[step-1,1]))
      lat_candidate <- track[step-1,2]+ (sigma * rnorm(1)) + (mot_y_new * (target_y - track[step-1,2]))
      i=i+1
      # How to select candidate destination, this is as you originally had it.
      if(i>90){ # Avoid infinite loop
        print("Edge Case 2")
        #print("Best Coords")
        #print(best_coordinates)
        #print("Last coords")
        #print(track[step-1,1:2])
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
    }
    pt=SpatialPoints(cbind(lon_candidate,lat_candidate))
    proj4string(pt)=proj4string(NOAM)                                                           # Change proj4string(NOAM) by proj4string(sp_poly)
    
    if(is.na(over(pt,NOAM,fn=NULL)$OBJECTID)){ #Birds can't stop over ocean (they must be over  # Change is.na(over(pt,NOAM,fn=NULL)$OBJECTID) by (is.na(over(pt, sp_poly, fn = NULL))); $OBJECTID is specific to NOAM shapefile
      # North America)
      print("EDGE Case 3")                                                                      # Edge Case 3 means the agent is not over the shapelayer
      #print(pt)
      track[step:days,1]=NA
      track[step:days,2]=NA
      return(track)
      #break
    }
    # Second searching step: now that we've added a destination for this step (and added some)
    # randomness, we want to simulate bird finding best location in close proximity to where it
    # ended up (small scale searching behavior, whereas earlier is larger scale from evolutinary
    # memory.)
    
    neig <- adjacent(my_rast,
                     cellFromXY(my_rast, matrix(c(lon_candidate, #put step in brackets here
                                                  lat_candidate), 1,2)),
                     directions=8, pairs=FALSE )
    # Get cell numbers for adjacent cells
    options <- data.frame() # Create blank dataframe
    for (i in 1:length(neig)){
      options[i,1]<-neig[i]   # ith row first column is each neighboring cell
      options[i,2]<- my_rast[neig[i]]
    }
    
    option <- c(options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ],
                options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ])
    
    if (is.null(option)){ # Ignore--edge case error handling
      track[step:days,1]=NA
      track[step:days,2]=NA
      break
    }
    new_cell <- sample(option,1)
    new_coords <- xyFromCell(my_rast,new_cell) #put step in brackets here
    track[step,1] <- new_coords[1]
    track[step,2] <- new_coords[2]
    
    
    if(is.na(my_rast[new_cell])){
      failures=failures
    }
    
    # Want to simulate bird mortality. If a bird's actual environmental value, as determined
    # by its landing point, is more than 50% off from its optimal, that counts as a "failure"
    # for that day. 5 or more consecutive failures leads to death.
    
    else if(abs(my_rast[new_cell])>.50*optimum){
      failures=failures+1
    }
    else{
      failures=0
    }
    
    if(failures>4 & mortality==TRUE){
      print('Agent died')
      track[(step+1):days,1]=NA #Bird died, rest of points are N/A
      track[(step+1):days,2]=NA
      break
    }
  }    
  return(track)
}

################################################################################
# This would be the "corrected" function (for test purposes)

moveSIM_helper_JFL <- function (sp, env, days, sigma, dest_x, dest_y, mot_x, mot_y,sp_poly,
                                search_radius,optimum,direction,single_rast,mortality) {
  track <- data.frame()
  track[1,1] <- sp@x  # 1st row 1st col is input x coord
  track[1,2] <- sp@y  # 1st row 2nd col is input y coord
  
  # We recognize that morphological characteristics of a species may affect the speed at which
  # they move. Thus we added these parameters to species class and had them affect motivation
  # (preliminary).
  # Idea: Bigger birds can fly further/faster
  #mot_x_new=mot_x+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  #mot_y_new=mot_y+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  
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
  failures=0
  in_box=FALSE
  
  
  for (step in 2:days) { # These are days
    if(single_rast){
      my_rast=env[[1]]
    }
    else{
      my_rast=env[[step]]
    }
    
    lon_candidate<--9999
    lat_candidate<--9999
    
    # Birds search area is a semicircle of radius (bird can't move North)
    
    test=circle.polygon(track[step-1,1], track[step-1,2], search_radius,
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
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps),proj4string=crs(sp_poly))                                # change crs(NOAM) by crs(sp_poly)
    my_bool=tryCatch(!is.null(intersect(my_rast,sps)), error=function(e) return(FALSE))
    #print(my_bool)
    if(my_bool){
      my_rast=crop(my_rast,extent(sps))
      my_rast<-mask(my_rast,sps,inverse=FALSE)
    }
    if(dest_x!=999 & dest_y!=999){
      pt=SpatialPoints(cbind(dest_x,dest_y))
      proj4string(pt)=proj4string(sp_poly)                                                 # change proj4string(NOAM) by proj4string(sp_poly)
    }
    
    # We are simulating birds that were captured at a study site in Mexico (-99.11, 19.15).
    # We didn't want to force birds there from the start, but if this study site falls
    # within the search area, we want birds to head in that direction.
    if(dest_x==999 & dest_y==999){
      cell_num=which.min(abs(my_rast))
      if (length(which.min(abs(my_rast)))==0){ #Ignore--edge case error handling
        print("Edge Case 1")
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      cell_num=sample(cell_num,1) # There may be ties so we need to sample 1
      
      best_coordinates=xyFromCell(my_rast,cell_num)
      
    }
    
    else if(!is.na(over(pt,sps,fn=NULL)[1]))
    {best_coordinates=c(dest_x,dest_y)
    in_box=TRUE}
    
    else if(in_box==TRUE){
      best_coordinates=c(dest_x,dest_y)
    }
    
    else{
      # If it doesn't fall within, then just take environmental cell
      # within search area that has minimal distance from optimal value
      cell_num=which.min(abs(my_rast))
      if (length(which.min(abs(my_rast)))==0){ #Ignore--edge case error handling
        print("Edge Case 1")
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      cell_num=sample(cell_num,1) # There may be ties so we need to sample 1
      
      best_coordinates=xyFromCell(my_rast,cell_num)
    }
    #browser()
    target_x=best_coordinates[1]
    target_y=best_coordinates[2]
    i=1
    while (is.na(extract(my_rast, matrix(c(lon_candidate,lat_candidate),1,2)))) {
      lon_candidate <- track[step-1,1]+ (sigma * rnorm(1)) + (mot_x_new * (target_x - track[step-1,1]))
      lat_candidate <- track[step-1,2]+ (sigma * rnorm(1)) + (mot_y_new * (target_y - track[step-1,2]))
      i=i+1
      # How to select candidate destination, this is as you originally had it.
      if(i>90){ # Avoid infinite loop
        print("Edge Case 2")
        #print("Best Coords")
        #print(best_coordinates)
        #print("Last coords")
        #print(track[step-1,1:2])
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
    }
    pt=SpatialPoints(cbind(lon_candidate,lat_candidate))
    proj4string(pt)=proj4string(sp_poly)                                                           # Change proj4string(NOAM) by proj4string(sp_poly)
    if(is.na(over(pt, sp_poly, fn = NULL))){ #Birds can't stop over ocean (they must be over       # Change is.na(over(pt,NOAM,fn=NULL)$OBJECTID) by (is.na(over(pt, sp_poly, fn = NULL))); $OBJECTID is specific to NOAM shapefile
      # North America)
      print("EDGE Case 3")                                                                         # Edge Case 3 means the agent is not over the shapelayer
      #print(pt)
      track[step:days,1]=NA
      track[step:days,2]=NA
      return(track)
      #break
    }
    # Second searching step: now that we've added a destination for this step (and added some)
    # randomness, we want to simulate bird finding best location in close proximity to where it
    # ended up (small scale searching behavior, whereas earlier is larger scale from evolutinary
    # memory.)
    
    neig <- adjacent(my_rast,
                     cellFromXY(my_rast, matrix(c(lon_candidate, #put step in brackets here
                                                  lat_candidate), 1,2)),
                     directions=8, pairs=FALSE )
    # Get cell numbers for adjacent cells
    options <- data.frame() # Create blank dataframe
    for (i in 1:length(neig)){
      options[i,1]<-neig[i]   # ith row first column is each neighboring cell
      options[i,2]<- my_rast[neig[i]]
    }
    
    option <- c(options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ],
                options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ])
    
    if (is.null(option)){ # Ignore--edge case error handling
      track[step:days,1]=NA
      track[step:days,2]=NA
      break
    }
    new_cell <- sample(option,1)
    new_coords <- xyFromCell(my_rast,new_cell) #put step in brackets here
    track[step,1] <- new_coords[1]
    track[step,2] <- new_coords[2]
    
    
    if(is.na(my_rast[new_cell])){
      failures=failures
    }
    
    # Want to simulate bird mortality. If a bird's actual environmental value, as determined
    # by its landing point, is more than 50% off from its optimal, that counts as a "failure"
    # for that day. 5 or more consecutive failures leads to death.
    
    else if(abs(my_rast[new_cell])>.50*optimum){
      failures=failures+1
    }
    else{
      failures=0
    }
    
    if(failures>4 & mortality==TRUE){
      print('Agent died')
      track[(step+1):days,1]=NA #Bird died, rest of points are N/A
      track[(step+1):days,2]=NA
      break
    }
  }    
  return(track)
}




