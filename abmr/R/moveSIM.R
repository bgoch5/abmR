#' Run model for one replicate
#'
#' Runs agent based modeling for one replicate of a single species. Used by function `generations`
#' to run more replicates of more species.
#'
#' @import raster
#' @import sp
#' @import rgdal
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
#' @param current_gen Fed into function by function genSIM; if using this function
#' alone use 1
#' @param search_radius Radius of semicircle to South of current location to search for next timestep (in km)
#'
#' @return A nx2 dataset containing longitude and latitude points for all n timesteps
#' @examples
#' my_results = moveSIM(sp = wiwa.pop, env = ndvi_raster, n = 27, sigma = 0.6,
#' dest_x = -100, dest_y = 25, mot_x = 0.9, mot_y = 0.9, search_radius = 200, current_gen = 1)
#' @export

moveSIM <- function (sp, env, days, sigma, dest_x, dest_y, mot_x, mot_y,sp_poly,current_gen,
                     search_radius) {
  track <- data.frame()
  track[1,1] <- sp@x  # 1st row 1st col is input x coord
  track[1,2] <- sp@y  # 1st row 2nd col is input y coord

  # We recognize that morphological characteristics of a species may affect the speed at which
  # they move. Thus we added these parameters to species class and had them affect motivation
  # (preliminary).
  # Idea: Bigger birds can fly further/faster
  mot_x_new=mot_x+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  mot_y_new=mot_y+(sp@mass-7.5)*.01+(sp@wing-15.5)*.01
  failures=0
  in_box=FALSE

  for (step in 2:days) { # These are days
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

    # Birds search area is a semicircle of radius (bird can't move North)

    test=circle.polygon(track[step-1,1], track[step-1,2], search_radius,
                        units = "km")
    test=data.frame(test)
    test=subset(test,test[,2]<=track[step-1,2])
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps),proj4string=crs(NOAM))
    my_bool=tryCatch(!is.null(intersect(my_rast,sps)), error=function(e) return(FALSE))
    if(my_bool){
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
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
      best_coordinates=xyFromCell(my_rast,cell_num)
    }
    target_x=best_coordinates[1]
    target_y=best_coordinates[2]
    i=1
    while (is.na(extract(env[[step]], matrix(c(lon_candidate,lat_candidate),1,2)))) {
      lon_candidate <- track[step-1,1]+ (sigma * rnorm(1)) + (mot_x_new * (target_x - track[step-1,1]))
      lat_candidate <- track[step-1,2]+ (sigma * rnorm(1)) + (mot_y_new * (target_y - track[step-1,2]))
      i=i+1
      # How to select candidate destination, this is as you originally had it.
      if(i>20){ # Avoid infinite loop
        track[step:days,1]=NA
        track[step:days,2]=NA
        break
      }
    }

    pt=SpatialPoints(cbind(lon_candidate,lat_candidate))
    proj4string(pt)=proj4string(NOAM)

    if(is.na(over(pt,NOAM,fn=NULL)$OBJECTID)){ #Birds can't stop over ocean (they must be over
      # North America)
      track[step:days,1]=NA
      track[step:days,2]=NA
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
      track[step:days,1]=NA
      track[step:days,2]=NA
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
      track[step:days,1]=NA #Bird died, rest of points are N/A
      track[step:days,2]=NA
      break
    }
  }
  return(track)
}
