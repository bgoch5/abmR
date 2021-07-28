#'
#' Run model for one replicate
#'
#' Runs agent based modeling for one replicate of a single species. Used by function `generations`
#' to run more replicates of more species.
#' @import raster sp rgdal 
#' @importFrom swfscMisc circle.polygon
#' @param sp A species object
#' @param env Raster, should represent NDVI or your environmental variable of interest
#' @param days Integer, how many days (timesteps), would you like to model
#' @param sigma Numeric, amount of random error
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction
#' @param mot_y Numeric, movement motivation in y direction
#' @param current_gen Fed into function by function genSIM; if using this function
#' alone use 1
#' @param search_radius Radius of semicircle to South of current location to search for next timestep (in km)
#' @param mortality Logical, should low energy levels result in death? Default T.
#' @return A nx2 dataset containing longitude and latitude points for all n timesteps
#' @keywords internal
#' @export

moveSIM_helper <- function(sp, env, days, sigma, dest_x, dest_y, mot_x, mot_y,
                           search_radius, optimum, n_failures, fail_thresh, direction, single_rast, mortality) {
  track <- data.frame()
  track[1, 1] <- sp@x
  track[1, 2] <- sp@y
  track[1:days, 3] <- "Alive"
  track[1:days, 4] <- "Alive"

  # We recognize that morphological characteristics of a species may affect the speed at which
  # they move. Thus we added these parameters to species class and had them affect motivation.

  if (length(sp@morphpar1 == 1) & length(sp@morphpar2 == 1)) {
    mot_x_new <- (mot_x + (sp@morphpar1 - sp@morphpar1mean) / sp@morphpar1sd * .1 * ifelse(sp@morphpar1sign == "Pos", 1, -1)
      + (sp@morphpar2 - sp@morphpar2mean) / sp@morphpar2sd * .1 * ifelse(sp@morphpar2sign == "Pos", 1, -1))

    mot_y_new <- (mot_y + (sp@morphpar1 - sp@morphpar1mean) / sp@morphpar1sd * .1 * ifelse(sp@morphpar1sign == "Pos", 1, -1)
      + (sp@morphpar2 - sp@morphpar2mean) / sp@morphpar2sd * .1 * ifelse(sp@morphpar2sign == "Pos", 1, -1))
  }
  else {
    mot_x_new <- mot_x
    mot_y_new <- mot_y
  }

  failures <- 0
  in_box <- FALSE

  for (step in 2:days) {
    if (single_rast) {
      my_rast <- env[[1]]
    }
    else {
      my_rast <- env[[(step - 1)]]
    }

    lon_candidate <- -9999
    lat_candidate <- -9999

    # Birds search area is a semicircle of radius
    test <- swfscMisc::circle.polygon(track[step - 1, 1], track[step - 1, 2], search_radius,
      units = "km"
    )
    test <- data.frame(test)
    if (direction == "S") {
      test <- subset(test, test[, 2] <= track[step - 1, 2])
    }
    else if (direction == "N") {
      test <- subset(test, test[, 2] >= track[step - 1, 2])
    }
    else if (direction == "E") {
      test <- subset(test, test[, 1] >= track[step - 1, 1])
    }
    else if (direction == "W") {
      test <- subset(test, test[, 1] <= track[step - 1, 1])
    }
    else if (direction == "R") {
      test <- test
    }
    p <- Polygon(test)
    ps <- Polygons(list(p), 1)
    sps <- SpatialPolygons(list(ps),
      proj4string = crs(env)
    )
    my_bool <- tryCatch(!is.null(intersect(my_rast, sps)), error = function(e) {
      return(FALSE)
    })
    # print(my_bool)
    if (my_bool) {
      my_rast <- crop(my_rast, extent(sps))
      my_rast <- mask(my_rast, sps, inverse = FALSE)
    }
    if (dest_x != 999 & dest_y != 999) {
      pt <- SpatialPoints(cbind(dest_x, dest_y))
      proj4string(pt) <- proj4string(env)
    }

    if (direction == "R") {
      random <- sampleRandom(my_rast, 1, xy = TRUE)
      track[step, 1] <- random[1]
      track[step, 2] <- random[2]
    }
    else {

      # We are simulating birds that were captured at a study site in Mexico (-99.11, 19.15).
      # We didn't want to force birds there from the start, but if this study site falls
      # within the search area, we want birds to head in that direction.
      if (dest_x == 999 & dest_y == 999) {
        cell_num <- which.min(abs(my_rast))
        if (length(which.min(abs(my_rast))) == 0) {
          print("Can't find any non-NA cells. Agent stopped.")
          track[step:days, 1] <- NA
          track[step:days, 2] <- NA
          track[step:days, 3] <- "Stopped"
          track[1:days, 4] <- "Stopped"
          return(track)
        }
        cell_num <- sample(cell_num, 1) # There may be ties so we need to sample 1

        best_coordinates <- xyFromCell(my_rast, cell_num)
      }

      else if (!is.na(over(pt, sps, fn = NULL)[1])) {
        best_coordinates <- c(dest_x, dest_y)
        in_box <- TRUE
      }

      else if (in_box == TRUE) {
        best_coordinates <- c(dest_x, dest_y)
      }

      else {
        # If it doesn't fall within, then just take environmental cell
        # within search area that has minimal distance from optimal value
        cell_num <- which.min(abs(my_rast))
        if (length(which.min(abs(my_rast))) == 0) { # Ignore--edge case error handling
          print("Can't find any non-NA cells. Agent stopped.")
          track[step:days, 1] <- NA
          track[step:days, 2] <- NA
          track[step:days, 3] <- "Stopped"
          track[1:days, 4] <- "Stopped"
          break
        }
        cell_num <- sample(cell_num, 1) # There may be ties so we need to sample 1
        best_coordinates <- xyFromCell(my_rast, cell_num)
      }

      target_x <- best_coordinates[1]
      target_y <- best_coordinates[2]
      i <- 1
      while (is.na(extract(my_rast, matrix(c(lon_candidate, lat_candidate), 1, 2)))) {
        lon_candidate <- track[step - 1, 1] + (sigma * rnorm(1)) + (mot_x_new * (target_x - track[step - 1, 1]))
        lat_candidate <- track[step - 1, 2] + (sigma * rnorm(1)) + (mot_y_new * (target_y - track[step - 1, 2]))
        i <- i + 1
        if (i > 90) { # Avoid infinite loop
          print("Can't find any non-NA cells. Agent stopped.")
          track[step:days, 1] <- NA
          track[step:days, 2] <- NA
          track[step:days, 3] <- "Stopped"
          track[1:days, 4] <- "Stopped"
          return(track)
        }
      }
      pt <- SpatialPoints(cbind(lon_candidate, lat_candidate))
      proj4string(pt) <- proj4string(env)

      if (is.na(over(pt, sps, fn = NULL))) { # Birds can't stop over ocean (they must be over
        # North America)
        print("Best coordinates not in search region, agent stopped")
        track[step:days, 1] <- NA
        track[step:days, 2] <- NA
        track[step:days, 3] <- "Stopped"
        track[1:days, 4] <- "Stopped"
        return(track)
        # break
      }
      # Second searching step: now that we've added a destination for this step (and added some)
      # randomness, we want to simulate bird finding best location in close proximity to where it
      # ended up (small scale searching behavior, whereas earlier is larger scale from evolutinary
      # memory.)

      neig <- adjacent(my_rast,
        cellFromXY(my_rast, matrix(c(
          lon_candidate, # put step in brackets here
          lat_candidate
        ), 1, 2)),
        directions = 8, pairs = FALSE
      )
      # Get cell numbers for adjacent cells
      options <- data.frame() # Create blank dataframe
      for (i in 1:length(neig)) {
        options[i, 1] <- neig[i] # ith row first column is each neighboring cell
        options[i, 2] <- my_rast[neig[i]]
      }

      option <- c(
        options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1],
        options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1]
      )

      if (is.null(option) | length(option) == 0) { # Ignore--edge case error handling
        print("Can't find any non-NA cells. Agent stopped.")
        track[step:days, 1] <- NA
        track[step:days, 2] <- NA
        track[step:days, 3] <- "Stopped"
        track[1:days, 4] <- "Stopped"
        return(track)
      }
      new_cell <- sample(option, 1)
      new_coords <- xyFromCell(my_rast, new_cell) # put step in brackets here
      track[step, 1] <- new_coords[1]
      track[step, 2] <- new_coords[2]


      if (is.na(my_rast[new_cell])) {
        failures <- failures
      }

      # Want to simulate bird mortality. If a bird's actual environmental value, as determined
      # by its landing point, is more than 50% off from its optimal, that counts as a "failure"
      # for that day. 5 or more consecutive failures leads to death.

      else if (abs(my_rast[new_cell]) > fail_thresh * optimum) {
        failures <- failures + 1
      }
      else {
        failures <- 0
      }

      if (failures > n_failures & mortality == TRUE) {
        print("Agent died")
        track[(step + 1):days, 1] <- NA # Bird died, rest of points are N/A
        track[(step + 1):days, 2] <- NA
        track[(step + 1):days, 3] <- "Died"
        track[1:days, 4] <- "Died"
        return(track)
      }
    }
  }
  return(track)
}
