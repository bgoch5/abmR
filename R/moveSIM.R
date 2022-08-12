#' Runs agent-based model (ABM) movement simulations based on environmental data
#'
#' Here, agent mortality occurs when agent fails to achieve suitable raster
#' values at least n_failures+1 timesteps in a row. Agent energy stores are not
#' dynamic, so movement speed isn't directly affected by quality of raster cells
#' achieved. Results may be analyzed with moveVIZ(). Relies on underlying
#' function moveSIM_helper(), which is not to be used alone.
#'
#' For each timestep, agents can have status "Alive",
#' "Stopped", or "Died". All agents start alive and may stop if, on a particular timestep,
#' there are no non-NA raster values in the search region. This often occurs when agents
#' are searching over an ocean or a large lake, for example. Once an agent stops, they
#' remain stopped for the rest of the run. Similarly, once an agent dies, they retain
#' this status for all subsequent timesteps. All timesteps with agent status "Stopped"
#' or "Died" will have lat/lon=NA, so as to not affect subsequent analyses.
#'
#' Arguments mortality, n_failures, and fail_thresh interact with each other. If
#' mortality = F, values for n_failures and fail_thresh are ignored. If mortality=T, fail_thresh
#' determines what constitutes a failure, and n_failures indicates how many failures are allowed before
#' death. Note: If n_failures=days, this is equivalent to mortality=F.
#'
#' @import raster sp rgdal
#' @importFrom methods as  setClass
#' @importFrom stats na.omit rbinom rnorm
#' @importFrom geosphere distHaversine
#' @importFrom utils write.csv
#' @param replicates Integer, desired number of replicates per run. Default 100.
#' @param days Integer, how many days (timesteps) would you like to model? Range (1,nlayers(env_rast))
#' @param env_rast Rasterstack or Rasterbrick with number of layers >= days
#' @param search_radius Radius of semicircle search regions (in km). Default 375.
#' @param sigma Numeric, randomness parameter, range (-Inf, Inf). Default 0.1.
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction, range (0,1], default 1.
#' @param mot_y Numeric, movement motivation in y direction, range (0,1], default 1.
#' @param modeled_species Object of class "species"
#' @param optimum Numeric, optimal environmental value
#' @param fail_thresh What percentage deviation from optimum leads to death? E.g. default of
#' .50 means 50 percent or greater deviation from optimum on a particular step constitutes failure.
#' @param n_failures How many failures are allowable before agent experiences death (at n_failures+1). What constitutes
#' a failure is determined by fail_thresh, range (1,days]. Default 4.
#' @param direction Character, movement direction, one of "N","S","E","W", or "R" (Random). Default "S".
#' @param mortality Logical, should low energy levels result in death? Default T.
#' @param write_results Logical, save results to csv? Default F.
#' @param single_rast Logical, are you using a one-layer raster for all timesteps? Default F.
#'
#' @return
#' Under "results", a (days+1 * replicates) row X 7 column dataframe containing data on agent_id, day, longitude, latitude,
#' current agent status (Alive, Stopped, or Died), distance traveled from last timestep (in km), and final status.
#' Using tidy_results() provides a cleaner display of results.
#'
#' Under "run_params", a record of function parameters used as well as missing_pct
#' and mortality_pct. missing_pct corresponds to the percent of rows in the results dataframe
#' missing information on lon/lat, which occurs when the agent has "died" or "stopped". mortality_pct
#' refers to the percentage of agents in the run that died.
#'
#' @examples
#' # Define species object
#' pop1 <- as.species(
#'   x = -98.7, y = 34.7)
#'  
#' # Run function
#' EX2 <- moveSIM( replicates = 3, days = 10, env_rast = ex_raster,
#' search_radius = 300, sigma = .1, dest_x = -108.6, dest_y = 26.2,
#'  mot_x = .8, mot_y = .8, modeled_species = pop1, optimum = .6,
#'   n_failures = 5, fail_thresh = .40, direction = "R",
#'   write_results = FALSE, single_rast = TRUE, mortality = TRUE)
#' # View Results in Clean Format
#' tidy_results(EX2, type = "results")
#' tidy_results(EX2, type = "run_params")
#' 
#' @export

moveSIM <- function(replicates = 100,
                    days,
                    modeled_species,
                    env_rast,
                    optimum,
                    dest_x,
                    dest_y,
                    mot_x,
                    mot_y,
                    search_radius = 375,
                    direction = "S",
                    sigma = 0.1,
                    mortality = TRUE,
                    fail_thresh = .5,
                    n_failures = 4,
                    single_rast = FALSE,
                    write_results = FALSE) {
  days <- days + 1
  sp <- modeled_species
  if (nlayers(env_rast) == 1 & single_rast == FALSE) {
    stop("Single layer environmental raster with single_rast=FALSE specified.
         Please check your raster or change to single_rast=TRUE. Exiting
         function")
  }
  if (nlayers(env_rast) != 1 & single_rast == TRUE) {
    warning("Multiple layer environmental raster with single_rast=TRUE specified.
    Using only first layer of raster")
  }

  my_min <- minValue(env_rast[[1]])
  my_max <- maxValue(env_rast[[1]])
  if (!(optimum > my_min | optimum < my_max)) {
    warning("Optimum specified is not in range of env_raster values. Consider changing.")
  }

  if (n_failures < 1 | n_failures > days) {
    warning("Value for n_failures is not in interval (1,days]. Please check your work.")
  }

  if (sp@x < xmin(env_rast) | sp@x > xmax(env_rast)
  | sp@y < ymin(env_rast) | sp@y > ymax(env_rast)) {
    stop("Species origin point outside env raster extent")
  }
  
  if (mot_x <0 | mot_y<0) {
    stop("mot_x and mot_y must be in range (0, 1]")
  }

  if (direction != "R") {
    my_env <- env_rast - optimum
  }
  else {
    my_env <- env_rast
    message("Direction=R specified--Raster will be ignored")
  }

  long <- data.frame(
    lon = numeric(), lat = numeric(), day = numeric(),
    agent_id = character()
  )

  for (i in 1:replicates) {
    Species <- moveSIM_helper(
      sp = modeled_species, env = my_env, days = days, sigma = sigma,
      dest_x = dest_x, dest_y = dest_y, mot_x = mot_x, mot_y = mot_y, search_radius = search_radius, optimum = optimum,
      n_failures = n_failures, fail_thresh = fail_thresh, direction = direction, mortality = mortality,
      single_rast = single_rast
    )
    names(Species) <- c("lon", "lat", "curr_status", "final_status")
    Species$day <- 0:(nrow(Species) - 1)
    number <- sprintf("%02d", i)
    Species$agent_id <- paste("Agent", number, sep = "_")
    Species$distance <- NA
    for (j in 2:nrow(Species)) {
      Species$distance[j] <- distHaversine(Species[(j - 1), 1:2], Species[j, 1:2]) / 1000
    }
    if (nrow(Species) == days) {
      long <- rbind(long, Species)
    }
    if (i %% 5 == 0) {
      message(paste0("Number of agents processed: ", i))
    }
  }

  col_order <- c("agent_id", "day", "lon", "lat", "curr_status", "distance", "final_status")
  long <- long[, col_order]

  if (write_results) {
    currentDate <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
    file_name <- paste("moveSIM_results_", currentDate, ".csv", sep = "")
    write.csv(long, file_name)
  }

  missing_pct <- sum(is.na(long$lon)) / nrow(long) * 100
  mortality_pct <- round(sum(long[, "final_status"] == "Died") / 28, 0) / replicates * 100

  params <- data.frame(
    replicates = replicates, days = (days - 1),
    env_raster = deparse(substitute(env_raster)),
    search_radius = search_radius, sigma = sigma, dest_x = dest_x,
    dest_y = dest_y, mot_x = mot_x, mot_y = mot_y,
    modeled_species = deparse(substitute(modeled_species)),
    optimum = optimum, n_failures = n_failures, direction = direction, mortality = mortality, write_results = write_results,
    single_rast = single_rast, missing_pct = missing_pct, mortality_pct = mortality_pct
  )
  return(list(results = long, run_params = params))
}
