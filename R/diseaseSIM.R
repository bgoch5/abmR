#'
#' Runs more advanced Brownian / Ornstein Uhlenbeck agent-based model for
#' multiple replicates, with disease features.
#'
#' Here, agent mortality occurs when agent reaches energy = 0. Agent energy
#' stores are dynamic, and affect search area as a multiplier, so movement
#' is directly affected by the quality of raster cells achieved. Results
#' may be visualized with `energyVIZ()`. Relies on underlying function
#' `energySIM_helper`, which is not to be used alone.
#'
#' For each timestep, agents can have status "Alive",
#' "Stopped", or "Died". All agents start alive and may stop if, on a particular
#' timestep, there are no non-NA raster values in the search region. This often
#' occurs when agents are searching over an ocean or a large lake, for example.
#' Once an agent stops, they remain stopped for the rest of the run. Similarly,
#' once an agent dies, they retain this status for all subsequent timesteps.
#' All timesteps with agent status "Stopped" or "Died" will have lat/lon=NA,
#' so as to not affect subsequent analyses.
#' 
#' @import raster sp rgdal
#' @param replicates Integer, desired number of replicates per run, default 100.
#' @param days Integer, How many days (timesteps) would you like to model? Range (1,nlayers(env_rast))
#' @param env_rast Rasterstack or Rasterbrick with number of layers >= days
#' @param search_radius Radius of semicircle search regions (in km). Default 375.
#' @param sigma Numeric, randomness parameter, range (-Inf, Inf). Default 0.1.
#' @param dest_x Numeric, destination x coordinate (longitude)
#' @param dest_y Numeric, destination y coordinate (latitude)
#' @param mot_x Numeric, movement motivation in x direction, range (0,1], default 1.
#' @param mot_y Numeric, movement motivation in y direction, range (0,1], default 1.
#' @param modeled_species Object of class "species"
#' @param optimum_lo Numeric, optimal environmental value (low)
#' @param optimum_hi Numeric, optimal environmental value (high)
#' @param init_energy Numeric, initial energy in interval (0,100]
#' @param direction Character, movement direction, one of "N","S","E","W", or "R" (Random). Default "S".
#' @param mortality Logical, should low energy levels result in death? Default T.
#' @param energy_adj Numeric, Vector of length 11 representing desired energy gain/penalty corresponding to achieved env values
#' in optimum range (1st element), and within 10, 20,  ..., 80, 90, and 90+ percent (11th element) of the average of optimum hi and optimum lo.
#' Recommend using default which is decreasing and symmetric about zero but can modify if desired.
#' @param write_results Logical, save results to csv? Default F.
#' @param single_rast Logical, are you using a one-layer raster for all timesteps?. Default F.
#' @param disease_loc Dataframe of x,y coordinates specifying the location of disease clusters.
#' @param disease_radius Numeric (0,Infty), what is the maximum distance from the source point is capable
#' capable of causing disease?
#' @param disease_mortality Numeric (0,1] What proportion of agents who pass directly
#' over disease source (maximum disease load) will experience mortality? Assume mortality rate then linearly
#' decreases to 0 at disease_radius.
#' @param disease_energy_interact Numeric (0,100] Below what energy level should all agents
#' exposed to disease die? 
#' @return
#' Under "results", a (days+1 X replicates) rows X 9 column dataframe containing data on agent_id, day, longitude, latitude,
#' current agent status (Alive, Stopped, or Died), energy, change in energy from last time_step,
#' distance traveled from last timestep (in km), and final status.
#' Using tidy_results() provides a cleaner display of results.
#'
#' Under "run_params", a record of function parameters used as well as missing_pct
#' and mortality_pct. missing_pct corresponds to the percent of rows in the results dataframe
#' missing information on lon/lat, which occurs when the agent has "died" or "stopped". mortality_pct
#' refers to the percentage of agents in the run that died.
#'
#' @examples
#' # Define species object
#' pop1 <- as.species(x=-100, y=55,
#' morphpar1=15, morphpar1mean=16, morphpar1sd=2,morphpar1sign="Pos",
#' morphpar2=19,morphpar2mean=18,morphpar2sd=1,morphpar2sign="Pos")
#' 
#' # Define disease locations
#' x=c(-99.475, -98.700, -102.725,-108.325,-98.700,-94.625,-95.175,-100.425, -103.675
#' ,-102.475, -100.925, -100.025, -101.425,  -99.125,  -99.275)
#' y=c(32.225, 34.700, 34.575, 34.425, 34.700, 34.375, 32.525, 32.175, 31.675, 29.575, 28.225, 
#'  26.275, 26.075, 25.025, 24.825)
#'  disease_points=data.frame(x=x,y=y)
#'
#' # Run function
#' EX1 <- diseaseSIM(replicates=15,days=27,env_rast=ndvi_raster, 
#' search_radius=400,
#' sigma=.1, dest_x=999, dest_y=999, mot_x=.9, mot_y=.9,
#' modeled_species=pabu.pop.new, optimum_lo=.6,optimum_hi=.8,init_energy=100, 
#' direction="S",write_results=TRUE,single_rast=FALSE,mortality = TRUE,
#' energy_adj=c(30,25,20,5,0,-5,-5,-10,-20,-25,-30),disease_loc=disease_points,
#' disease_energy_interact = 60, disease_mortality=.5,disease_radius=300)
#'
#' # View Results in Clean Format
#' tidy_results(EX1, type = "results")
#' tidy_results(EX1, type = "run_params")
#' @export

diseaseSIM <- function(replicates = 100,
                      days,
                      modeled_species,
                      env_rast,
                      optimum_lo,
                      optimum_hi,
                      dest_x,
                      dest_y,
                      mot_x,
                      mot_y,
                      search_radius = 375,
                      direction = "S",
                      sigma = .1,
                      mortality = TRUE,
                      init_energy = 100,
                      energy_adj = c(25, 20, 15, 10, 5, 0, -5, -10, -15, -20, -25),
                      single_rast = FALSE,
                      write_results = FALSE,
                      disease_loc,
                      disease_radius=5,
                      disease_mortality=.5,
                      disease_energy_interact=40) {
  days <- days + 1
  sp <- modeled_species
  if (optimum_hi < optimum_lo) {
    print("Error: opt_hi smaller than opt_lo--please check your work")
    stop()
  }
  
  if (init_energy < 0 | init_energy > 100) {
    print("Error: Initial Energy should be such that 0<init_energy<=100")
    stop()
  }
  
  if (length(energy_adj) != 11) {
    cat("Error: Supplied energy_adj vector not of required length 11. Please check your work
        and consult documentation for more info on energy_adj")
    stop()
  }
  
  if (nlayers(env_rast) == 1 & single_rast == FALSE) {
    cat("Error: Single layer environmental raster with single_rast=FALSE specified.
         Please check your raster or change to single_rast=TRUE. Exiting
         function")
    stop()
  }
  
  if (nlayers(env_rast) != 1 & single_rast == TRUE) {
    cat("Warning: Multiple layer environmental raster with single_rast=TRUE specified.
    Using only first layer of raster")
  }
  
  if (length(sp@morphpar1) == 1 & length(sp@morphpar2) == 1) {
    if (sp@morphpar1 > sp@morphpar1mean + 3.5 * sp@morphpar1sd | sp@morphpar1 < sp@morphpar1mean - 3.5 * sp@morphpar1sd) {
      cat("Error: Specified value for morphpar1 is greater than 3.5 SDs away from the
         specified mean, which is extremely unlikely. Consider adjusting morphpar1, morphpar1mean,
         or morphpar1SD. Function terminated.")
      stop()
    }
    
    if (sp@morphpar2 > sp@morphpar2mean + 3.5 * sp@morphpar2sd | sp@morphpar2 < sp@morphpar2mean - 3.5 * sp@morphpar2sd) {
      cat("Error: Specified value for morphpar2 is greater than 3.5 SDs away from the
         specified mean, which is extremely unlikely. Consider adjusting morphpar2, morphpar2mean,
         or morphpar2SD. Function terminated.")
      stop()
    }
  }
  
  my_min <- minValue(env_rast[[1]])
  my_max <- maxValue(env_rast[[1]])
  if (!(my_min <= optimum_hi && my_max >= optimum_lo)) {
    cat("Warning: Optimum range specified does not overlap with range of env_raster values.
       Consider changing.")
  }
  
  
  if (modeled_species@x < xmin(env_rast) | modeled_species@x > xmax(env_rast)
      | modeled_species@y < ymin(env_rast) | modeled_species@y > ymax(env_rast)) {
    print("Error: Species origin point outside env raster extent")
    stop()
  }
  
  if (!is.data.frame(disease_loc)) {
    cat("Error: Supplied input for disease_loc is not a data frame. Please convert to
        this input format.")
    stop()
  }
  
  if (disease_mortality<0 | disease_mortality>1) {
    cat("Error: Supplied disease_mortality not in interval (0,1]. Please check your work.")
    stop()
  }
  
  if (disease_radius <= 0) {
    cat("Error: Supplied disease_radius not in interval (0,Infty). Please check your work.")
    stop()
  }
  
  if (disease_energy_interact<0 | disease_energy_interact>100) {
    cat("Error: Supplied disease_energy_interact not in interval (0,100]. Please check your work.")
    stop()
  }

  optimum <- (optimum_hi + optimum_lo) / 2
  
  if (direction != "R") {
    my_env <- env_rast - optimum
  }
  else {
    my_env <- env_rast
    print("Direction=R specified--Raster will be ignored")
  }
  
  long <- data.frame(
    lon = numeric(), lat = numeric(), energy = numeric(),
    curr_status = character(), final_status = character()
  )
  
  for (i in 1:replicates) {
    Species <- diseaseSIM_helper(
      sp = modeled_species,
      env_orig = env_rast,
      env_subtract = my_env,
      days = days,
      sigma = sigma,
      dest_x = dest_x,
      dest_y = dest_y,
      mot_x = mot_x,
      mot_y = mot_y,
      search_radius = search_radius,
      optimum_lo = optimum_lo,
      optimum_hi = optimum_hi,
      init_energy = init_energy,
      direction = direction,
      mortality = mortality,
      energy_adj = energy_adj,
      single_rast = single_rast,
      disease_loc=disease_loc,
      disease_radius=disease_radius,
      disease_mortality=disease_mortality,
      disease_energy_interact=disease_energy_interact
    )
    names(Species) <- c("lon", "lat", "energy", "curr_status", "final_status")
    Species$day <- 0:(nrow(Species) - 1)
    number <- sprintf("%02d", i)
    Species$agent_id <- paste("Agent", number, sep = "_")
    Species$distance <- NA
    Species$delta_energy <- NA
    for (j in 2:nrow(Species)) {
      Species$distance[j] <- distHaversine(Species[(j - 1), 1:2], Species[j, 1:2]) / 1000
      Species$delta_energy[j] <- Species[j, 3] - Species[(j - 1), 3]
    }
    if (nrow(Species) == days) {
      long <- rbind(long, Species)
    }
    if (i %% 5 == 0) {
      print(paste0("Number of agents processed: ", i))
    }
  }
  
  
  col_order <- c(
    "agent_id", "day", "lon", "lat", "curr_status", "energy",
    "delta_energy", "distance", "final_status"
  )
  long <- long[, col_order]
  
  if (write_results) {
    currentDate <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
    file_name <- paste("energySIM_results_", currentDate, ".csv", sep = "")
    write.csv(long, file_name)
  }
  
  missing_pct <- sum(is.na(long$lon)) / nrow(long) * 100
  mortality_pct <- length(which(long$energy == 0)) / replicates * 100
  
  params <- data.frame(
    replicates = replicates, days = (days - 1),
    env_raster = deparse(substitute(env_raster)),
    search_radius = search_radius, sigma = sigma, dest_x = dest_x,
    dest_y = dest_y, mot_x = mot_x, mot_y = mot_y,
    modeled_species = deparse(substitute(modeled_species)),
    optimum_lo = optimum_lo, optimum_hi = optimum_hi, init_energy = init_energy,
    direction = direction, mortality = mortality, write_results = write_results,
    single_rast = single_rast, missing_pct = missing_pct,
    mortality_pct = mortality_pct, disease_loc=deparse(substitute(disease_loc)),disease_radius=disease_radius,
    disease_mortality=disease_mortality,disease_energy_interact=disease_energy_interact
  )
  return(list(results = long, run_params = params))
}
