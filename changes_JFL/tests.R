######## ABMr TESTs #########

# moveSIM_helper_JFL ###########################################################

# You should change moveSIM to use moveSIM_helper_JFL instead of moveSIM_helper

library (raster)  
library(dismo)  
library(oceanmap)
library(abmr)

N_pop=species(x=10,y=80)

gc <- raster(nrows = 100, ncols = 100, xmn = 0, xmx = 100, ymn = 0, ymx = 100)  
gc[] <- rep(1:100, 100)
poly <- extent2SpatialPolygon(gc)

testing=moveSIM(replicates=1,days=20,env_rast=gc, search_radius=500,
                sigma=.4, dest_x=90, dest_y=0, mot_x=1, mot_y=1,modeled_species=N_pop,
                my_shapefile=poly,optimum=100,direction="S",write_results=TRUE,single_rast=TRUE)