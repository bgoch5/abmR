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
poly <- extent2SpatialPolygon(extent(gc))

plot(gc)

testing=moveSIM(replicates=1,days=20,env_rast=gc, search_radius=300,
                sigma=.4, dest_x=90, dest_y=0, mot_x=1, mot_y=1,modeled_species=N_pop,
                my_shapefile=poly,optimum=50,direction="S",write_results=F,single_rast=TRUE, mortality = F)

lines(testing$results[,1:2])
# When using radius too small (>= 50) Error in sample.int(length(x), size, replace, prob) :
# It should be substitute by something like "radious too small!")
# When using a lot of days Edge Case 2 (it's OK)

for (i in 1:10){
  gc <- stack(gc, gc-(i*10))
}

testing=moveSIM(replicates=1,days=10,env_rast=gc, search_radius=400,
                sigma=.4, dest_x=90, dest_y=0, mot_x=1, mot_y=1,modeled_species=N_pop,
                my_shapefile=poly,optimum=20,direction="S",write_results=F,single_rast=F, mortality = F)


# energySIM_helper_JFL #########################################################
# test
enerTest <- energySIM(replicates=5,days=27,env_rast=ndvi_raster, search_radius=400,
              sigma=.1, dest_x=-108.6, dest_y=26.2, mot_x=.9, mot_y=.9,
              modeled_species=ben_test_pop, my_shapefile=NOAM,
              optimum_lo=.6,optimum_hi=.8,init_energy=100,
              direction="S",write_results=FALSE,single_rast=FALSE,mortality = TRUE)












