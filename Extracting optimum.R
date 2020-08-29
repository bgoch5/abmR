#install.packages("ebirdst")
library(ebirdst)
library(rgdal)
library(raster)

##################### Can also do occurence, here do for Wilson's Warbler ###################
test=ebirdst_runs
getwd()
WIWA_path <- ebirdst_download(species = "wlswar")# wlswar, orcwar
OCWA_path <- ebirdst_download(species = "orcwar")

WIWA_occur <- load_raster("occurrence", path = WIWA_path)
date_vector <- parse_raster_dates(WIWA_occur)


rappdirs::user_data_dir("ebirdst")

# August 26th thru September 20th?
WIWA_1=WIWA_occur[[34]]
WIWA_2=WIWA_occur[[35]]
WIWA_3=WIWA_occur[[36]]
WIWA_4=WIWA_occur[[37]]


setwd("./Input_Data")
NOAM=readOGR(layer="NOAM",dsn=".")
ndvi_raster=stack("NDVI_2013.gri")
e=extent(ndvi_raster)
ndvi_composite=raster("NDVI_2013_composite.gri")
ndvi_composite=crop(ndvi_composite,e)

my_crs=crs(ndvi_composite)

memory.limit(size=100000)

proj_WIWA_1=projectRaster(WIWA_1,crs=my_crs,method="ngb")

plot(proj_WIWA_1)
plot(NOAM,add=TRUE)

proj_WIWA_2=projectRaster(WIWA_2,crs=my_crs,method="ngb")
proj_WIWA_3=projectRaster(WIWA_3,crs=my_crs,method="ngb")
proj_WIWA_4=projectRaster(WIWA_4,crs=my_crs,method="ngb")

proj_WIWA_1=crop(proj_WIWA_1,extent(ndvi_composite))
proj_WIWA_2=crop(proj_WIWA_2,extent(ndvi_composite))
proj_WIWA_3=crop(proj_WIWA_3,extent(ndvi_composite))
proj_WIWA_4=crop(proj_WIWA_4,extent(ndvi_composite))


plot(proj_WIWA_1)
plot(NOAM,add=TRUE)

plot(proj_WIWA_2)
plot(NOAM,add=TRUE)

plot(proj_WIWA_3)
plot(NOAM,add=TRUE)

plot(proj_WIWA_4)
plot(NOAM,add=TRUE)
extent(proj_WIWA_1)

plot(proj_WIWA_1)
plot(proj_WIWA_2)
plot(proj_WIWA_3)
plot(proj_WIWA_4)

# there's an overlay function
# r3 = overlay(r1, r2, fun=function(x,y){return(x+y)})