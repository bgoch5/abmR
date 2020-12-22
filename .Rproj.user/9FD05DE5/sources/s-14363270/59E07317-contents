# Converting
library(raster)
library(ncdf4)
setwd("C:\\Users\\BGOCHANOUR\\Documents\\NDVI_Data")
r <- raster('AVHRR-Land_v005_AVH13C1_NOAA-19_20130826_c20170408052048.nc', varname="NDVI")
plot(r)

files=list.files(path=".")
my_list=list()

files
for (file in files){
  layer=raster(file,varname="NDVI")
  my_list=append(my_list,layer)
}

my_list

my_stack=stack(my_list)
plot(my_stack[[1]])
e=extent(c(-15,63,0,60))
extent(my_stack)
my_stack_test=crop(my_stack,e)
plot(my_stack_test[[1]])
ndvi_composite=mean(my_stack_test)

writeRaster(my_stack_test,"NDVI_2013_Europe.gri")
writeRaster(ndvi_composite,"NDVI_2013_Europe_composite.gri")

plot(my_stack[[1]])

class(layer)
extent(my_stack)
plot(layer)
