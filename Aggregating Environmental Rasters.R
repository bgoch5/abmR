# Constructing environmental Raster stacks to work with abmr
# NDVI Example

# Benjamin Gochanour
# Last Updated: April 14, 2021


# I. Download data
# NOAA AVHRR NDVI Data can be downloaded from
# https://www.ncei.noaa.gov/data/avhrr-land-normalized-difference-vegetation-index/access/

# Here use dates
# 8/26/2013-9/21/2013.

# II. Compile Data Together

library(raster)
library(ncdf4)
setwd("C:\\Users\\BGOCHANOUR\\Documents\\NDVI_Data")

files=list.files(path=".")
my_list=list()

files # Print list of files, make sure these are in order by date
      # so they stack in correct order. If they are not, re-order them.

for (file in files){
  layer=raster(file,varname="NDVI")
  my_list=append(my_list,layer)
}
my_list # Dates look in order here as well
my_stack=stack(my_list)

plot(my_stack[[1]]) # Examine first layer

# III. Crop to regions of interest and write results

# Create European datafiles
e=extent(c(-15,63,0,60))
EU=crop(my_stack,e)
EU_composite=mean(EU)

writeRaster(EU,"NDVI_2013_Europe.gri")
writeRaster(EU_composite,"NDVI_2013_Europe_composite.gri")

# Create North America datafiles
e2=extent(c(-160,-40,10,90))
NorthAmerica=crop(my_stack,e2)
NA_composite=mean(NorthAmerica)

writeRaster(NorthAmercia,"NDVI_2013_NA.gri")
writeRaster(NA_composite,"NDVI_2013_NA_composite.gri")