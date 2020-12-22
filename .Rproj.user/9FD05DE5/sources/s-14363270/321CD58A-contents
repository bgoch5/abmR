#'
#' Imports example data
#'
#' @import googledrive
#' @import purrr
#' @import rgdal
#' @import raster
#' @examples
#' get_ex_data()
#' @export

get_ex_data=function(){
  folder_url <- "https://drive.google.com/drive/folders/1hV1kmtw4u8fh1crAsxI7zBaySqgjTZaJ"
  folder <- drive_get(as_id(folder_url))
  my_files <- drive_ls(folder)
  my_wd=getwd()  # Record working directory for later use
  
  # Download files, will save to hard drive
  walk(my_files$id, ~ drive_download(as_id(.x),overwrite=TRUE))
  
  # Read in Europe Raster (composite is average of all days for plotting)
  ndvi_raster_Europe=stack(paste0(my_wd,"/NDVI_2013_Europe.gri"))
  ndvi_raster_Europe_composite=stack(paste0(my_wd,"/NDVI_2013_Europe_composite.gri"))
  
  # Read in NA Raster
  ndvi_raster_NA=stack(paste0(my_wd,"/NDVI_2013.gri"))
  ndvi_raster_composite=stack(paste0(my_wd,"/NDVI_2013_composite.gri"))
  
  # Read in NOAM
  NOAM=readOGR(dsn=my_wd,layer="NOAM")
  
}