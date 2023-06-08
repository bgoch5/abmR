#'
#' Downloads example NDVI data
#'
#' Warning: this function will download to your hard drive
#' (to a location specified by your current working directory) the below
#' files, totaling approximately 620 MB. Please do not attempt to use if you
#' have insufficient hard drive space
#' or Random Access Memory (RAM). Objects to be downloaded are listed under
#' "details". The first time you use this function, you will be directed to your
#' browser and required to sign in to your Google account to connect to the 
#' Tidyverse API. If you use the function a second time, you may simply follow
#' the prompts and enter a number corresponding to the previous accounts listed.
#'
#' @import googledrive purrr raster
#' @details
#'  \itemize{
#'  \item{NDVI_2013_NA:}{  A raster stack containing daily NDVI data for North America,
#' on a .05 x .05 degree grid. Data runs from 8/26/2013-9/21/2013.}
#'  \item{NDVI_2013_NA_composite:}{  Single layer raster formed by taking mean of NDVI_2013_NA}
#'  \item{NDVI_2013_Europe:}{  A raster stack containing daily NDVI data for Europe,
#' on a .05 x .05 degree grid. Data runs from 8/26/2013-9/21/2013.}
#'  \item{NDVI_2013_Europe_composite:}{  Single layer raster formed by taking mean of NDVI_2013_Europe}
#'  }
#' @source
#' Vermote, Eric; NOAA CDR Program. (2019): NOAA Climate Data Record (CDR)
#' of AVHRR Normalized Difference Vegetation Index (NDVI), Version 5.
#' NOAA National Centers for Environmental Information.
#' https://doi.org/10.7289/V5ZG6QH9. Accessed 12/26/2020.
#' @return No return value, called to download files to hard drive.
#' @export

get_ex_data <- function() {
  folder_url <- "https://drive.google.com/drive/folders/1hV1kmtw4u8fh1crAsxI7zBaySqgjTZaJ"
  folder <- drive_get(as_id(folder_url))
  my_files <- drive_ls(folder)
  my_wd <- getwd() # Record working directory for later use

  # Download files, will save to hard drive
  walk(my_files$id, ~ drive_download(as_id(.x), overwrite = TRUE))
}
