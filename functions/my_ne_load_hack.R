# this is my hack of the original ne_load function to avoid readOGR console messages
# cengel - oct 2017
my_ne_load_hack <- function (scale = 110, type = "countries", category = c("cultural", 
                                                                      "physical", "raster"), destdir = tempdir(), file_name = NULL) 
{
  category <- match.arg(category)
  #returnclass <- match.arg(returnclass)
  if (is.null(file_name)) {
    file_name <- rnaturalearth::ne_file_name(scale = scale, type = type, 
                              category = category)
  }
  error_msg <- paste0("the file ", file_name, " seems not to exist in your local folder ", 
                      destdir, "\nDid you download it using ne_download()?")
  if (category == "raster") {
    file_tif <- file.path(destdir, file_name, paste0(file_name, 
                                                     ".tif"))
    if (!file.exists(file_tif)) 
      stop(error_msg)
    rst <- raster::raster(file_tif)
    return(rst)
  }
  else {
    if (!file.exists(file.path(destdir, paste0(file_name, 
                                               ".shp")))) 
      stop(error_msg)
    sp_object <- rgdal::readOGR(destdir, file_name, encoding = "UTF-8", 
                                stringsAsFactors = FALSE, use_iconv = TRUE, verbose = FALSE)
    sp_object@data[sp_object@data == "-99" | sp_object@data == 
                     "-099"] <- NA
    return(sp_object)
  }
}