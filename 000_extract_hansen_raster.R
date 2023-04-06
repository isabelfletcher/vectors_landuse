library(raster)
library(sp)
library(dplyr)
library(rgdal)

############################################################
# Extract rasters for each year of data (2001-2019)
files <- list.files("data/hansen/input/all", pattern = ".tif", full.names = TRUE)

# change from 1:length(files), for each grid/raster
j <- 1

for (i in 1:19){
  
  # for each timepoint (2001-2019) create a yearly raster
  ras <- raster::raster(files[j])
  
  ras_name <- gsub(".tif", "", gsub("data/hansen/input/all/Hansen_GFC-2019-v1.7_lossyear_", "", files[j]))
  
  ras_sub <- ras == i
  writeRaster(ras_sub, filename = paste0("data/hansen/ras_input/", ras_name,  i, ".tif"), overwrite=T)
  
  rm(ras_sub)
  rm(ras)
}
