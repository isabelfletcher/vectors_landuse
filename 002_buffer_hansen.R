library(raster)
library(sp)
library(dplyr)
library(rgdal)

# This script extracts deforestation data, based on the abundance data points, from Hansen rasters of annual forest loss

############################################################
# Extract rasters for each year of data

data <- read.csv("data/hansen/input/abundance_data.csv")

pts <- data[c(4,5)]

pts <- SpatialPoints(pts, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
proj4string(pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

files <- list.files("data/hansen/input/all", pattern = ".tif", full.names = TRUE)

j <- 1

for (i in 1:19){
  
  ras <- raster::raster(files[j])
  
  ras_name <- gsub(".tif", "", gsub("data/hansen/input/all/Hansen_GFC-2019-v1.7_lossyear_", "", files[j]))
  
  ras_sub <- ras == i
  writeRaster(ras_sub, filename = paste0("data/hansen/ras_input/", ras_name,  i, ".tif"), overwrite=T)
  
  rm(ras_sub)
  rm(ras)
}

########################################################################################################################
# Extract data 

# Mosquito flight ranges https://www.sciencedirect.com/science/article/pii/S0075951113001011 
# take average of aedes and anopheles sp mentioned in table 4: 89+541.9/2 = 315.5m = 320m

data <- read.csv("data/hansen/input/abundance_data.csv")
pts_df <- data[c(4,5)]
pts <- SpatialPoints(pts_df, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
proj4string(pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# rasters to extract from
ras_grids <- c("00N_040W", "00N_050W", "00N_060W", "00N_070W", "00N_080W", "00N_090W",
               "10N_060W", "10N_070W", "10S_040W", "10S_050W", "20N_090W", "20S_050W", 
               "30S_060W")

for (j in ras_grids) {
  
  files <- list.files("data/hansen/ras_input", pattern = j, full.names = TRUE)
  
  df <- NULL

  for (i in files){
  
  ras <- raster::raster(i)
  
  year <- gsub(".tif", "", gsub(".*W", "", i))
  
  ras_name <- gsub("W.*", "", gsub(".tif", "", gsub("data/hansen/ras_input/", "", i)))
  
  df_ras <- raster::extract(ras, pts, df = TRUE, buffer = 320, fun = sum) %>% 
    mutate(year = year) %>%
    dplyr::rename(point=ID)
  names(df_ras)[2] <- "cells_deforestation"
  
  df <- rbind(df, df_ras)
  }

  df$lon <- pts_df$lon
  df$lat <- pts_df$lat
  write.csv(df, paste0("data/hansen/buffer_output/320/", ras_name, ".csv"))
  
}

########################################################################################################################
# combine data
files <- list.files("data/hansen/buffer_output/320", pattern = ".csv", full.names = TRUE)

deforest_data <- NULL

for (i in 1:length(files)){
  
  deforest_data <- read.csv(files[i]) %>% dplyr::select(-X) %>%
  dplyr::filter(!is.na(cells_deforestation)) %>% subset(year != 0) %>% # remove when updated
  dplyr::group_by(year, lon, lat) %>%
  # remove repeats
  dplyr::summarise(cells_deforestation = max(cells_deforestation)) %>% 
  dplyr::ungroup() %>%
  rbind(deforest_data) %>% as.data.frame() %>%
    mutate(year = as.factor(year)) %>% 
    mutate(year = recode(year,
                         "1" = "2001",
                         "2" = "2002",
                         "3" = "2003",
                         "4" = "2004",
                         "5" = "2005",
                         "6" = "2006",
                         "7" = "2007",
                         "8" = "2008",
                         "9" = "2009",
                         "10" = "2010",
                         "11" = "2011",
                         "12" = "2012",
                         "13" = "2013",
                         "14" = "2014",
                         "15" = "2015",
                         "16" = "2016",
                         "17" = "2017",
                         "18" = "2018",
                         "19" = "2019"))
  
}

########################################################################################################################
# Combine with rest of dataframe and temporally match data
data <- read.csv("data/hansen/input/abundance_data.csv") %>%
# subset to years of deforestation data (2000-2019)
  subset(sample_start_year > 2000)

sites <- unique(data$site_number)

data_all <- NULL 

for (i in 1:length(sites)){
  
  df <- data %>% subset(site_number == sites[i])
  study_year <- as.numeric(df$sample_start_year[1])
  
  data_all <- deforest_data %>% subset(lon == df$lon[1] & lat == df$lat[1]) %>% 
  # subset to last 5 years
  subset(year == study_year | year == study_year-1 |
           year == study_year-2 | year == study_year-3|
           year == study_year-4) %>% 
  dplyr::group_by(lon, lat) %>%
  dplyr::summarise(cells_deforestation = sum(cells_deforestation)) %>%
  dplyr::ungroup() %>%
  mutate(study_number = df$study_number[1],
         site_number  = df$site_number[1]) %>%
    rbind(data_all)

}

########################################################################################################################
# calculate proportion of deforestation
# raster resolution  = 25m x 25m
# cell area = 25x25 = 625m
cell_area <- 625
# buffer area (320 m) = pi * 320^2
buffer_area <- pi * 320^2

# deforestation = ncells deforested x cell area/total buffer area

data_all <- data_all %>% mutate(deforest_prop = (cells_deforestation * cell_area)/buffer_area) %>% 
  mutate(deforest_prop = deforest_prop * 100) 


########################################################################################################################
# combine with rest of dataframe 
data <- read.csv("data/hansen/input/abundance_data.csv") %>%
  subset(site_number %in% unique(data_all$site_number)) %>% 
  merge(data_all, by = c("site_number", "study_number", "lon", "lat"), all = TRUE) 

data %>% 
  write.csv("data/inla_input/deforestation_data.csv")
