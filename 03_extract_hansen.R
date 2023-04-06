
# Extract point-level deforestation data for each abundance study location from Hansen rasters (2000-2019)

################################################################################

pacman::p_load("raster", "rgdal", "sp")

### 1. Subset abundance data to primary and secondary vegetation sites
data <- read.csv("data/inla_input/abundance.csv") %>%
  subset(land_use == "primary vegetation" | land_use == "secondary vegetation") 

# subset to studies with unique coords for each site
studies <- NULL
for (i in unique(data$study_number)){
  
  df <- data %>% subset(study_number == i)
  df$coords <- paste(df$lat, df$lon, sep = " , ")
  
  if (length(unique(df$coords)) > 1){
    
    studies <- c(studies, i)
    
  }
  
}

# write to file
data %>% subset(study_number %in% studies) %>% write.csv("data/hansen/input/abundance_data.csv")

### 2. Extract points over deforestation rasters (provided as yearly rasters per grid)
# read in data and create sp dataframe
data <- read.csv("data/hansen/input/abundance_data.csv")
pts_df <- data[c(4,5)]
pts <- SpatialPoints(pts_df, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
proj4string(pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# list rasters to extract from
ras_grids <- c("00N_040W", "00N_050W", "00N_060W", "00N_070W", "00N_080W", "00N_090W",
               "10N_060W", "10N_070W", "10S_040W", "10S_050W", "20N_090W", "20S_050W", 
               "30S_060W")

for (j in ras_grids) {
  
  # ras files
  files <- list.files("data/hansen/ras_input", pattern = j, full.names = TRUE)
  
  df <- NULL
  
  for (i in files){
    
    ras <- raster::raster(i)
    
    year <- gsub(".tif", "", gsub(".*W", "", i))
    
    ras_name <- gsub("W.*", "", gsub(".tif", "", gsub("data/hansen/ras_input/", "", i)))
    
    # including buffer of 320m
    # raster values are either 1 (deforestation) or 0 (no deforestation)
    df_ras <- raster::extract(ras, pts, df = TRUE, buffer = 320, fun = sum) %>% 
      mutate(year = year) %>%
      dplyr::rename(point=ID)
    names(df_ras)[2] <- "cells_deforestation"
    
    df <- rbind(df, df_ras)
  }
  
  df$lon <- pts_df$lon
  df$lat <- pts_df$lat
  # write to file
  write.csv(df, paste0("data/hansen/buffer_output/320/", ras_name, ".csv"), row.names = FALSE)
  
}

### 3. Combine all raster data and recode years
files <- list.files("data/hansen/buffer_output/320", pattern = ".csv", full.names = TRUE)

deforest_data <- NULL

for (i in 1:length(files)){
  
  deforest_data <- read.csv(files[i]) %>% 
    # remove NA and zero loss i.e. year = 0
    dplyr::filter(!is.na(cells_deforestation)) %>% subset(year != 0) %>% # remove when updated
    dplyr::group_by(year, lon, lat) %>%
    # remove repeats (total number of cells deforested)
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

### 4. Temporally match site data to deforestation data
data <- read.csv("data/hansen/input/abundance_data.csv") %>%
  # subset to years of deforestation data (2000-2019)
  subset(sample_start_year > 2000)

sites <- unique(data$site_number)

data_all <- NULL 

# for each site calculate total deforestation in last 5 years (total number of cells)
for (i in 1:length(sites)){
  
  df <- data %>% subset(site_number == sites[i])
  study_year <- as.numeric(df$sample_start_year[1])
  
  # match by lon lat
  data_all <- deforest_data %>% subset(lon == df$lon[1] & lat == df$lat[1]) %>% 
    # subset to last 5 years
    subset(year == study_year | year == study_year-1 |
             year == study_year-2 | year == study_year-3|
             year == study_year-4) %>% 
    dplyr::group_by(lon, lat) %>%
    # total deforestation in last 5 years
    dplyr::summarise(cells_deforestation = sum(cells_deforestation)) %>%
    dplyr::ungroup() %>%
    mutate(study_number = df$study_number[1],
           site_number  = df$site_number[1]) %>%
    rbind(data_all)
  
}

### 5. Calculate proportion of deforestation within buffer (320m)
# raster resolution  = 25m x 25m
# cell area = 25x25 = 625m
cell_area <- 625
# buffer area (320 m) = pi * 320^2
buffer_area <- pi * 320^2

# deforestation = ncells deforested x cell area/total buffer area

data_all <- data_all %>% mutate(deforest_prop = (cells_deforestation * cell_area)/buffer_area) %>% 
  mutate(deforest_prop = deforest_prop * 100) 

### 6. Combine with rest of dataframe and write to file
data <- read.csv("data/hansen/input/abundance_data.csv") %>%
  subset(site_number %in% unique(data_all$site_number)) %>% 
  merge(data_all, by = c("site_number", "study_number", "lon", "lat"), all = TRUE) 

data %>% 
  write.csv("data/inla_input/deforestation_data.csv")
