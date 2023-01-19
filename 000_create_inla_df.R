pacman::p_load("dplyr")

# Read in data
# subset to unassigned land use categories and where sampling effort not available
data <- read.csv("data/data.csv")

data <- data %>% subset(land_use != "cannot decide",
                        land_use_intensity != "cannot decide") %>%
  subset(sampling_effort != "not specified") %>%
  subset(sampling_effort != "NA") %>%
  dplyr::select(-lat_lon_source, -location_description, -climate_description, -water_body, -land_cover,
                -land_use_source, -intervention, -study_date, -sub_sites_sampled,
                -container, -container_size, -breeding_habitat, -baited, -bait_type, -identification_method, -coords) %>% 
  dplyr::mutate(sampling_effort = as.numeric(sampling_effort))

# For each study and observation assign a study sample, where sampling methodology is consistent
load("functions/assign_study_sample.RData")

df <- NULL
for (i in unique(data$study_number)){
  
  df <- rbind(df,
              data %>% subset(study_number == i) %>%
                assign_study_sample())
}

data <- df

# Scale sampling effort
load("functions/scale_sampling_effort.RData")

#data <- scale_sampling_effort(data) %>%
  # adjust measurement by sampling effort
 # mutate(measurement_adj = as.numeric(measurement)/as.numeric(rescaled_sampling_effort))

# test not rescaling sampling effort
data <- data %>%
  # adjust measurement by sampling effort
  mutate(measurement_adj = as.numeric(measurement)/as.numeric(sampling_effort))

# check 
which(is.na(data$measurement_adj))
na_df <- data[is.na(data$measurement_adj),]
unique(na_df$measurement)
unique(na_df$rescaled_sampling_effort)

# set NAs to zero
data$measurement_adj[is.na(data$measurement_adj)] <- 0

# Subset to Anopheles and Aedes, combine land use categories

# Non-aggregated land use categories
data_all <- data %>%
  # subset to Anopheles and Aedes spp.
  subset(genus == "Anopheles" | genus == "Aedes") %>% 
  mutate(species_name = paste0(genus, " ", species)) %>%
  # create single categorical variable of land use and intensity
  mutate(LUI = paste0(land_use, "-", land_use_intensity)) 

# check for studies with only one LUI category
for (i in unique(data_all$study_number)){
  
  df <- data_all %>% subset(study_number == i)
  
  if(length(unique(df$LUI))==1){
    print(i)
  }
  
}


# remove where study only has one land use type
data_lu <- NULL 
for (i in unique(data_all$study_number)){
  
  df <- data_all %>% subset(study_number == i)
  
  if(length(unique(df$LUI))>1){
    data_lu <- rbind(data_lu, df)
  }
  
}

# Write to file
data_lu %>% write.csv("data/data_comb.csv")

# Combined land use categories
data <- data %>%
# subset to Anopheles and Aedes spp.
  subset(genus == "Anopheles" | genus == "Aedes") %>% 
  mutate(species_name = paste0(genus, " ", species)) %>%
  # combine cropland, pasture and plantation into single land use
  mutate(land_use = ifelse(land_use == "pasture" | 
                             land_use == "agriculture" | 
                             land_use == "plantation" |
                             land_use == "cropland", "managed", land_use)) %>%
  # combine light and intense land use into single category 
  mutate(land_use_intensity_agg = ifelse(land_use_intensity == "light" | 
                                     land_use_intensity == "intense", 
                                     "substantial", land_use_intensity)) %>%
    # create single categorical variable of land use and intensity
  mutate(LUI = paste0(land_use, "-", land_use_intensity_agg)) 

# check for studies with only one LUI category
for (i in unique(data$study_number)){
  
  df <- data %>% subset(study_number == i)
  
  if(length(unique(df$LUI))==1){
    print(i)
  }
  
}


# remove where study only has one land use type
data_lu <- NULL 
for (i in unique(data$study_number)){
  
  df <- data %>% subset(study_number == i)
  
  if(length(unique(df$LUI))>1){
    data_lu <- rbind(data_lu, df)
  }
  
}
data <- data_lu

# Write to file
data %>% write.csv("data/inla_input/abundance.csv")

####################################################################
# Species richness
data <- read.csv("data/inla_input/abundance.csv")

# Calculate site-level species richness, excluding zero records
richness <- data %>% subset(species_studied == "multiple") %>% 
  subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample)

richness_data <- NULL
for (i in unique(richness$site_number)){
  richness_data <- rbind(richness_data,
     richness %>% subset(site_number == i) %>%
      dplyr::group_by(site_number, study_block, study_number, study_sample) %>%
       dplyr::summarise(richness = length(unique(species_name))))
}

# Combine with all data, excluding zero records
df <- 
  data %>% subset(species_studied == "multiple") %>% 
  subset(measurement_adj > 0) %>% 
  dplyr::group_by(study_number, site_number, study_block, study_sample) %>%
  dplyr::summarise(reference = unique(reference),
                   lat = unique(lat),
                   lon = unique(lon),
                   country = unique(country),
                   biome     = unique(biome),
                   land_use  = unique(land_use),
                   land_use_specific  = unique(land_use_specific),
                   land_use_intensity  = unique(land_use_intensity),
                   LUI  = unique(LUI)) 

df_all <-
  merge(df, richness_data, by = c("site_number", "study_block", "study_sample", "study_number"))
unique(df_all$richness)

## Add Amazon points
amazon <- rgdal::readOGR("data/amazon_poly", "amapoly_ivb")

pts <- df_all %>% dplyr::select(lon, lat) %>% 
  SpatialPoints(proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Extract over amazon
amazon_pts <- over(pts, amazon)
amazon_pts$amazon <- 1
# label points outside polygon
amazon_pts$amazon[is.na(amazon_pts$AREA)] <- 0

df_all %>% 
  mutate(amazon = amazon_pts$amazon) %>% 
  # Write to file
  write.csv("data/inla_input/richness.csv")

df_all %>% 
  mutate(amazon = amazon_pts$amazon) %>% 
  write.csv("data/inla_input/richness.csv")

##### Anopheline species richness 
# Calculate site-level species richness, excluding zero records
richness <- data %>% subset(species_studied == "multiple") %>% 
  subset(genus == "Anopheles") %>%
  subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample)

richness_data <- NULL
for (i in unique(richness$site_number)){
  richness_data <- rbind(richness_data,
                         richness %>% subset(site_number == i) %>%
                           dplyr::group_by(site_number, study_block, study_number, study_sample) %>%
                           dplyr::summarise(richness = length(unique(species_name))))
}

# Combine with all data, excluding zero records
df <- 
  data %>% subset(species_studied == "multiple") %>% 
  subset(measurement_adj > 0) %>% subset(genus == "Anopheles") %>%
  dplyr::group_by(study_number, site_number, study_block, study_sample) %>%
  dplyr::summarise(reference = unique(reference),
                   lat = unique(lat),
                   lon = unique(lon),
                   country = unique(country),
                   biome     = unique(biome),
                   land_use  = unique(land_use),
                   land_use_specific  = unique(land_use_specific),
                   land_use_intensity  = unique(land_use_intensity),
                   LUI  = unique(LUI)) 

df_all <-
  merge(df, richness_data, by = c("site_number", "study_block", "study_sample", "study_number"))
unique(df_all$richness)

## Add Amazon points
amazon <- rgdal::readOGR("data/amazon_poly", "amapoly_ivb")

pts <- df_all %>% dplyr::select(lon, lat) %>% 
  SpatialPoints(proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Extract over amazon
amazon_pts <- over(pts, amazon)
amazon_pts$amazon <- 1
# label points outside polygon
amazon_pts$amazon[is.na(amazon_pts$AREA)] <- 0

df_all %>% 
  mutate(amazon = amazon_pts$amazon) %>% 
  # Write to file
  write.csv("data/inla_input/ano_richness.csv")

df_all %>% 
  mutate(amazon = amazon_pts$amazon) %>% 
  write.csv("data/inla_input/ano_richness.csv")

##### Aedes species richness 
# Calculate site-level species richness, excluding zero records
richness <- data %>% subset(species_studied == "multiple") %>% 
  subset(genus == "Aedes") %>%
  subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample)

richness_data <- NULL
for (i in unique(richness$site_number)){
  richness_data <- rbind(richness_data,
                         richness %>% subset(site_number == i) %>%
                           dplyr::group_by(site_number, study_block, study_number, study_sample) %>%
                           dplyr::summarise(richness = length(unique(species_name))))
}


# Combine with all data, excluding zero records
df <- 
  data %>% subset(species_studied == "multiple") %>% 
  subset(measurement_adj > 0) %>% subset(genus == "Aedes") %>%
  dplyr::group_by(study_number, site_number, study_block, study_sample) %>%
  dplyr::summarise(reference = unique(reference),
                   lat = unique(lat),
                   lon = unique(lon),
                   country = unique(country),
                   biome     = unique(biome),
                   land_use  = unique(land_use),
                   land_use_specific  = unique(land_use_specific),
                   land_use_intensity  = unique(land_use_intensity),
                   LUI  = unique(LUI)) 

df_all <-
  merge(df, richness_data, by = c("site_number", "study_block", "study_sample", "study_number"))
unique(df_all$richness)

## Add Amazon points
amazon <- rgdal::readOGR("data/amazon_poly", "amapoly_ivb")

pts <- df_all %>% dplyr::select(lon, lat) %>% 
  SpatialPoints(proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Extract over amazon
amazon_pts <- over(pts, amazon)
amazon_pts$amazon <- 1
# label points outside polygon
amazon_pts$amazon[is.na(amazon_pts$AREA)] <- 0

df_all %>% 
  mutate(amazon = amazon_pts$amazon) %>% 
# Write to file
  write.csv("data/inla_input/aed_richness.csv")
