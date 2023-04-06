
# Descriptive statistics in main text

################################################################################

pacman::p_load("dplyr", "tibble", "sp", "stringr")

# Load functions
load("functions/extract_estimates2.RData")
load("functions/extract_deforest_estimates.RData")
################################################################################

## Dataset
data <- read.csv("data/inla_input/abundance.csv") 

# Total number of records 
data %>% nrow()

# Total number of sites 
length(unique(data$site_number))

# Total number of studies
length(unique(data$study_number)) 

# Number of sites by land use
# combine cropland, pasture and plantation into single land use
data %>% mutate(land_use = ifelse(land_use == "pasture" | 
                                    land_use == "plantation" |
                                    land_use == "cropland", "managed", land_use)) %>%
  dplyr::group_by(land_use) %>%
  dplyr::summarise(n = length(unique(site_number)))
# 292/632

# Number of records by land use
# combine cropland, pasture and plantation into single land use
data %>% mutate(land_use = ifelse(land_use == "pasture" | 
                                    land_use == "plantation" |
                                    land_use == "cropland", "managed", land_use)) %>%
  subset(land_use == "primary vegetation") %>% 
  nrow()
# 3835/10244 = 37%

# Number of countries
length(unique(data$country))

# Number of records in Brazil
data %>% subset(country == "Brazil") %>% nrow()
round(data %>% subset(country == "Brazil") %>% nrow()/data %>% nrow() *100)

# Count no. Amazonian vs. non-Amazonian sites
amazon <- rgdal::readOGR("data/amazon_poly", "amapoly_ivb")

pts <- data %>% dplyr::select(lon, lat) %>% 
  SpatialPoints(proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Extract over amazon
amazon_pts <- over(pts, amazon)
amazon_pts$amazon <- 1
# label points outside polygon
amazon_pts$amazon[is.na(amazon_pts$AREA)] <- 0

# Combine with data
data <- data %>% 
  mutate(amazon = amazon_pts$amazon)
write.csv(amazon_pts, "data/input/americas/inla_input/amazon_pts.csv")

data %>%
  dplyr::group_by(amazon, site_number) %>%
  dplyr::summarise(amazon = unique(amazon)) %>%
  dplyr::count(amazon) 
#431/632

# Count no. Atlantic forest sites
atlantic_forest <- rgdal::readOGR("data/atlantic_forest_shp/limites_wgs94", "limite_ma_wwf_200ecprregions_wgs84")

pts <- data %>% dplyr::select(lon, lat) %>% 
  SpatialPoints(proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Extract over atlantic forest
atlantic_forest_pts <- over(pts, atlantic_forest)
atlantic_forest_pts$atlantic_forest <- 1
# label points outside polygon
atlantic_forest_pts$atlantic_forest[is.na(atlantic_forest_pts$G200_BIOME)] <- 0

# Combine with data
data <- data %>% 
  mutate(atlantic_forest = atlantic_forest_pts$atlantic_forest)

data %>%
  dplyr::group_by(atlantic_forest, site_number) %>%
  dplyr::summarise(atlantic_forest = unique(atlantic_forest)) %>%
  dplyr::count(atlantic_forest) 
#111/632

# Number of biomes
unique(data$biome)

# Number of forested sites
data %>% subset(biome == unique(data$biome)[1] |
                  biome == unique(data$biome)[2] |
                  biome == unique(data$biome)[7]) %>% 
  dplyr::summarise(n = length(unique(site_number)))
# 609/632

# Total number of species
length(unique(data$species))

# Number of Aedes species
data %>% subset(genus == "Aedes") %>% 
  dplyr::summarise(n = length(unique(species)))
# 33/91 

# Number of Anopheles species
data %>% subset(genus == "Anopheles") %>% 
  dplyr::summarise(n = length(unique(species)))
# 58/91

################################################################################

## Effect of land use change on species richness and abundance
##### difference in Anopheles and Aedes richness in urban areas
files <- list.files("models/richness", pattern = "_lui", full.names = TRUE)
files <- c(files[1], files[2])

data <- NULL
for (i in files){
  
  load(i)
  
  if (grepl("anopheles", i)==TRUE){
    genus_name <- "Anopheles"
  } else {
    genus_name <- "Aedes"
  }
  
  df <- extract_estimates2(mod, "intensity") %>% mutate(genus = genus_name)
  data <- rbind(data, df)
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

data %>% subset(land_use == "urban") %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

##### difference in Anopheles abundance in urban and managed
files <- list.files("models/abundance", pattern = "_lui", full.names = TRUE)
files <- c(files[1], files[2])

data <- NULL
for (i in files){
  
  load(i)
  
  if (grepl("anopheles", i)==TRUE){
    genus_name <- "Anopheles"
  } else {
    genus_name <- "Aedes"
  }
  
  df <- extract_estimates2(mod, "intensity") %>% mutate(genus = genus_name)
  data <- rbind(data, df)
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

data %>% subset(genus == "Anopheles") %>% 
  subset(land_use == "managed" | land_use == "urban") %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

# Total species richness
files <- list.files("models/richness", pattern = "_lui", full.names = TRUE)
load(files[3])

data <- extract_estimates2(mod, "intensity")

# Set primary minimal category as baseline
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 

# Difference from intercept
data[,2:4] = (exp(data[2:4]) - 1)*100

data %>% mutate(mean = round(mean, 0),
                lci = round(lci, 1),
                uci = round(uci, 1)) %>%
  subset(land_use == "urban")

################################################################################

## Species-specific abundance responses

files <- list.files("models/abundance/species", pattern = "_lui_mod", full.names = TRUE)

data <- NULL
for (i in files){
  
  load(i)
  
  # Create df of estimates
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% 
                  mutate(species = gsub("_", " ", gsub("_lui_mod.RData", "" , gsub("models/abundance/species/", "", i)))))
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

# Aedes
data %>% subset(species == "aedes aegypti") %>%
  dplyr::select(land_use, species, mean, lci, uci) %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

data %>% subset(species == "aedes albopictus") %>%
  dplyr::select(land_use, species, mean, lci, uci) %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

data %>% subset(species == "aedes aegypti" | species == "aedes albopictus") %>%
  subset(land_use == "secondary vegetation" | land_use == "primary vegetation-substantial") %>% 
  dplyr::select(land_use, species, mean, lci, uci) %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

data %>% subset(species == "aedes aegypti" | species == "aedes albopictus") %>%
  subset(land_use == "urban") %>% dplyr::select(land_use, species, mean, lci, uci) %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

data %>% subset(species == "aedes scapularis" | species == "aedes serratus") %>%
  subset(land_use == "managed" | land_use == "primary vegetation-substantial") %>% 
  dplyr::select(land_use, species, mean, lci, uci) %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

#### Anopheles
data %>% subset(species == "anopheles albitarsis") %>%
  subset(land_use == "managed") %>% dplyr::select(land_use, species, mean, lci, uci) %>%
  mutate(mean = round(mean, 0),
         lci = round(lci, 1),
         uci = round(uci, 1))

################################################################################

# Effect of deforestation
# Richness models
files <- list.files("models/richness", pattern = "_l_deforestation_mod", full.names = TRUE)

data_rich <- NULL
for (i in files){
  
  load(i)
  
  if (grepl("anopheles", i)==TRUE){
    genus_name <- "Anopheles"
  } else {if (grepl("aedes", i)==TRUE){
    genus_name <- "Aedes"
  } else{
    genus_name <- "Total"
  }
  }
  df <- extract_deforest_estimates(mod, "linear") %>% mutate(genus = genus_name,
                                                             model = "richness")
  data_rich <- rbind(data_rich, df)
}

files <- list.files("models/abundance", pattern = "_l_deforestation_mod", full.names = TRUE)

data_abun <- NULL
for (i in files){
  
  load(i)
  
  if (grepl("anopheles", i)==TRUE){
    genus_name <- "Anopheles"
  } else {if (grepl("aedes", i)==TRUE){
    genus_name <- "Aedes"
  } else{
    genus_name <- "Total"
  }
  }
  df <- extract_deforest_estimates(mod, "linear") %>% mutate(genus = genus_name,
                                                             model = "abundance")
  data_abun <- rbind(data_abun, df)
}

df <-data_rich %>% rbind(data_abun) %>% mutate(mean = round(mean, 2),
                                          lci  = round(lci, 2),
                                          uci  = round(uci, 2)) 
df

# increase in richness with every % increase in deforestation
round((1-exp(subset(df, df$model == "richness" & genus == "Anopheles")$mean))*100)

## Species abundance models
files <- list.files("models/abundance/species", pattern = "_l_deforestation_mod", full.names = TRUE)

data <- NULL
for (i in files){
  
  load(i)
  
  df <- extract_deforest_estimates(mod, "linear") %>% mutate(species = gsub("_l.*", "", gsub(".*species/", "", i))) %>%
    mutate(species = as.factor(str_to_sentence(gsub("_", " ", species))))
  data <- rbind(data, df)
}

df <- data %>% mutate(mean = round(mean, 2),
                lci  = round(lci, 2),
                uci  = round(uci, 2)) 

df %>% dplyr::arrange(-mean)

# increase in abundance with unit increase in deforestation
round((exp(subset(df, df$species == "Anopheles darlingi")$mean)-1)*100)
round((exp(subset(df, df$species == "Anopheles albitarsis")$mean)-1)*100)
round((exp(subset(df, df$species == "Aedes serratus")$mean)-1)*100)

################################################################
# Discussion
load("functions/aggregate_lui.RData")

# number of an. darlingi secondary vegetation sites
data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") %>% mutate(species_name = tolower(species_name)) %>% 
  subset(species_name == "anopheles darlingi" & LUI == "secondary vegetation")
length(unique(data$site_number))

# number of secondary vegetation sites
secondary_sites <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") %>% 
  subset(LUI == "secondary vegetation")
all_sites <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") 
length(unique(secondary_sites$site_number))/length(unique(all_sites$site_number)) *100
