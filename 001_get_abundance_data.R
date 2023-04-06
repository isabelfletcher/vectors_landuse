pacman::p_load("dplyr", "rnaturalearth", "sf", "ggplot2", "sp")

# Subset abundance data to primary and secondary vegetation sites
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
  
data %>% subset(study_number %in% studies) %>% write.csv("data/hansen/input/abundance_data.csv")
  
# Visualise points on map
world <- ne_countries(scale = "medium", returnclass = "sf")

data <- read.csv("data/hansen/input/abundance_data.csv") %>% 
  dplyr::select(lon, lat)

sp_data <- read.csv("data/hansen/input/abundance_data.csv") %>% 
  dplyr::select(lon, lat) %>%
  SpatialPoints(proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>%
  st_as_sf()
  
ggplot(data = world) +
    geom_sf() +
    geom_sf(data = sp_data, colour = "red",
            size = 1.2, alpha = 1) +
    coord_sf(xlim = c(-120, -30), ylim = c(-57, 35)) +
    theme_light() 
