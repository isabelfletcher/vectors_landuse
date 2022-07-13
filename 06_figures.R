
# Plot results of analysis

################################################################################

pacman::p_load("dplyr", "gplot2", "tidyr", "rgdal", "sf", 
               "rnaturalearth", "rgeos", "scales",
               "raster", "stringr", "ggpubr",
               "cowplot", "gtable", "gridExtra",
               "grid", "scales", "tibble", "forcats",
               "tidyverse", "ggtext")

load("functions/extract_estimates.RData")
load("functions/extract_estimates2.RData")
load("functions/plot_estimates.RData")
load("functions/plot_estimates_cross.RData")
load("functions/extract_deforest_estimates.RData")
load("functions/aggregate_lui.RData")

################################################################################
############# Figure 1

############### Map of sites 
world <- ne_countries(scale = "medium", returnclass = "sf")

data <- read.csv("data/inla_input/abundance.csv") %>% 
  dplyr::select(lat, lon, land_use, land_use_intensity_agg) %>%
  unique() %>%
  mutate(land_use = factor(land_use, levels = c("primary vegetation", "secondary vegetation", "managed", "urban")),
         land_use_intensity_agg = factor(land_use_intensity_agg, 
                                     levels = c("minimal", "substantial")))

amazon <- rgdal::readOGR("data/amazon_poly", "amapoly_ivb")
atlantic_forest <- rgdal::readOGR("data/atlantic_forest_shp/limites_wgs94", "limite_ma_wwf_200ecprregions_wgs84")

p <- ggplot(data=world) +
  geom_sf(colour="white") +
  geom_sf(data = st_as_sf(amazon), colour = "transparent", fill = "#0a9396", alpha = 0.3) +
  geom_sf(data = st_as_sf(atlantic_forest), colour = "transparent", fill = "#0a9396", alpha = 0.3) +
  coord_sf(xlim = c(-120, -30), ylim = c(-57, 35), expand = FALSE) +
  geom_point(aes(x = lon, y = lat, fill = land_use), data = data, shape = 21, colour = "white", alpha = 0.7) +
  theme_light() +
  labs(subtitle = "A") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.23,0.6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.subtitle = element_text(face="bold", size = 13)) +
  scale_fill_manual("Land use", values = c("#99D939", "#5BA39D", "#E69F00","#4A4F97"))

# Number of studies and by land use
n_studies <- read.csv("data/inla_input/abundance.csv") %>%
  dplyr::summarise(n = length(unique(study_number)))

n_sites <- read.csv("data/inla_input/abundance.csv")   %>%
  subset(genus == "Anopheles" | genus == "Aedes") %>%
  subset(land_use != "cannot decide",
         land_use_intensity != "cannot decide") %>%
  dplyr::summarise(n = length(unique(site_number)))

n_lu <- read.csv("data/inla_input/abundance.csv")  %>%
  dplyr::group_by(land_use) %>%
  dplyr::summarise(n = length(unique(site_number)))
n_lu <- n_lu[c(2,3,1,4),]


# list of text justifications
x_just <- c(0.035, 0.035, 0.125, 0.125, 0.125, 0.125)

table <- data.frame(metric = c(rep("Number of studies", length(n_studies$n)),
                               rep("Number of sites", length(n_sites$n)),
                               c(n_lu$land_use)),
                    n = c(n_studies$n, n_sites$n, n_lu$n)) %>%
  tableGrob(theme = ttheme_minimal(base_size = 8,
                                   padding = unit(c(4, 4), "mm"),
                                   core=list(fg_params=list(hjust=0, x=as.vector(x_just)))), rows = NULL, cols = NULL)
g_tab <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 1.5)),
                         t = 1, b = nrow(table), l = 1, r = ncol(table))
g_tab <- gtable_add_grob(g_tab,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 1, lty=3)),
                         t = 2, b = nrow(table), l = 1, r = ncol(g_tab))
g_tab <- gtable_add_grob(g_tab,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 1, lty=3)),
                         t = 3, b = nrow(table), l = 1, r = ncol(g_tab))

map <- ggdraw() +
  draw_plot(p) +
  draw_plot(g_tab, 0.17, 0.15, width = 0.2, height = 0.2) 

############### Data distribution
# Site distribution by biome
data <- read.csv("data/inla_input/abundance.csv") 

# Add Amazonian vs. non-Amazonian sites
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

p1 <- 
data %>%
  drop_na(biome) %>%
  # relabel biomes
  mutate(biome = dplyr::recode(biome, "Tropical & Subtropical Grasslands, Savannas & Shrublands" = "Grassland & shrubland",
                               "Tropical & Subtropical Moist Broadleaf Forests" = "Forest",
                               "Tropical & Subtropical Coniferous Forests"= "Forest",
                               "Tropical & Subtropical Dry Broadleaf Forests"= "Forest",
                               "Tropical & Subtropical Grasslands, Savannas & Shrublands" = "Grassland & shrubland",
                               "Tropical & Subtropical Moist Broadleaf Forests"= "Forest",
                               "Temperate Grasslands, Savannas & Shrublands" = "Grassland & shrubland")) %>%
  dplyr::group_by(amazon, site_number) %>%
  dplyr::summarise(biome = unique(biome)) %>%
  dplyr::group_by(amazon) %>%
  dplyr::count(biome) %>%
  dplyr::mutate(amazon = dplyr::recode(amazon, "0" = "extra-Amazonian",
                                       "1" = "Amazonian")) %>% 
  ggplot() +
  geom_bar(aes(x = amazon, fill = biome, y = n), stat = "identity") +
  xlab("") +
  ylab("Number of sites") +
  labs(subtitle = "B") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 500)) +
  scale_fill_manual("Ecoregion", values = c("#304e52", "#97b669", "#4099af")) +
  theme_light() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.70),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.subtitle = element_text(face="bold", size = 13)) 

################
# Proportion of sp. richness
total_richness <- read.csv("data/inla_input/abundance.csv") %>%
  # exclude zero records
  subset(measurement > 0) %>%
  mutate(species_name = paste(genus, species, sep = " ")) %>% 
  dplyr::summarise(richness = length(unique(species_name)))

# richness per genus
cols <- brewer.pal(11, "PiYG") 
cols <- c("#5d1b3f", "#4099af")

p2 <- read.csv("data/inla_input/abundance.csv") %>%
  subset(genus == "Anopheles" | genus == "Aedes") %>%
  # exclude zero records
  subset(measurement > 0) %>%
  mutate(species_name = paste(genus, species, sep = " ")) %>%
  dplyr::group_by(genus) %>% 
  dplyr::summarise(richness = length(unique(species_name))) %>%
  mutate(prop_richness = richness/total_richness$richness *100,
         prop_richness = richness/total_richness$richness *100) %>%
  ggplot() +
  geom_bar(aes(x = genus, fill = genus, y = prop_richness), stat = "identity", alpha = 0.6) +
  xlab("") +
  ylab("Proportion of total\nspecies richness (%)") +
  labs(subtitle = "C") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_fill_manual(" ", values = c(cols[1], cols[2])) +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = "italic"),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.subtitle = element_text(face="bold", size = 13))

# Plot 
p <- grid.arrange(p1, p2, nrow = 2)
p <- grid.arrange(map, p, ncol = 2)
ggsave("figures/main/figure1", p, device = "tiff", width = 8, height = 4.5, units = "in", dpi = 500)

################################################################################
############# Figure 2
# Anopheles and Aedes abundance
files <- list.files("models/abundance", pattern = "_lui_mod", full.names = TRUE)
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

cols <- c("#5d1b3f", "#4099af")

# Plot
abundance_plot <- plot_estimates(data, "intensity", "abundance") +
  scale_fill_manual("",
                    values = c(cols[1], cols[2])) +
  scale_colour_manual("",
                      values = c(cols[1], cols[2])) +
  theme(legend.position = "bottom") +
  #guides(shape = "none") +
  ylim(c(-60, 40))

# Anopheles and Aedes species richness
files <- list.files("models/richness", pattern = "_lui_mod", full.names = TRUE)
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

cols <- c("#5d1b3f", "#4099af")

# Plot
richness_plot <- plot_estimates(data, "intensity", "richness") +
  scale_fill_manual("",
                    values = c(cols[1], cols[2])) +
  scale_colour_manual("",
                      values = c(cols[1], cols[2])) +
  theme(legend.position = "bottom") +
  #guides(shape = "none") +
  ylim(c(-60, 40))

ggpubr::ggarrange(richness_plot,
                  abundance_plot,
                  labels = c("A", "B"),
                  common.legend = TRUE,
                  legend = "bottom")
ggsave("figures/main/figure2", device = "tiff", width = 8.5, height = 4, units = "in", dpi = 500) 

################################################################################
############# Figure 3
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

cols <- c("#479A97", "#a37800", "#ffbe0c", "#6741ad")

# Plot
ggpubr::ggarrange(
  data %>% mutate(g = gsub(" .*", "", species)) %>% 
    subset(g == "aedes") %>% 
    mutate(species = str_to_sentence(species)) %>%
    plot_estimates("intensity", "abundance") +
    facet_wrap(~species, scales = "fixed") +
               #scales = "free_y") +
    scale_fill_manual("",
                      values = cols) +
    scale_colour_manual("",
                        values = cols) +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "italic")),
  
  data %>% mutate(g = gsub(" .*", "", species)) %>% 
    subset(g == "anopheles") %>% 
    mutate(species = str_to_sentence(species)) %>%
    plot_estimates("intensity", "abundance") +
    facet_wrap(~species, scales = "fixed") +
               #, scales = "free_y") +
    scale_fill_manual("",
                      values = cols) +
    scale_colour_manual("",
                        values = cols) +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "italic")),
  nrow = 2, labels = c("", ""),
  common.legend = TRUE, legend = "bottom")

ggsave("figures/main/figure3", device = "tiff", width = 7, height = 8, units = "in", dpi = 500) 

################################################################################
############# Figure 4

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
  data_rich <- rbind(data_rich, df) %>% subset(genus != "Total")
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
  data_abun <- rbind(data_abun, df) %>% subset(genus != "Total")
}

cols <- c("#5d1b3f", "#4099af")

p1 <- 
  data_rich %>% rbind(data_abun) %>% 
  ggplot(aes(genus, mean, colour = genus, group = model)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(position=position_dodge(width=0.4), aes(pch=model), alpha = 0.6, size = 3) +
  geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.4), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
  theme_light() + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(face = "italic")) +
  xlab("") + ylab("Estimate") +
    scale_fill_manual("",
                      values = c(cols[1], cols[2])) +
    scale_colour_manual("",
                        values = c(cols[1], cols[2])) +
  coord_flip() +
    guides(colour = "none") 

# Species abundance models
files <- list.files("models/abundance/species", pattern = "_l_deforestation_mod", full.names = TRUE)

data <- NULL
for (i in files){
  
  load(i)
  
  df <- extract_deforest_estimates(mod, "linear") %>% mutate(species = gsub("_l.*", "", gsub(".*species/", "", i))) %>%
    mutate(species = as.factor(str_to_sentence(gsub("_", " ", species))))
  data <- rbind(data, df)
}

# Order species
data$species <- fct_reorder(data$species, data$mean, min)

cols <- c("#802417", "#c06636", "#ce9344", "#e8b960", "#646e3b", "#2b5851", "#508ea2", "#17486f")
cols <- c("#bd3106", "#d9700e", "#e9a00e", "#eebe04", "#5b7314", "#c3d6ce", "#89a6bb", "#454b87")
  
p2 <- 
  data %>% 
  ggplot(aes(species, mean, colour = species, group = species)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(position=position_dodge(width=0.6), alpha = 0.6, size = 3) +
  geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
  theme_light() + 
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic")) +
  xlab("") + ylab("Estimate") +
  coord_flip() +
scale_fill_manual(values=cols) +
 scale_colour_manual(values=cols)


ggarrange(p1, p2, nrow = 2,
          labels = c("A", "B"))
ggsave("figures/main/figure4", device = "tiff", width = 6, height = 6.5, units = "in", dpi = 500)


################################################################################################################################################################
# Appendix

################################################################################
############# Figure S2
# Model residuals
load("functions/extract_residuals.RData")

# Plot model residuals
# Abundance
files <- list.files("models/abundance", pattern = "_lui_mod", full.names = TRUE)

residuals <- NULL
for (i in files){
  load(i)
  
  if (grepl("anopheles", i) == TRUE){
    genus <- "Anopheles"
  }else{if(grepl("aedes", i)==TRUE){
    genus <- "Aedes"
  }else{
    genus <- "Total"
  }
  }
  
  residuals <- rbind(residuals,
                     extract_residuals(mod, genus, "abundance"))
}

residuals$response_fitted[residuals$response_fitted == "-Inf"] <- 0
residuals$response_value[residuals$response_value == "-Inf"] <- 0

# lm for fit line
fit <- lm(residuals$response_fitted ~ residuals$response_value)

abun_residuals <- 
  residuals %>% mutate(genus = paste(genus, "abundance", sep = " ")) %>% 
  ggplot(aes(response_value, response_fitted)) +
  geom_point() +
  theme_light() +
  xlab("Observed") +
  ylab("Fitted") +
  facet_wrap(~genus) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "grey40")) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], col='darkred') 

# Species richness
files <- list.files("models/richness", full.names = TRUE, pattern = "_lui_mod")

residuals <- NULL
for (i in files){
  load(i)
  
  if (grepl("anopheles", i) == TRUE){
    genus <- "Anopheles"
  }else{if(grepl("aedes", i)==TRUE){
    genus <- "Aedes"
  }else{
    genus <- "Total"
  }
  }
  
  residuals <- rbind(residuals,
                     extract_residuals(mod, genus, "richness"))
}

rich_residuals <- residuals %>% mutate(genus = paste(genus, "richness", sep = " ")) %>%
  ggplot(aes(response_value, response_fitted)) +
  geom_point() +
  theme_light() +
  xlab("Observed") +
  ylab("Fitted") +
  facet_wrap(~genus) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "grey40")) +
  geom_abline(intercept = 0, slope = 1, col='darkred') 

ggpubr::ggarrange(abun_residuals,
                  rich_residuals,
                  labels = c("A", "B"),
                  nrow = 2)
ggsave("figures/supplementary/figureS2", device = "tiff", width = 6, height = 5, units = "in", dpi = 500)

################################################################################
############# Figure S3
# Most common study countries
world <- ne_countries(scale = "medium", returnclass = "sf")

data <- read.csv("data/input/data.csv") %>%
  # Count studies per country
  dplyr::group_by(study_number) %>%
  dplyr::summarise(country = unique(country)) %>% 
  dplyr::group_by(country) %>%
  dplyr::count() %>% mutate(n = as.numeric(n)) %>%
  dplyr::rename(admin = country) %>%
  # recode French Guiana to match map
  mutate(admin = recode(admin, "French Guiana" = "France"))

world_data <- merge(world, data, by = c("admin"), all = TRUE)

# Add coordinates for labels
world_data$x <- 0
world_data$y <- 0
world_data$x[world_data$admin == "Brazil"] <- -55
world_data$y[world_data$admin == "Brazil"] <- -9
world_data$x[world_data$admin == "Ecuador"] <- -85
world_data$y[world_data$admin == "Ecuador"] <- -0.5
world_data$x[world_data$admin == "France"] <- -51.5
world_data$y[world_data$admin == "France"] <- 8
world_data$x[world_data$admin == "Peru"] <- -75
world_data$y[world_data$admin == "Peru"] <- -10
world_data$x[world_data$admin == "Venezuela"] <- -65
world_data$y[world_data$admin == "Venezuela"] <- 7
world_data$x[world_data$admin == "Mexico"] <- -100
world_data$y[world_data$admin == "Mexico"] <- 20
world_data$x[world_data$admin == "Colombia"] <- -74
world_data$y[world_data$admin == "Colombia"] <- 3
world_data$x[world_data$admin == "Bolivia"] <- -64
world_data$y[world_data$admin == "Bolivia"] <- -17
world_data$x[world_data$admin == "Argentina"] <- -64
world_data$y[world_data$admin == "Argentina"] <- -35
world_data$x[world_data$admin == "Puerto Rico"] <- -62
world_data$y[world_data$admin == "Puerto Rico"] <- 22
world_data$x[world_data$admin == "Costa Rica"] <- -91
world_data$y[world_data$admin == "Costa Rica"] <- 8
world_data$x[world_data$admin == "Trinidad and Tobago"] <- -71
world_data$y[world_data$admin == "Trinidad and Tobago"] <- 15
world_data$x[world_data$admin == "Saint Kitts and Nevis"] <- -51
world_data$y[world_data$admin == "Saint Kitts and Nevis"] <- 16

world_data <- world_data %>% mutate(label = paste0(admin, "\n", n))

ggplot(data = world_data) +
  geom_sf() +
  geom_sf(aes(fill = n), size = 0.2) +
  coord_sf(xlim = c(-120, -30), ylim = c(-57, 35), expand = FALSE) +
  theme_light() +
  xlab("") + ylab("") +
  theme(plot.subtitle = element_text(face = "bold", size = 14)) +
  scale_fill_gradient2(name = "Number of\nstudies",
                       low = "#F1E4EB", mid = "#CDA2B8", high = "#693A48",
                       na.value = "grey80", midpoint = 14,
                       guide = guide_colourbar(ticks = FALSE)) +
  geom_label(data = world_data, aes(x,y, 
                                    label = gsub("France", "French Guiana", label)),
             size = 2.5, label.padding = unit(0.15, "lines"), fill = "white")
ggsave("figures/supplementary/figureS3", device = "tiff", width = 7, height = 7, units = "in", dpi = 500)

################################################################################
############# Figure S4 - total abundance and richness responses

# Total species richness
files <- list.files("models/richness", pattern = "_lui_mod", full.names = TRUE)
load(files[3])

data <- extract_estimates2(mod, "intensity")

# Set primary minimal category as baseline
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 

# Difference from intercept
data[,2:4] = (exp(data[2:4]) - 1)*100

cols <- c("#479A97", "#a37800", "#ffbe0c", "#6741ad")

# Plot
richness_plot <- plot_estimates(data, "intensity", "richness") +
  scale_fill_manual("",
                    values = cols) +
  scale_colour_manual("",
                      values = cols) +
  ylim(-60, 30)

# Total species abundance
files <- list.files("models/abundance", pattern = "_lui", full.names = TRUE)
load(files[5])

data <- extract_estimates2(mod, "intensity")

# Set primary minimal category as baseline
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 

# Difference from intercept
data[,2:4] = (exp(data[2:4]) - 1)*100

cols <- c("#479A97", "#a37800", "#ffbe0c", "#6741ad")

# Plot
abundance_plot <- plot_estimates(data, "intensity", "abundance") +
  scale_fill_manual("",
                    values = cols) +
  scale_colour_manual("",
                      values = cols) +
  ylim(-60, 30)

ggpubr::ggarrange(richness_plot,
                  abundance_plot,
                  common.legend = TRUE,
                  legend = "bottom",
                  labels = c("A", "B"))
ggsave("figures/supplementary/figureS4", device = "tiff", width = 8.5, height = 4, units = "in", dpi = 500) 

################################################################################
############# Figure S5
# Geographical cross-validation
##### Aedes
# Abundance
files <- list.files("models/abundance/cross_validation", pattern = "brazil_lui", full.names = TRUE)
load(files[1])
data <- extract_estimates2(mod, "intensity") %>% mutate(model = "without Brazil")

# Add model not Brazil data
load("models/abundance/aedes_lui_mod.RData")
df <- extract_estimates2(mod, "intensity") %>% mutate(model = "with Brazil") 
data <- data %>% rbind(df)

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

aedes_abun_plot <- data %>% 
  plot_estimates_cross("intensity", "abundance") +
  scale_fill_manual("",
                    values = c("white", "grey70")) +
  scale_colour_manual("",
                      values = c("black", "grey70")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic")) +
    labs(subtitle = "C") +
  guides(shape="none")

# Richness
files <- list.files("models/richness/cross_validation", pattern = "brazil_lui", full.names = TRUE)
load(files[1])
data <- extract_estimates2(mod, "intensity") %>% mutate(model = "without Brazil")

# Add model not Brazil data
load("models/richness/aedes_lui_mod.RData")
df <- extract_estimates2(mod, "intensity") %>% mutate(model = "with Brazil") 
data <- data %>% rbind(df)

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

aedes_rich_plot <- data %>% 
  plot_estimates_cross("intensity", "richness") +
  scale_fill_manual("",
                    values = c("white", "grey70")) +
  scale_colour_manual("",
                      values = c("black", "grey70")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic")) +
  labs(subtitle = "A") +
  guides(shape="none")

##### Anopheles
# Abundance
files <- list.files("models/abundance/cross_validation", pattern = "brazil_lui", full.names = TRUE)
load(files[2])
data <- extract_estimates2(mod, "intensity") %>% mutate(model = "without Brazil")

# Add model not Brazil data
load("models/abundance/anopheles_lui_mod.RData")
df <- extract_estimates2(mod, "intensity") %>% mutate(model = "with Brazil") 
data <- data %>% rbind(df)

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

anopheles_abun_plot <- data %>% 
  plot_estimates_cross("intensity", "abundance") +
  scale_fill_manual("",
                    values = c("white", "grey70")) +
  scale_colour_manual("",
                      values = c("black", "grey70")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic")) +
  labs(subtitle = "D") +
  guides(shape="none")

# Richness
files <- list.files("models/richness/cross_validation", pattern = "brazil_lui", full.names = TRUE)
load(files[2])
data <- extract_estimates2(mod, "intensity") %>% mutate(model = "without Brazil")

# Add model not Brazil data
load("models/richness/anopheles_lui_mod.RData")
df <- extract_estimates2(mod, "intensity") %>% mutate(model = "with Brazil") 
data <- data %>% rbind(df)

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

anopheles_rich_plot <- data %>% 
  plot_estimates_cross("intensity", "richness") +
  scale_fill_manual("",
                    values = c("white", "grey70")) +
  scale_colour_manual("",
                      values = c("black", "grey70")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic")) +
  labs(subtitle = "B") +
  guides(shape="none")

ggpubr::ggarrange(aedes_rich_plot, anopheles_rich_plot,
                  aedes_abun_plot, anopheles_abun_plot,
                  nrow = 2, ncol = 2,
                  common.legend = TRUE,
                  legend = "bottom")
ggsave("figures/supplementary/figureS5", device = "tiff", width = 7, height = 6, units = "in", dpi = 500)  

################################################################################
############# Figure S6
# Biome sensitivity analysis
# Anopheles and Aedes abundance
files <- list.files("models/abundance/cross_validation/biome", pattern = "anopheles", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/abundance/cross_validation/biome/anopheles_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

anopheles_abundance_plot <-
  data %>% mutate(model = gsub("_", " ", model)) %>%
  plot_estimates_cross("intensity", "abundance") +
  labs(subtitle = "D") +
  theme(plot.subtitle = element_text(face="bold")) +
  guides(colour = guide_legend(nrow = 6),
         shape = guide_legend(nrow = 3)) +
  scale_colour_discrete(labels = function(x) str_wrap(x, width = 20))


files <- list.files("models/abundance/cross_validation/biome", pattern = "aedes", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/abundance/cross_validation/biome/aedes_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

aedes_abundance_plot <-
  data %>% mutate(model = gsub("_", " ", model)) %>%
  plot_estimates_cross("intensity", "abundance") +
  labs(subtitle = "C") +
  theme(plot.subtitle = element_text(face="bold")) +
  guides(colour = guide_legend(nrow = 6),
         shape = guide_legend(nrow = 3)) +
  scale_colour_discrete(labels = function(x) str_wrap(x, width = 20))


# Anopheles and Aedes richness
files <- list.files("models/richness/cross_validation/biome", pattern = "anopheles", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/richness/cross_validation/biome/anopheles_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

anopheles_richness_plot <-
  data %>% mutate(model = gsub("_", " ", model)) %>%
  plot_estimates_cross("intensity", "richness") +
  labs(subtitle = "B") +
  theme(plot.subtitle = element_text(face="bold"))+
  guides(colour = guide_legend(nrow = 6),
         shape = guide_legend(nrow = 3)) +
  scale_colour_discrete(labels = function(x) str_wrap(x, width = 20))


files <- list.files("models/richness/cross_validation/biome", pattern = "aedes", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/richness/cross_validation/biome/aedes_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

aedes_richness_plot <-
  data %>% mutate(model = gsub("_", " ", model)) %>%
  plot_estimates_cross("intensity", "richness") +
  labs(subtitle = "A") +
  theme(plot.subtitle = element_text(face="bold")) +
  guides(colour = guide_legend(nrow = 6),
         shape = guide_legend(nrow = 3)) +
  scale_colour_discrete(labels = function(x) str_wrap(x, width = 20))

# Combine
biome_sensitivity <- 
  ggpubr::ggarrange(aedes_richness_plot, anopheles_richness_plot,
                    aedes_abundance_plot, anopheles_abundance_plot,
                    common.legend = TRUE,
                    legend = "right")

# map of sites
biome_shp <- rgdal::readOGR("data/wwf_biome", "wwf_terr_ecos")
definitions <- read.csv("data/wwf_biome/biome_definitions.csv") %>%
  dplyr::rename(BIOME = biome)
biome_data <- merge(biome_shp, definitions, by = c("BIOME"))

# Plot study sites
sp_data <- read.csv("data/inla_input/abundance.csv") %>% 
dplyr::select(lat, lon, land_use, land_use_intensity_agg) %>%
  unique() %>%
  mutate(land_use = factor(land_use, levels = c("primary vegetation", "secondary vegetation", "managed", "urban")),
         land_use_intensity_agg = factor(land_use_intensity_agg, 
                                         levels = c("minimal", "substantial"))) 

map <- merge(biome_shp, biome_data, by = c("BIOME")) %>%
  subset(BIOME == 1 | BIOME == 2 | BIOME == 3 | BIOME == 4 |
           BIOME == 5 | BIOME == 6 | BIOME == 7 | BIOME == 8 |
           BIOME == 9 | BIOME == 10 | BIOME == 11 | BIOME == 12 |
           BIOME == 13 | BIOME == 14) %>% 
  st_as_sf() %>%
  ggplot() + geom_sf(aes(fill = as.factor(definition)), size = 0,
                     colour = "transparent") +
  theme_light() +
  labs(subtitle = "E") + 
  theme(legend.position = "bottom",
        legend.title = element_blank(), 
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.2, "cm"),
        plot.subtitle = element_text(face="bold", hjust = -0.3)) +
  geom_point(aes(x = lon, y = lat), data = sp_data, colour = "red", size=0.45) +
  coord_sf(xlim = c(-120, -30), ylim = c(-57, 35), expand = FALSE) +
  guides(fill=guide_legend(nrow=5,byrow=TRUE)) +
  scale_fill_manual("Biome",
                    values = c("#B3D4B7", "#FCD57A", "#5BA187",
                               "#4099AF", "#D16E3F", "#B98320", "#A5C790",
                               "#75A95E", "#97B669", "#317A22",
                               "#2D5823", "#628746", "#304E52",
                               "#C1E1DD"),
                    labels = function(x) str_wrap(x, width = 50))

ggpubr::ggarrange(biome_sensitivity,
                  map,
                  nrow = 2)
ggsave("figures/supplementary/figureS6", device = "tiff", width = 8, height = 10, units = "in", dpi = 500)

################################################################################
############# Figure S7
# K-fold sensitivity analysis

# Anopheles and Aedes abundance
files <- list.files("models/abundance/cross_validation/kfold", pattern = "anopheles", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/abundance/cross_validation/kfold/anopheles_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

anopheles_abundance_plot <-
  data %>% 
  plot_estimates_cross("intensity", "abundance") +
  labs(subtitle = "D") +
  theme(plot.subtitle = element_text(face="bold"))

files <- list.files("models/abundance/cross_validation/kfold", pattern = "aedes", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/abundance/cross_validation/kfold/aedes_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

aedes_abundance_plot <-
  data %>% 
  plot_estimates_cross("intensity", "abundance") +
  labs(subtitle = "C") +
  theme(plot.subtitle = element_text(face="bold"))

# Anopheles and Aedes richness
files <- list.files("models/richness/cross_validation/kfold", pattern = "anopheles", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/richness/cross_validation/kfold/anopheles_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

anopheles_richness_plot <-
  data %>% 
  plot_estimates_cross("intensity", "richness") +
  labs(subtitle = "B") +
  theme(plot.subtitle = element_text(face="bold"))

files <- list.files("models/richness/cross_validation/kfold", pattern = "aedes", full.names = TRUE)

data <- NULL
for (i in files){
  load(i)
  
  model_name <- gsub("_lui_mod.RData", "", gsub("models/richness/cross_validation/kfold/aedes_", "", i))
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% mutate(model = model_name))
  
}

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

aedes_richness_plot <-
  data %>% 
  plot_estimates_cross("intensity", "richness") +
  labs(subtitle = "A") +
  theme(plot.subtitle = element_text(face="bold"))

# Combine
ggpubr::ggarrange(aedes_richness_plot, anopheles_richness_plot,
                  aedes_abundance_plot, anopheles_abundance_plot,
                  common.legend = TRUE,
                  legend = "bottom")
ggsave("figures/supplementary/figureS7", device = "tiff", width = 8, height = 6, units = "in", dpi = 500)


################################################################################
############# Figure S8

# Abundance without influential species
files <- list.files("models/abundance/cross_validation/species", full.names = TRUE)

data <- NULL
for (i in files[1:8]){
  
  load(i)
  
  # Create df of estimates
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% 
                  mutate(species = gsub("_", " ", gsub("_lui_mod.RData", "" , gsub("models/abundance/cross_validation/species/", "", i))),
                         model = "without")) 
}

# Add model not omitting species
load("models/abundance/aedes_lui_mod.RData")
df <- extract_estimates2(mod, "intensity") %>% mutate(model = "with") 
df <- do.call("rbind", replicate(8, df, simplify = FALSE)) %>% mutate(species = c(rep(c("aedes aegypti", "aedes albopictus",
                                                                                      "aedes scapularis", "aedes serratus"), each = 5),
                                                                      rep(c("anopheles albimanus", "anopheles albitarsis",
                                                                            "anopheles darlingi", "anopheles nuneztovari"), each = 5)))
data <- data %>% rbind(data, df)

# Set primary min to 0
data[data$land_use == "primary vegetation-minimal", 2:4] <- 0 # intercept

# difference from intercept
data[ , 2:4 ] = (exp(data[ 2:4]) - 1)*100

data %>% mutate(g = gsub(" .*", "", species)) %>% 
  mutate(species = str_to_sentence(species)) %>%
  plot_estimates_cross("intensity", "abundance") +
  facet_wrap(~species, scales = "free_y", ncol=2) +
  scale_fill_manual("",
                    values = c("white", "grey70")) +
  scale_colour_manual("",
                      values = c("black", "grey70")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic"))
ggsave("figures/supplementary/figureS8", device = "tiff", width = 7, height = 8, units = "in", dpi = 500)  




