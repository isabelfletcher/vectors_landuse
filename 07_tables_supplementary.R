
# Supplementary tables

################################################################################

pacman::p_load("dplyr", "tibble")


# Load functions
load("functions/aggregate_lui.RData")

################################################################################
############# Table S4
# Land use categories used in models

# Read in data and reset factor levels
data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all")

# number of sites
sites <- data %>% dplyr::group_by(LUI) %>%
  dplyr::summarise(n = length(unique(site_number))) %>%
  arrange(LUI)

# number of records
records <- data %>% dplyr::group_by(LUI) %>%
  count() %>%
  arrange(LUI) 

merge(sites, records, by = c("LUI")) %>%
  arrange(LUI) %>% print(row.names = F)

################################################################################
############# Table S6
load("functions/aggregate_lui.RData")
# Abundance models
data <- read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>% aggregate_lui("all") %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>% dplyr::ungroup()

# Total abundance
# number of sites
data %>% 
  dplyr::summarise(n = length(unique(site_number)))

# number of records
data %>%
  count() 

# Aedes abundance
data <- read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>% aggregate_lui("all") %>%
  subset(genus == "Aedes") %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>% dplyr::ungroup()

# number of sites
data %>% 
  dplyr::summarise(n = length(unique(site_number)))

# number of records
data %>%
  count() 

# Anopheles abundance
data <- read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>% aggregate_lui("all") %>%
  subset(genus == "Anopheles") %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>% dplyr::ungroup()

# number of sites
data %>% 
  dplyr::summarise(n = length(unique(site_number)))

# number of records
data %>% 
  count() 

# Richness models
data <- read.csv("data/input/americas/inla_input/richness_unscaled_sampling_effort.csv") %>% aggregate_lui("all")

# Total richness
# number of sites
data %>% 
  dplyr::summarise(n = length(unique(site_number)))

# number of records
data %>%
  count() 

# Aedes richness
data <- read.csv("data/input/americas/inla_input/aed_richness_unscaled_sampling_effort.csv") %>% aggregate_lui("all")

# number of sites
data %>% 
  dplyr::summarise(n = length(unique(site_number)))

# number of records
data %>% 
  count() 

# Anopheles richness
data <- read.csv("data/input/americas/inla_input/ano_richness_unscaled_sampling_effort.csv") %>% aggregate_lui("all")
# number of sites
data %>% 
  dplyr::summarise(n = length(unique(site_number)))

# number of records
data %>% 
  count() 


################################################################################
############# Table S7

# Number of sites
# Species included in species-specific abundance models
# Total
read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>%
  dplyr::group_by(study_number, site_name, site_number, study_block, study_sample, LUI, genus, species) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(genus, species) %>%
  dplyr::summarise(n = length(unique(site_number))) %>%
  dplyr::arrange(genus, species) %>%
  subset(species == "albopictus" | species == "aegypti" |
           species == "serratus" | species == "scapularis" | 
           species == "albimanus" | species == "albitarsis" |
           species == "darlingi" | species == "nuneztovari")

# By LUI
read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>% aggregate_lui("all") %>%
  dplyr::group_by(study_number, site_name, site_number, study_block, study_sample, LUI, genus, species) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(genus, species, LUI) %>%
  dplyr::summarise(n = length(unique(site_number))) %>%
  dplyr::arrange(genus, species, LUI) %>%
  subset(species == "albopictus" | species == "aegypti" |
           species == "serratus" | species == "scapularis" | 
           species == "albimanus" | species == "albitarsis" |
           species == "darlingi" | species == "nuneztovari") %>%
  dplyr::ungroup() %>%
  subset(species == "aegypti") %>%
  dplyr::select(n) %>% 
  print.data.frame(row.names = F) 


# Number of records
# Total
read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>%
  dplyr::group_by(study_number, site_name, site_number, study_block, study_sample, LUI, genus, species) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(genus, species) %>%
  dplyr::count() %>%
  dplyr::arrange(genus, species) %>%
  subset(species == "albopictus" | species == "aegypti" |
           species == "serratus" | species == "scapularis" | 
           species == "albimanus" | species == "albitarsis" |
           species == "darlingi" | species == "nuneztovari")

# By LUI
read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>% aggregate_lui("all") %>%
  dplyr::group_by(study_number, site_name, site_number, study_block, study_sample, LUI, genus, species) %>%
  dplyr::summarise(abundance = sum(measurement_adj)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(genus, species, LUI) %>%
  dplyr::count() %>%
  dplyr::arrange(genus, species, LUI) %>%
  subset(species == "albopictus" | species == "aegypti" |
           species == "serratus" | species == "scapularis" | 
           species == "albimanus" | species == "albitarsis" |
           species == "darlingi" | species == "nuneztovari") %>%
  dplyr::ungroup() %>%
  subset(species == "nuneztovari") %>%
  dplyr::select(n) %>% 
  print.data.frame(row.names = F) 

################################################################################
############# Table S8
# Richness model selection
read.csv("figures/supplementary/richness_model_selection.csv") %>% 
  subset(Variable == "land use intensity") %>% subset(Genus == "Anopheles") %>%
  dplyr::select(WAIC) %>% 
  print.data.frame(row.names = F)

################################################################################
############# Table S9
# abundance model selection
read.csv("figures/supplementary/abundance_model_selection.csv") %>% 
  subset(Variable == "land use intensity") %>% subset(Genus == "Anopheles") %>%
  dplyr::select(WAIC) %>% 
  print.data.frame(row.names = F)

################################################################################
############# Table S10
# Deforestation models
# Total abundance
data <- read.csv("data/input/americas/inla_input/deforestation_data.csv") %>%
  mutate(abundance = measurement_adj) %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, biome) %>%
  dplyr::summarise(abundance = sum(abundance),
                   deforestation = unique(cells_deforestation))
# Sites
length(unique(data$site_number))
# Records
nrow(data)

# Aedes abundance
data <- read.csv("data/input/americas/inla_input/deforestation_data.csv") %>%
  mutate(abundance = measurement_adj) %>% subset(genus == "Aedes") %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, biome) %>%
  dplyr::summarise(abundance = sum(abundance),
                   deforestation = unique(cells_deforestation))
# Sites
length(unique(data$site_number))
# Records
nrow(data)

# Anopheles abundance
data <- read.csv("data/input/americas/inla_input/deforestation_data.csv") %>%
  mutate(abundance = measurement_adj) %>% subset(genus == "Anopheles") %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, biome) %>%
  dplyr::summarise(abundance = sum(abundance),
                   deforestation = unique(cells_deforestation))
# Sites
length(unique(data$site_number))
# Records
nrow(data)

# Total richness
load("functions/create_richness_df.RData")
data <- create_deforest_richness_df("total")
# Sites
length(unique(data$site_number))
# Records
nrow(data)

# Aedes richness
data <- create_deforest_richness_df("Aedes")
# Sites
length(unique(data$site_number))
# Records
nrow(data)

# Aedes richness
data <- create_deforest_richness_df("Anopheles")
# Sites
length(unique(data$site_number))
# Records
nrow(data)

################################################################################
############# Table S11
# Species included in abundance models - number of records per species
read.csv("data/input/americas/inla_input/abundance_unscaled_sampling_effort.csv") %>%
  dplyr::group_by(genus, species) %>%
  dplyr::count() %>%
  dplyr::arrange(genus, species) %>% as.data.frame() %>%
  #write.csv("figures/americas/supplementary/species_list_count_abun.csv", row.names = FALSE) %>%
  dplyr::select(n) %>% #subset(genus == "Aedes") %>%
  print.data.frame(row.names = F) #%>% nrow()

################################################################################
############# Table S12
# Parameter estimates-richness models
load("functions/extract_estimates2.RData")

# Total richness
files <- list.files("models/americas/richness", pattern = "v2_lui", full.names = TRUE)
load(files[3])

extract_estimates2(mod, "intensity") %>% dplyr::select(mean, land_use, lci, uci) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>% 
  print.data.frame(row.names = F)

# Aedes richness
files <- list.files("models/americas/richness", pattern = "v2_lui", full.names = TRUE)
load(files[1])

extract_estimates2(mod, "intensity") %>% dplyr::select(mean, land_use, lci, uci) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>% 
  print.data.frame(row.names = F)

# Anopheles richness
files <- list.files("models/americas/richness", pattern = "v2_lui", full.names = TRUE)
load(files[2])

extract_estimates2(mod, "intensity") %>% dplyr::select(mean, land_use, lci, uci) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>% 
  print.data.frame(row.names = F)

################################################################################
############# Table S13
# Parameter estimates-abundance models
load("functions/extract_estimates2.RData")

# Total abundance
files <- list.files("models/americas/abundance/site_level", pattern = "v2_lui", full.names = TRUE)
load(files[3])

extract_estimates2(mod, "intensity") %>% dplyr::select(mean, land_use, lci, uci) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>% 
  print.data.frame(row.names = F)

# Aedes abundance
files <- list.files("models/americas/abundance/site_level", pattern = "v2_lui", full.names = TRUE)
load(files[1])

extract_estimates2(mod, "intensity") %>% dplyr::select(mean, land_use, lci, uci) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>% 
  print.data.frame(row.names = F)

# Anopheles abundance
files <- list.files("models/americas/abundance/site_level", pattern = "v2_lui", full.names = TRUE)
load(files[2])

extract_estimates2(mod, "intensity") %>% dplyr::select(mean, land_use, lci, uci) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>% 
  print.data.frame(row.names = F)

################################################################################
############# Table S14
# Parameter estimates- species abundance models
files <- list.files("models/americas/abundance/species", pattern = "v2", full.names = TRUE)

data <- NULL
for (i in files){
  
  load(i)
  
  # Create df of estimates
  data <- rbind(data,
                extract_estimates2(mod, "intensity") %>% 
                  mutate(species = gsub("_", " ", gsub("_lui_v2_lui_mod.RData", "" , gsub("models/americas/abundance/species/site_level/", "", i)))))
}

data %>% dplyr::select(mean, land_use, lci, uci, species) %>%
  dplyr::mutate(mean = round(mean, 2),
                lci = round(lci, 2),
                uci = round(uci, 2)) %>%
  dplyr::select(uci) %>%
  print.data.frame(row.names = F)
