
# Run models of genus and species-level abundance and richness models across land-use change

################################################################################
pacman::p_load("dplyr", "INLA")

# Load functions
load("functions/aggregate_lui.RData")
load("functions/run_mod.RData")

# Read in data, reset factor levels and aggregate LUI classes
data_abun <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all")
data_rich <- read.csv("data/inla_input/richness.csv") %>% aggregate_lui("all")
data_rich_ano <- read.csv("data/inla_input/ano_richness.csv") %>% aggregate_lui("all")
data_rich_aed <- read.csv("data/inla_input/aed_richness.csv") %>% aggregate_lui("all")

################################################################################
# Species abundance - Land use intensity
################################################################################
# Total species abundance
df_inla <- data_abun %>%
  mutate(abundance = measurement_adj) %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
  dplyr::summarise(abundance = sum(abundance))

length(unique(df_inla$study_number))
length(unique(df_inla$site_number))

run_mod(df_inla=df_inla, 
        land_use_var = "intensity", 
        response="abundance", 
        filename="total",
        species="Y")

######################################################################
# Anopheles and Aedes abundance

for (i in c("Anopheles", "Aedes")){
  
  if (i == "Anopheles"){
    df_inla <- data_abun %>% subset(genus == "Anopheles") %>%
      mutate(abundance = measurement_adj) %>%
      dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
      dplyr::summarise(abundance = sum(abundance))
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{
    df_inla <- data_abun %>% subset(genus == "Aedes") %>%
      mutate(abundance = measurement_adj) %>%
      dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
      dplyr::summarise(abundance = sum(abundance))
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0(tolower(i)),
          species="Y")
}

######################################################################
# Species richness - Land use intensity
######################################################################
# Total species richness
df_inla <- data_rich

length(unique(df_inla$study_number))
length(unique(df_inla$site_number))

run_mod(df_inla=df_inla, 
        land_use_var = "intensity", 
        response="richness", 
        filename="total",
        species="N")

######################################################################
# Anopheles and Aedes species richness

for (i in c("Anopheles", "Aedes")){
  
  if (i == "Anopheles"){
    df_inla <- data_rich_ano
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{
    df_inla <- data_rich_aed
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0(tolower(i)),
          species="N")
}

################################################################################
# Species-level abundance models
# Read in data and reset factor levels
# Read in data, reset factor levels and aggregate classes
data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") %>% mutate(species_name = tolower(species_name)) 

vector_list <- read.csv("data/vector_list/input/anopheles_sinka.csv") %>%
  mutate(Species_name = tolower(Species_name))

# Select species with most mentions
anopheles <- data %>% subset(species_name %in% vector_list$Species_name) %>%
  #dplyr::group_by(study_number) %>%
  #dplyr::summarise(species = unique(species_name)) %>% 
  dplyr::group_by(species) %>%
  dplyr::count() %>% mutate(n = as.numeric(n)) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(-n) 

# select dominant malaria vectors -An. marajoara estimates have large uncertainty
anopheles <- rbind(anopheles[1,], anopheles[2,],
                   anopheles[3,], anopheles[4,]) %>%
  mutate(species_name = paste("anopheles", species),
         genus = "Anopheles")
# check if all four associated with human disease

aedes <- data %>% subset(genus == "Aedes") %>%
  dplyr::group_by(species_name) %>%
  dplyr::count() %>% mutate(n = as.numeric(n)) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(-n) %>%
  dplyr::slice(1:4) 

################################################################################
# Land use intensity
################################################################################

for (i in c(aedes$species_name, anopheles$species_name)){
  
  for (j in unique(data$study_number)){
    
    df <- data %>% subset(species_name == i & study_number == j)
    
    # check studies do not include where species is not recorded
    if (max(df$measurement, na.ignore=T)<1){
      print(paste0("species not recorded", j))
    }
  }
  
  df_inla <- data %>% subset(species_name == i) %>% 
    mutate(abundance = measurement_adj) %>%
    dplyr::group_by(study_number, site_name, site_number, study_block, study_sample, LUI, biome) %>%
    dplyr::summarise(abundance = sum(abundance))
  
  length(unique(df_inla$study_number))
  length(unique(df_inla$site_number))
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("species/", gsub(" ", "_", i)),
          species="N")
}
