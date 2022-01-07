
# Run models of abundance, species richness and deforestation

################################################################################
pacman::p_load("dplyr", "INLA", "tibble", "ggplot2", "scales", "ggpubr")

# Load functions
load("functions/run_deforest_mod.RData")
load("functions/extract_deforest_estimates.RData")
load("functions/create_richness_df.RData")

# Read in data
data_abun <- read.csv("data/inla_input/deforestation_data.csv") %>%
  mutate(deforestation = deforest_prop)

################################################################################
# Species richness - Land use intensity
################################################################################
richness_models <- c("total", "Aedes", "Anopheles")

for (i in richness_models){
  
  df_inla <- create_deforest_richness_df(i) %>% mutate(model = i)
  
  df_inla$deforestation <- scale(df_inla$deforestation)[,1]

  run_deforest_mod(df_inla=df_inla, 
                   deforest_var = "linear", 
                   response="richness", 
                   filename=tolower(i),
                   species="N")
  
}

################################################################################
# Total abundance - Land use intensity
################################################################################
# Total species abundance
df_inla <- data_abun %>%
  mutate(abundance = measurement_adj) %>%
  dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, biome) %>%
  dplyr::summarise(abundance = sum(abundance),
                   deforestation = unique(cells_deforestation))

df_inla$deforestation <- scale(df_inla$deforestation)[,1]

run_deforest_mod(df_inla=df_inla, 
        deforest_var = "linear", 
        response="abundance", 
        filename="total",
        species="Y")


######################################################################
# Anopheles and Aedes abundance

for (i in c("Anopheles", "Aedes")){
  
  if (i == "Anopheles"){
    df_inla <- data_abun %>% subset(genus == "Anopheles") %>%
      mutate(abundance = measurement_adj) %>%
      dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, biome) %>%
      dplyr::summarise(abundance = sum(abundance),
                       deforestation = unique(cells_deforestation))
    
    df_inla$deforestation <- scale(df_inla$deforestation)[,1]
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{
    df_inla <- data_abun %>% subset(genus == "Aedes") %>%
      mutate(abundance = measurement_adj) %>%
      dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, biome) %>%
      dplyr::summarise(abundance = sum(abundance),
                       deforestation = unique(cells_deforestation))
    
    df_inla$deforestation <- scale(df_inla$deforestation)[,1]
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
  }
  
  run_deforest_mod(df_inla=df_inla, 
                   deforest_var = "linear", 
                   response="abundance", 
                   filename=paste0(tolower(i)),
                   species="Y")
}

######################################################################
# Species abundance
sort(unique(data_abun$species_name))

# Species with largest no. records
data_abun %>% 
  dplyr::group_by(species_name) %>%
  dplyr::count() %>% mutate(n = as.numeric(n)) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(-n) %>% as.data.frame() 

species_names <- c("Aedes aegypti", "Aedes albopictus", # species picked on basis of uncertainty
                   "Aedes scapularis", "Aedes serratus",
                    "Anopheles albitarsis", 
                   "Anopheles darlingi",
                   "Anopheles mattogrossensis",
                   "Anopheles nuneztovari") 

for (i in c(species_names)){
  
  for (j in unique(data$study_number)){
    
    df <- data_abun %>% subset(species_name == i & study_number == j)
    
    # check studies do not include where species is not recorded
    if (max(df$measurement, na.ignore=T)<1){
      print(paste0("species not recorded", j))
    }
  }
  
  df_inla <- data_abun %>% subset(species_name == i) %>% 
    mutate(abundance = measurement_adj) %>%
    dplyr::group_by(study_number, site_name, site_number, study_block, study_sample, biome) %>%
    dplyr::summarise(abundance = sum(abundance),
                     deforestation = unique(deforestation))
  
  df_inla$deforestation <- scale(df_inla$deforestation)[,1]
  
  run_deforest_mod(df_inla=df_inla, 
          deforest_var = "linear", 
          response="abundance", 
          filename=paste0("species/", gsub(" ", "_", tolower(i))),
          species="N")
}
