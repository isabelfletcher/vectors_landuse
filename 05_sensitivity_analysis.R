
# Run sensitivity analysis of abundance and species richness models, by leaving out influential data points

################################################################################
pacman::p_load("dplyr", "INLA", "sp", "rgdal", "dismo")

# Load functions
load("functions/aggregate_lui.RData")
load("functions/run_mod.RData")
load("functions/create_richness_df.RData")

# Read in data, reset factor levels and aggregate classes
data_abun <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all")
data_rich <- read.csv("data/inla_input/richness.csv") %>% aggregate_lui("all")
data_rich_ano <- read.csv("data/inla_input/ano_richness.csv") %>% aggregate_lui("all")
data_rich_aed <- read.csv("data/inla_input/aed_richness.csv") %>% aggregate_lui("all")

################################################################################
# Formulate a model without Brazil

# Species abundance 
################################################################################

# Total species abundance
df_inla <- data_abun %>% subset(country != "Brazil") %>% 
  mutate(abundance = measurement_adj) 

length(unique(df_inla$study_number))
length(unique(df_inla$site_number))

run_mod(df_inla=df_inla, 
        land_use_var = "intensity", 
        response="abundance", 
        filename="cross_validation/total_w_brazil",
        species="Y")

######################################################################
# Anopheles and Aedes abundance

for (i in c("Anopheles", "Aedes")){
  
  if (i == "Anopheles"){
    df_inla <- data_abun %>% subset(genus == "Anopheles") %>%
      subset(country != "Brazil") %>% 
      mutate(abundance = measurement_adj)
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{
    df_inla <- data_abun %>% subset(genus == "Aedes") %>%
      subset(country != "Brazil") %>% 
      mutate(abundance = measurement_adj)
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("cross_validation/", tolower(i), "_w_brazil"),
          species="Y")
}


######################################################################
# Species richness 
######################################################################
# Total species richness
df_inla <- data_rich %>% subset(country != "Brazil")

length(unique(df_inla$study_number))
length(unique(df_inla$site_number))

run_mod(df_inla=df_inla, 
        land_use_var = "intensity", 
        response="richness", 
        filename="cross_validation/total_w_brazil",
        species="N")

######################################################################
# Anopheles and Aedes species richness

for (i in c("Anopheles", "Aedes")){
  
  if (i == "Anopheles"){
    df_inla <- data_rich_ano %>% subset(country != "Brazil")
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{ 
    df_inla <- data_rich_aed %>% subset(country != "Brazil")
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0("cross_validation/", tolower(i), "_w_brazil"),
          species="N")
}

################################################################################
# We then formulate a model without influential species 

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
  
################################### Abundance
# Total abundance
for (i in c(anopheles$species_name, aedes$species_name)){
  
  # exclude species from model
  df_inla <- data_abun %>% subset(species_name != str_to_sentence(i)) %>% 
    mutate(abundance = measurement_adj)
  
  length(unique(df_inla$study_number))
  length(unique(df_inla$site_number))
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("cross_validation/species/total_", tolower(gsub(" ", "_", i))),
          species="Y")
}

# Anopheles and Aedes abundance
for (i in c(anopheles$species_name, aedes$species_name)){
  
  if (grepl("anopheles", i)==TRUE){
    df_inla <- data_abun %>% subset(genus == "Anopheles") %>%
      subset(species_name != str_to_sentence(i)) %>% 
      mutate(abundance = measurement_adj)
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{
    df_inla <- data_abun %>% subset(genus == "Aedes") %>%
      subset(species_name != str_to_sentence(i)) %>% 
      mutate(abundance = measurement_adj)
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("cross_validation/species/", #gsub("[ ].*", "_", tolower(i)), 
                          gsub(" ", "_", tolower(i))),
          species="Y")
}

################################### Richness
# Total richness
for (i in c(anopheles$species_name, aedes$species_name)){
  
  # exclude species from model
  df_inla <- create_richness_df(i,"total")
  
  length(unique(df_inla$study_number))
  length(unique(df_inla$site_number))
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0("cross_validation/species/total_", tolower(gsub(" ", "_", i))),
          species="N")
}

# Anopheles and Aedes richness
for (i in c(anopheles$species_name, aedes$species_name)){
  
  if(grepl("anopheles", i)){
    df_inla <- create_richness_df(i, "Anopheles")
  } else {
    df_inla <- create_richness_df(i, "Aedes")
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0("cross_validation/species/", gsub("[ ].*", "_", tolower(i)), gsub(" ", "_", tolower(i))),
          species="N")
}



################################################################################
# # Formulate a model without each biome
data_abun <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all")


biomes <- c(unique(data_abun$biome)[1], unique(data_abun$biome)[2], unique(data_abun$biome)[3],
            unique(data_abun$biome)[5], unique(data_abun$biome)[6], unique(data_abun$biome)[7])

# Total abundance
##########################################

for (i in 1:length(biomes)){
  
  df_inla <- data_abun %>% subset(biome != biomes[i]) %>% 
    mutate(abundance = measurement_adj)
  
  length(unique(df_inla$study_number))
  length(unique(df_inla$site_number))
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("cross_validation/biome/total_", tolower(gsub(" ", "_", biomes[i]))),
          species="Y")
}

######################################################################
# Anopheles and Aedes abundance

for (i in 1:length(biomes)){
  
  for (j in c("Anopheles", "Aedes")){
  
  if (grepl("Anopheles", j)==TRUE){
    df_inla <- data_abun %>% subset(genus == "Anopheles") %>%
      subset(biome != biomes[i]) %>% 
      mutate(abundance = measurement_adj)
    
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  } else{
    df_inla <- data_abun %>% subset(genus == "Aedes") %>%
      subset(biome != biomes[i]) %>% 
      mutate(abundance = measurement_adj)
    
    
    length(unique(df_inla$study_number))
    length(unique(df_inla$site_number))
    
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("cross_validation/biome/", tolower(j), "_", tolower(gsub(" ", "_", biomes[i]))),
          species="Y")
  
  }
  
}

################################### Richness

# Total species richness
for (i in 1:length(biomes)){
  
  df_inla <- data_rich %>% subset(biome != biomes[i]) 
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0("cross_validation/biome/total_", tolower(gsub(" ", "_", biomes[i]))),
          species="N")
}

######################################################################
# Anopheles and Aedes species richness

for (i in 1:length(biomes)){
  
 for (j in c("Anopheles", "Aedes")){
   
  if (grepl("Anopheles", j)==TRUE){
    df_inla <- data_rich_ano %>% subset(biome != biomes[i]) 
  } else{
    df_inla <- data_rich_aed %>% subset(biome != biomes[i]) 
  }
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0("cross_validation/biome/", tolower(j), "_", tolower(gsub(" ", "_", biomes[i]))),
          species="N")
 }
  
}

################################################################################
# Formulate a model holding out 12.5% of data at a time 
# Randomly remove studies to test sensitivity
# Excludes 12.5% of studies per iteration
library(dismo)

studies <- unique(data_abun$study_number)
kfold  <- data.frame(study_number = studies, kfold = kfold(studies, 8)) %>% 
 write.csv("data/inla_input/kfold_23_04_01.csv", row.names =FALSE)

kfold <- read.csv("data/inla_input/kfold_23_04_01.csv")

# Abundance
##########################################
# Total abundance

for (i in 1:8){
  
  df_inla <- data_abun %>% merge(kfold, by = "study_number") %>% 
    subset(kfold != i) %>% 
    mutate(abundance = measurement_adj)
  
  length(unique(df_inla$study_number))
  length(unique(df_inla$site_number))
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0("cross_validation/kfold/total_", i),
          species="Y")
}

######################################################################
# Anopheles and Aedes abundance

for (i in 1:8){
  
  for (j in c("Anopheles", "Aedes")){
    
    if (grepl("Anopheles", j)==TRUE){
      df_inla <- data_abun %>% subset(genus == "Anopheles") %>%
        merge(kfold, by = "study_number") %>% 
        subset(kfold != i) %>% 
        mutate(abundance = measurement_adj)
      
      
      length(unique(df_inla$study_number))
      length(unique(df_inla$site_number))
      
    } else{
      df_inla <- data_abun %>% subset(genus == "Aedes") %>%
        merge(kfold, by = "study_number") %>% 
        subset(kfold != i) %>% 
        mutate(abundance = measurement_adj)
      
      
      length(unique(df_inla$study_number))
      length(unique(df_inla$site_number))
      
    }
    
    run_mod(df_inla=df_inla, 
            land_use_var = "intensity", 
            response="abundance", 
            filename=paste0("cross_validation/kfold/", tolower(j), "_", i),
            species="Y")
    
  }
  
}

################################### Richness

# Total species richness
for (i in 1:8){
  
  df_inla <- data_rich %>% merge(kfold, by = "study_number") %>% 
    subset(kfold != i) 
  
  run_mod(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0("cross_validation/kfold/total_", i),
          species="N")
}

######################################################################
# Anopheles and Aedes species richness

for (i in 1:8){
  
  for (j in c("Anopheles", "Aedes")){
    
    if (grepl("Anopheles", j)==TRUE){
      df_inla <- data_rich_ano %>% merge(kfold, by = "study_number") %>% 
        subset(kfold != i) 
    } else{
      df_inla <- data_rich_aed %>% merge(kfold, by = "study_number") %>% 
        subset(kfold != i) 
    }
    
    run_mod(df_inla=df_inla, 
            land_use_var = "intensity", 
            response="richness", 
            filename=paste0("cross_validation/kfold/", tolower(j), "_", i),
            species="N")
  }
  
}
