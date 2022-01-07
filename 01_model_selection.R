
# Run model selection procedure of genus and species-level abundance and richness models across land-use change

################################################################################
pacman::p_load("dplyr", "INLA")

# Load functions
load("functions/run_mod_selection.RData")
load("functions/aggregate_lui.RData")


# Read in data, reset factor levels and aggregate classes
data_abun <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all")
data_rich <- read.csv("data/inla_input/richness.csv") %>% aggregate_lui("all")
data_rich_ano <- read.csv("data/inla_input/ano_richness.csv") %>% aggregate_lui("all")
data_rich_aed <- read.csv("data/inla_input/aed_richness.csv") %>% aggregate_lui("all")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                  Species abundance - Land use intensity                ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#                       Total species abundance                ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_inla <- data_abun %>%
  mutate(abundance = measurement_adj) %>%
dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
  dplyr::summarise(abundance = sum(abundance))

length(unique(df_inla$study_number))
length(unique(df_inla$site_number))

run_mod_selection(df_inla=df_inla, 
        land_use_var = "intensity", 
        response="abundance", 
        filename="total",
        species="Y")

#                   Anopheles and Aedes abundance              ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  
  run_mod_selection(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="abundance", 
          filename=paste0(tolower(i)),
          species="Y")
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                    Species richness - Land use intensity                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#                       Total species richness                 ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_inla <- data_rich

length(unique(df_inla$study_number))
length(unique(df_inla$site_number))

run_mod_selection(df_inla=df_inla, 
        land_use_var = "intensity", 
        response="richness", 
        filename="total",
        species="N")

#                    Anopheles and Aedes richness              ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  
  run_mod_selection(df_inla=df_inla, 
          land_use_var = "intensity", 
          response="richness", 
          filename=paste0(tolower(i)),
          species="N")
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Model table                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#                           Abundance models                   ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
files <- list.files("models/abundance/model_selection", pattern = "*.RData", full.names = TRUE)

abundance_models <- NULL
for (i in files){
  load(i)
  
  if (grepl("lu_mod", i)==TRUE){
    variable <- "land use"
  } else {
    variable <- "land use intensity"
  }
  
  if (grepl("anopheles", i)==TRUE){
    genus <- "Anopheles"
  } else {
    if (grepl("aedes", i)==TRUE){
      genus <- "Aedes"
    } else {
      genus <- "Total"
    }
  }
  
  model    <- as.numeric(gsub(".*?([0-9]+).*", "\\1", i))
  
  abundance_models <- rbind(abundance_models,
                    data.frame(Variable = variable,
                               Genus    = genus,
                               Model    = model,
                               DIC      = mod$dic$dic,
                               WAIC     = mod$waic$waic))
}

abundance_models %>% dplyr::arrange(Variable,Genus, Model) %>%
  mutate(DIC = round(DIC, 2),
         WAIC = round(WAIC, 2)) %>%
  write.csv("figures/supplementary/abundance_model_selection.csv", row.names = FALSE)


#                          Richness models                     ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
files <- list.files("models/richness/model_selection", pattern = "*.RData", full.names = TRUE)

richness_models <- NULL
for (i in files){
  load(i)
  
  if (grepl("lu_mod", i)==TRUE){
    variable <- "land use"
  } else {
    variable <- "land use intensity"
  }
  
  if (grepl("anopheles", i)==TRUE){
    genus <- "Anopheles"
  } else {
    if (grepl("aedes", i)==TRUE){
      genus <- "Aedes"
    } else {
      genus <- "Total"
    }
  }
  
  model    <- as.numeric(gsub(".*?([0-9]+).*", "\\1", i))
  
  richness_models <- rbind(richness_models,
                            data.frame(Variable = variable,
                                       Genus    = genus,
                                       Model    = model,
                                       DIC      = mod$dic$dic,
                                       WAIC     = mod$waic$waic))
}

richness_models %>% dplyr::arrange(Variable,Genus, Model) %>%
  mutate(DIC = round(DIC, 2),
         WAIC = round(WAIC, 2)) %>%
  write.csv("figures/supplementary/richness_model_selection.csv", row.names = FALSE)

