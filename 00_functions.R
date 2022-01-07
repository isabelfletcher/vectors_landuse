## This script provides functions used for analysis and data visualisation

###########################################################################
# Set factor levels (for non-aggregated land use classes)
# data: data input
set_factor_levels <- function(data){
  
  data <- data %>% 
    mutate(land_use = factor(land_use, levels = c("primary vegetation", "secondary vegetation", "managed", "urban")),
           LUI      = factor(LUI, levels = c("primary vegetation-minimal", "primary vegetation-substantial", "secondary vegetation-minimal",
                                           "secondary vegetation-substantial", "managed-substantial", "urban-minimal", "urban-substantial")))
  
  return(data)
}
save(set_factor_levels, file = "functions/set_factor_levels.RData")

###########################################################################
# Aggregate land use classes
# data (df)
# all species or specific dataframe (chr)

aggregate_lui <- function(data, species){
  
  if(grepl("darlingi", species)==TRUE){
    
    data <- data %>%
      mutate(LU_cat = factor(LU_cat, levels=c("forest", "deforested/forest fringe", "managed", "urban")))
    
  } else{
  
    data$LUI[data$LUI == "urban-minimal"] <- "urban"
    data$LUI[data$LUI == "urban-substantial"] <- "urban"
    data$LUI[data$LUI == "secondary vegetation-minimal"] <- "secondary vegetation"
    data$LUI[data$LUI == "secondary vegetation-substantial"] <- "secondary vegetation"
    data$LUI[data$LUI == "managed-substantial"] <- "managed"
    data$LUI[data$LUI == "managed-light"] <- "managed"
    
    data <- data %>% mutate(LUI = factor(LUI, levels = c("primary vegetation-minimal", "primary vegetation-substantial", 
                                                                   "secondary vegetation", "managed", "urban")))
  }
  
  return(data)
}
save(aggregate_lui, file = "functions/aggregate_lui.RData")

###########################################################################
# Run INLA model selection
# df
# land use change variable (land use or land use intensity)
# abundance or richness
# to save file
# whether model includes a random effect term for species

run_mod_selection <- function(df_inla, land_use_var, response, filename, species){
  
  if(grepl("abundance", response)==TRUE){
    response <- "log(abundance+1)"
  }else{
    response <- "richness"
  }
  
  if(grepl("intensity", land_use_var)==TRUE){
    var <- "~ 1 + LUI +"
    mod_name <- "_lui"
  }else{
    var <- "~ 1 + land_use +"
    mod_name <- "_lu"
  }
  
  if (grepl("Y", species)==TRUE){
    formulae <- data.frame(model   = c(1:5),
                           formula = c("f(study_number, model = 'iid') + f(site_number, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid') + f(study_sample, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid') + f(study_sample, model = 'iid') + f(species_name, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid') + f(study_sample, model = 'iid') + f(species_name, model = 'iid') + f(biome, model = 'iid')"))         
  }else{
    # remove species from formula (for species-level models and richness models)
    formulae <- data.frame(model   = c(1:4),
                           formula = c("f(study_number, model = 'iid') + f(site_number, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid') + f(study_sample, model = 'iid')",
                                       "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_block, model = 'iid') + f(study_sample, model = 'iid') + f(biome, model = 'iid')"))         
    
  }
  
  for (i in 1:nrow(formulae)){
    
    formula <- formula(paste(response, var, 
                             formulae$formula[i]))
    
    if(grepl("abundance", response)==TRUE){
      mod <- inla(formula, data = df_inla,
                  family = "gaussian", 
                  verbose = FALSE,
                  control.predictor=list(compute = TRUE),
                  control.compute=list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
      save(mod, file = paste0("models/abundance/model_selection/", i, filename, mod_name, "_mod.RData"))
    }else{
      
      mod <- inla(formula, data = df_inla,
                  family = "poisson", E = 1,
                  verbose = FALSE,
                  control.predictor=list(compute = TRUE, link = 1),
                  #control.family = list(link = "log"),
                  control.compute=list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
      save(mod, file = paste0("models/richness/model_selection/", i, filename, mod_name, "_mod.RData"))
    }
    
    
  }
}
save(run_mod_selection, file = "functions/run_mod_selection.RData")


###########################################################################
# Run land use INLA model, based on best random effects structure 
# df_inla - dataframe with model data (df)
# land use var - land use or land use intensity (chr)
# response - abundance or richness (chr)
# filename - name and extended path where to save model (chr)
# species - include species in model? Y or N (chr)
run_mod <- function(df_inla, land_use_var, response, filename, species){
  
  if(grepl("abundance", response)==TRUE){
    response <- "log(abundance+1)"
  }else{
    response <- "richness"
  }
  
  if(grepl("intensity", land_use_var)==TRUE){
    var <- "~ 1 + LUI +"
    mod_name <- "_lui"
  }else{
    var <- "~ 1 + land_use +"
    mod_name <- "_lu"
  }
  
  if (grepl("Y", species)==TRUE){
  formula <- formula(paste(response, var, 
                           "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(species_name, model = 'iid') + f(study_sample, model = 'iid') + f(biome, model = 'iid')"))
  }else{
    # remove species from formula (for species-level models and richness models)
    formula <- formula(paste(response, var, 
                             "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_sample, model = 'iid') + f(biome, model = 'iid')"))
  }
  
  if(grepl("abundance", response)==TRUE){
    mod <- inla(formula, data = df_inla,
                family = "gaussian", 
                verbose = FALSE,
                control.predictor=list(compute = TRUE),
                control.compute=list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
    save(mod, file = paste0("models/abundance/", filename, mod_name, "_mod.RData"))
  }else{
    
    formula <- formula(paste(response, var, 
                             "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_sample, model = 'iid')"))
    
    mod <- inla(formula, data = df_inla,
                family = "poisson", E = 1,
                verbose = FALSE,
                control.predictor=list(compute = TRUE, link = 1),
                control.compute=list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
    save(mod, file = paste0("models/richness/", filename, mod_name, "_mod.RData"))
  }
  
}
save(run_mod, file = "functions/run_mod.RData")

###########################################################################
# Run deforestation INLA model
# df_inla - dataframe with model data (df)
# land use var - land use or land use intensity (chr)
# response - abundance or richness (chr)
# filename - name and extended path where to save model (chr)
# deforest_var - option to include deforestaton as linear or non linear variable (chr)
# species - include species in model? Y or N (chr)
run_deforest_mod <- function(df_inla, response, deforest_var, filename, species){
  
  if(grepl("abundance", response)==TRUE){
    response <- "log(abundance+1)"
  }else{
    response <- "richness"
  }
  
  if(grepl("non linear", deforest_var)==TRUE){
    var <- "~ 1 + f(deforestation, model = 'iid') +"
    mod_name <- "_nl_deforestation"
  }else{
    var <- "~ 1 + deforestation +"
    mod_name <- "_l_deforestation"
  }
  
  if (grepl("Y", species)==TRUE){
    formula <- formula(paste(response, var, 
                             "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(species_name, model = 'iid') + f(study_sample, model = 'iid') + f(biome, model = 'iid')"))
  }else{
    # remove species from formula (for species-level models and richness models)
    formula <- formula(paste(response, var, 
                             "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_sample, model = 'iid') + f(biome, model = 'iid')"))
  }
  
  if(grepl("abundance", response)==TRUE){
    mod <- inla(formula, data = df_inla,
                family = "gaussian", 
                verbose = FALSE,
                control.predictor=list(compute = TRUE),
                control.compute=list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
    save(mod, file = paste0("models/abundance/", filename, mod_name, "_mod.RData"))
  }else{
    
    formula <- formula(paste(response, var, 
                             "f(study_number, model = 'iid') + f(site_number, model = 'iid') + f(study_sample, model = 'iid')"))
    
    mod <- inla(formula, data = df_inla,
                family = "poisson", E = 1,
                verbose = FALSE,
                control.predictor=list(compute = TRUE, link = 1),
                #control.family = list(link = "log"),
                control.compute=list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
    save(mod, file = paste0("models/richness/", filename, mod_name, "_mod.RData"))
  }
  
}
save(run_deforest_mod, file = "functions/run_deforest_mod.RData")

###########################################################################

# Create df of species richness for cross-validation models, where some species are left out 

# species - which species to leave out (chr)
# genus  - total, aedes or anopheles species richness (chr)

create_richness_df <- function(species_out, genus){
  
  data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all")
  
  if(genus == "total"){
  
  richness <- data %>% 
    subset(species_studied == "multiple") %>% 
    mutate(species_name = tolower(species_name)) %>%
    subset(species_name != species_out) %>%
    subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample)
  
  richness_data <- NULL
  for (i in unique(richness$site_number)){
    richness_data <- rbind(richness_data,
                           richness %>% subset(site_number == i) %>%
                             dplyr::group_by(site_number, study_block, study_number, study_sample) %>%
                             dplyr::summarise(richness = length(unique(species_name))))
  }
  
  df <- 
    data %>% subset(species_studied == "multiple") %>% 
    mutate(species_name = tolower(species_name)) %>%
    subset(species_name != species_out) %>%
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
  
  } else if(genus == "Aedes"){
    
    richness <- data %>% subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(species_name != species_out) %>%
      subset(genus == "Aedes") %>%
      subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample)
    
    richness_data <- NULL
    for (i in unique(richness$site_number)){
      richness_data <- rbind(richness_data,
                             richness %>% subset(site_number == i) %>%
                               dplyr::group_by(site_number, study_block, study_number, study_sample) %>%
                               dplyr::summarise(richness = length(unique(species_name))))
    }
    
    df <- 
      data %>% subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(species_name != species_out) %>%
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
    
  } else {
    richness <- data %>% subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(species_name != species_out) %>%
      subset(genus == "Anopheles") %>%
      subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample)
    
    richness_data <- NULL
    for (i in unique(richness$site_number)){
      richness_data <- rbind(richness_data,
                             richness %>% subset(site_number == i) %>%
                               dplyr::group_by(site_number, study_block, study_number, study_sample) %>%
                               dplyr::summarise(richness = length(unique(species_name))))
    
    df <- 
      data %>% subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(species_name != species_out) %>%
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
    }
    
  }
  
  return(df_all)
}
save(create_richness_df, file = "functions/create_richness_df.RData")

###########################################################################
# Create richness df for deforestation analysis
create_deforest_richness_df <- function(genus){
  
  data <- read.csv("data/inla_input/deforestation_data.csv") 
  
  # Total species richness
  if(genus == "total"){
    
    richness <- data %>% 
      subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample, sample_start_year, cells_deforestation, biome)
    
    richness_data <- NULL
    for (i in unique(richness$site_number)){
      richness_data <- rbind(richness_data,
                             richness %>% subset(site_number == i) %>%
                               dplyr::group_by(site_number, study_block, study_number, study_sample, biome) %>%
                               dplyr::summarise(richness = length(unique(species_name)),
                                                deforestation = unique(cells_deforestation)))
    }
    
  } else if(genus == "Aedes"){
    
    richness <- data %>% subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(genus == "Aedes") %>%
      subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample, sample_start_year, cells_deforestation, biome)
    
    richness_data <- NULL
    for (i in unique(richness$site_number)){
      richness_data <- rbind(richness_data,
                             richness %>% subset(site_number == i) %>%
                               dplyr::group_by(site_number, study_block, study_number, study_sample, biome) %>%
                               dplyr::summarise(richness = length(unique(species_name)),
                                                deforestation = unique(cells_deforestation)))
    }
    
  } else {
    # Anopheles richness
    richness <- data %>% subset(species_studied == "multiple") %>% 
      mutate(species_name = tolower(species_name)) %>%
      subset(genus == "Anopheles") %>%
      subset(measurement_adj > 0) %>% dplyr::select(site_number, species_name, study_block, study_number, study_sample, sample_start_year, cells_deforestation, biome)
    
    richness_data <- NULL
    for (i in unique(richness$site_number)){
      richness_data <- rbind(richness_data,
                             richness %>% subset(site_number == i) %>%
                               dplyr::group_by(site_number, study_block, study_number, study_sample, biome) %>%
                               dplyr::summarise(richness = length(unique(species_name)),
                                                deforestation = unique(cells_deforestation)))
    }
  }
  
  return(richness_data)
}

save(create_deforest_richness_df, file = "functions/create_deforest_richness_df")

###########################################################################
# Extract fixed effects from land use models
# mod - model
# land use var - land use or land use intensity (chr)

extract_estimates <- function(mod, land_use_var){
  
  estimates <- mod$summary.fixed %>% as.data.frame() %>%
    rownames_to_column(var = "land_use") %>% dplyr::select(land_use, mean, `0.025quant`, `0.975quant`) %>%
    dplyr::rename(lci = `0.025quant`,
                  uci = `0.975quant`) 
  
  if (grepl("intensity", land_use_var) == TRUE){
    
    estimates <- estimates %>% 
      dplyr::mutate(land_use = factor(gsub("LUI", "", land_use),
                                      levels = c("(Intercept)","primary vegetation-substantial", 
                                                 "secondary vegetation-minimal",
                                                 "secondary vegetation-substantial", 
                                                 "managed-minimal", 
                                                 "managed-substantial", 
                                                 "urban-minimal",
                                                 "urban-substantial"))) %>%
      dplyr::mutate(land_use = recode(land_use, "(Intercept)" = "primary vegetation-minimal")) %>%
      dplyr::mutate(LU = unlist(lapply(strsplit(as.character(land_use), "-"), "[[", 1)),
                    intensity = unlist(lapply(strsplit(as.character(land_use), "-"), "[[", 2))) %>%
      # reorder factor levels
      dplyr::mutate(LU = factor(LU, levels = c("primary vegetation", "secondary vegetation",
                                               "managed", "urban"))) %>%
      dplyr::mutate(intensity = factor(intensity, levels = c("minimal", 
                                                             "substantial"))) 
  }else{
    estimates <- estimates %>% dplyr::mutate(land_use = factor(gsub("land_use", "", land_use),
                                                               levels = c("(Intercept)",
                                                                          "secondary vegetation",
                                                                          "managed", 
                                                                          "urban"))) %>%
      dplyr::mutate(land_use = recode(land_use, "(Intercept)" = "primary vegetation"))
    
  }
  return(estimates)
}
save(extract_estimates, file = "functions/extract_estimates.RData")

###########################################################################
# Extract fixed effects from land use models - version 2 (for aggregated land use classes)
# mod - model
# land use var - land use or land use intensity (chr)

extract_estimates2 <- function(mod, land_use_var){
  
  estimates <- mod$summary.fixed %>% as.data.frame() %>%
    rownames_to_column(var = "land_use") %>% dplyr::select(land_use, mean, `0.025quant`, `0.975quant`) %>%
    dplyr::rename(lci = `0.025quant`,
                  uci = `0.975quant`) 
  
  if (grepl("intensity", land_use_var) == TRUE){
    
    estimates <- estimates %>% 
      dplyr::mutate(land_use = factor(gsub("LUI", "", land_use),
                                      levels = c("(Intercept)","primary vegetation-substantial", 
                                                 "secondary vegetation",
                                                 "managed", 
                                                 "urban"))) %>%
      dplyr::mutate(land_use = recode(land_use, "(Intercept)" = "primary vegetation-minimal")) %>%
      dplyr::mutate(LU = c("primary vegetation", "primary vegetation", "secondary vegetation", "managed", "urban"),
                    intensity = c("minimal", "substantial", "combined", "combined", "combined")) %>%
      # reorder factor levels
      dplyr::mutate(LU = factor(LU, levels = c("primary vegetation", "secondary vegetation",
                                               "managed", "urban"))) %>%
      dplyr::mutate(intensity = factor(intensity, levels = c("minimal", 
                                                             "substantial",
                                                             "combined"))) 
  }else{
    estimates <- estimates %>% dplyr::mutate(land_use = factor(gsub("land_use", "", land_use),
                                                               levels = c("(Intercept)",
                                                                          "secondary vegetation",
                                                                          "managed", 
                                                                          "urban"))) %>%
      dplyr::mutate(land_use = recode(land_use, "(Intercept)" = "primary vegetation"))
    
  }
  return(estimates)
}
save(extract_estimates2, file = "functions/extract_estimates2.RData")

###########################################################################
# Extract fixed effects from deforestation models
# mod - model
# deforest var - linear or non-linear deforestation (chr)

extract_deforest_estimates <- function(mod, deforest_var){

  
  if (grepl("non linear", deforest_var) == TRUE){

    estimates <- mod$summary.random$deforestation %>% as.data.frame() %>%
      dplyr::rename(lci = `0.025quant`,
                    uci = `0.975quant`) %>%
      dplyr::select(ID, mean, lci, uci)
    
  }else{
    estimates <- mod$summary.fixed %>% as.data.frame() %>%
      rownames_to_column(var = "deforestation") %>% dplyr::select(deforestation, mean, `0.025quant`, `0.975quant`) %>%
      dplyr::rename(lci = `0.025quant`,
                    uci = `0.975quant`) %>%
      subset(deforestation == "deforestation")
    
  }
  return(estimates)
}
save(extract_deforest_estimates, file = "functions/extract_deforest_estimates.RData")

###########################################################################
# Plot fixed effects from land use models
# df - dataframe of model estimates to be plotted
# land_use_var - land use or land use intensity
# model response - abundance or richness

plot_estimates <- function(df, land_use_var, response){
  
  if (grepl("abundance", response) == TRUE){
    y_title <- "Abundance difference (%)"
  } else {
    y_title <- "Richness difference (%)"
  }
  
  if("genus" %in% names(df) == TRUE){
    
    if (grepl("intensity", land_use_var) == TRUE){
      p <- df %>% 
        ggplot(aes(LU, mean, split=intensity, col=genus)) + 
        geom_hline(yintercept=0, lty=2) +
        geom_point(position=position_dodge(width=0.6), aes(pch=intensity), alpha = 0.6, size = 3) +
        geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
        theme_light() + 
        theme(legend.title = element_blank()) +
        xlab("") + 
        ylab(y_title) +
        scale_x_discrete(labels = wrap_format(18))
    } else {
      p <- df %>% 
        ggplot(aes(land_use, mean, col=genus)) + 
        geom_hline(yintercept=0, lty=2) +
        geom_point(position=position_dodge(width=0.6), alpha = 0.6, size = 3) +
        geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
        theme_light() + 
        theme(legend.title = element_blank(),
              legend.position = "bottom") +
        xlab("") + 
        ylab(y_title) +
        scale_x_discrete(labels = wrap_format(18))
    }
    
  }else{
    if (grepl("intensity", land_use_var) == TRUE){
      p <- df %>% 
        ggplot(aes(LU, mean, split=intensity, col=LU)) + 
        geom_hline(yintercept=0, lty=2) +
        geom_point(position=position_dodge(width=0.6), aes(pch=intensity), alpha = 0.6, size = 3) +
        geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
        theme_light() + 
        theme(legend.title = element_blank(),
              legend.position = "bottom") +
        xlab("") + 
        ylab(y_title) + 
        guides(color = "none") +
        scale_x_discrete(labels = wrap_format(18))
    } else {
      p <- df %>% 
        ggplot(aes(land_use, mean, col=land_use)) + 
        geom_hline(yintercept=0, lty=2) +
        geom_point(position=position_dodge(width=0.6), alpha = 0.6, size = 3) +
        geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
        theme_light() + 
        theme(legend.position = "none") +
        xlab("") + 
        ylab(y_title) +
        guides(color = "none") +
        scale_x_discrete(labels = wrap_format(18))
    }
    
  }
  return(p)
}
save(plot_estimates, file = "functions/plot_estimates.RData")

###########################################################################
# Plot estimates for cross-validation
# df - dataframe of model estimates to be plotted
# land_use_var - land use or land use intensity
# model response - abundance or richness

plot_estimates_cross <- function(df, land_use_var, response){
  
  if (grepl("abundance", response) == TRUE){
    y_title <- "Abundance difference (%)"
  } else {
    y_title <- "Richness difference (%)"
  }
  
    if (grepl("intensity", land_use_var) == TRUE){
      p <- df %>% 
        ggplot(aes(LU, mean, split=intensity, col=model)) + 
        geom_hline(yintercept=0, lty=2) +
        geom_point(position=position_dodge(width=0.6), aes(pch=intensity), alpha = 0.6, size = 3) +
        geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
        theme_light() + 
        theme(legend.title = element_blank(),
              legend.position = "bottom",
              plot.subtitle = element_text(face="bold")) +
        xlab("") + 
        ylab(y_title) +
        scale_x_discrete(labels = wrap_format(18))
    } else {
      p <- df %>% 
        ggplot(aes(land_use, mean, col=model)) + 
        geom_hline(yintercept=0, lty=2) +
        geom_point(position=position_dodge(width=0.6), alpha = 0.6, size = 3) +
        geom_errorbar(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.6), width=0.2, lwd=0.5, alpha=1, show.legend = FALSE) + 
        theme_light() + 
        theme(legend.position = "none",
              plot.subtitle = element_text(face="bold")) +
        xlab("") + 
        ylab(y_title) +
        scale_x_discrete(labels = wrap_format(18))
    }
    
  return(p)
}
save(plot_estimates_cross, file = "functions/plot_estimates_cross.RData")

###########################################################################

# Extract fitted values and residuals from model

# mod - model to be plotted
# genus - mosquito genus (Aedes or Anopheles, or total)
# model response - abundance or richness

extract_residuals <- function(mod, genus, response){
  
  if (response == "abundance"){
    
    if (grepl("Anopheles", genus)==TRUE){
      
      data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") %>% 
        subset(genus == "Anopheles") %>%
        mutate(abundance = measurement_adj) %>%
        dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
        dplyr::summarise(abundance = sum(abundance)) %>%
      dplyr::ungroup() %>%
      mutate(response_value = log(abundance)+1,
             response_fitted = c(mod$summary.fitted.values$mean)) %>%
      mutate(residuals = response_fitted - response_value,
             genus = genus)
  
    }else{if(grepl("Aedes", genus)==TRUE){
      data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") %>% 
        subset(genus == "Aedes") %>%
        mutate(abundance = measurement_adj) %>%
        dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
        dplyr::summarise(abundance = sum(abundance)) %>%
        dplyr::ungroup() %>%
        mutate(response_value = log(abundance)+1,
               response_fitted = c(mod$summary.fitted.values$mean)) %>%
        mutate(residuals = response_fitted - response_value,
               genus = genus)
      
      
      
    }  else{
      data <- read.csv("data/inla_input/abundance.csv") %>% aggregate_lui("all") %>% 
        mutate(abundance = measurement_adj) %>%
        dplyr::group_by(site_number, study_block, study_number, study_sample, species_name, LUI, biome) %>%
        dplyr::summarise(abundance = sum(abundance)) %>%
        dplyr::ungroup() %>%
        mutate(response_value = log(abundance)+1,
               response_fitted = c(mod$summary.fitted.values$mean)) %>%
        mutate(residuals = response_fitted - response_value,
               genus = genus)
    }
    }
  } else {
    if(grepl("Anopheles", genus)==TRUE){
      data <- read.csv("data/inla_input/ano_richness.csv") %>%
        mutate(response_value = richness,
               response_fitted = mod$summary.fitted.values$mean) %>%
        mutate(residuals = response_fitted - response_value,
               genus = "Anopheles")
    } else {
      if(grepl("Aedes", genus)==TRUE){
        data <- read.csv("data/inla_input/aed_richness.csv") %>%
          mutate(response_value = richness,
                 response_fitted = mod$summary.fitted.values$mean) %>%
          mutate(residuals = response_fitted - response_value,
                 genus = "Aedes")
      } else{
        data <- read.csv("data/inla_input/richness.csv") %>%
          mutate(response_value = richness,
                 response_fitted = mod$summary.fitted.values$mean) %>%
          mutate(residuals = response_fitted - response_value,
                 genus = "Total")
      }
    }
  }
  return(data)
  
}
save(extract_residuals, file = "functions/extract_residuals.RData")
