
library(blockrand)
library(dplyr)
library(arsenal)
library(writexl)


# Define the custom function
RandCapGen <- function(n, id_prefix, block_prefix, block_size,
                       arms, ratio = NULL, strat_vars = NULL,
                       strat_labels = NULL, seed, project_acronym) {
  # Set the initial seed for generating seeds for each stratum
  set.seed(seed)
  
  # Automatically determine the number of levels based on the length of 'arms'
  num.levels <- length(arms)
  
  # Generate balanced levels if a ratio is provided
  if (!is.null(ratio)) {
    # Function to generate levels according to the specified ratio
    generate_levels_with_ratio <- function(arms, ratio) {
      unlist(mapply(rep, arms, ratio))
    }
    # Use the generated levels with ratio
    balanced_levels <- generate_levels_with_ratio(arms, ratio)
  } else {
    # If no ratio is provided, use the original arms
    balanced_levels <- arms
  }
  
  # Calculate block sizes by dividing each element in block_size by the number of balanced levels
  block.sizes <- block_size / length(balanced_levels)
  
  # Check if stratification variables are provided
  if (is.null(strat_vars) || is.null(strat_labels)) {
    # No stratification: Perform a single global block randomization
    rand_df <- blockrand(
      n = n,
      num.levels = length(balanced_levels),
      id.prefix = id_prefix,
      block.prefix = block_prefix,
      block.sizes = block.sizes,
      levels = balanced_levels
    )
    
    # Rename the treatment column to "rando_bras"
    colnames(rand_df)[colnames(rand_df) == "treatment"] <- "rando_bras"
    
    # Full dataset with only the relevant columns
    final_df <- rand_df[, c("id", "block.id", "block.size", "rando_bras")]
    
    # Simplified dataset with only "rando_bras"
    simplified_df <- final_df[, "rando_bras", drop = FALSE]
    
  } else {
    # Stratified randomization with combinations of stratification variables
    
    # Generate all combinations of stratification variables
    strat_combinations <- expand.grid(strat_vars)
    
    # Generate a unique seed for each stratum
    stratum_seeds <- sample(1:10000, nrow(strat_combinations), replace = FALSE)
    
    # Initialize an empty list to store randomization dataframes
    all_strata <- list()
    
    # Define the function to shuffle the blocks
    melanger_blocs <- function(df) {
      # Split the dataframe into blocks based on block.id
      blocs <- split(df, df$block.id)
      # Shuffle the blocks
      blocs_melanges <- sample(blocs)
      # Recombine the shuffled blocks into a final dataframe
      df_melange <- do.call(rbind, blocs_melanges)
      return(df_melange)
    }
    
    # Iterate through each combination of stratification variables
    for (i in 1:nrow(strat_combinations)) {
      # Set the seed for this specific stratum
      set.seed(stratum_seeds[i])
      
      # Get the current combination of stratification levels as a data.frame
      strat_comb <- strat_combinations[i, , drop = FALSE]
      
      # Retrieve the stratification variable name(s)
      strat_var_names <- names(strat_comb)
      
      # Create a custom label for the stratum based on the combination of
      # stratification variables and their labels
      if (ncol(strat_combinations) == 1) {
        stratum_label <- paste0(strat_labels[[strat_var_names]], "_", strat_comb[[1]])
      } else {
        stratum_label <- paste(
          mapply(function(var_name, value) {
            paste0(strat_labels[[var_name]], "_", value)
          }, strat_var_names, strat_comb),
          collapse = "|"
        )
      }
      
      # Define a custom block prefix with stratum order (S1, S2, ...) if stratification is present
      custom_block_prefix <- paste0("S", i, "_", block_prefix)
      custom_id_prefix <- paste0("S", i, "_", id_prefix)
      # Perform block randomization for this specific stratum
      rand_df <- blockrand(
        n = n,
        num.levels = length(balanced_levels),
        id.prefix = custom_id_prefix,
        block.prefix = custom_block_prefix,
        stratum = stratum_label,
        block.sizes = block.sizes,
        levels = balanced_levels
      )
      
      # Rename the treatment column to "rando_bras"
      colnames(rand_df)[colnames(rand_df) == "treatment"] <- "rando_bras"
      
      # Add stratification variables to the dataframe
      for (strat_var_name in strat_var_names) {
        if (tolower(strat_var_name) %in% c("center", "centre")) {
          new_var_name <- "redcap_data_access_group"
        } else {
          new_var_name <- paste0("rando_", strat_var_name)
        }
        
        rand_df[[new_var_name]] <- strat_comb[[strat_var_name]]
      }
      
      # Shuffle the blocks within the randomized dataframe
      rand_df_shuffled <- melanger_blocs(rand_df)
      
      # Append the shuffled randomization dataframe for this stratum to the list
      all_strata[[stratum_label]] <- rand_df_shuffled
    }
    
    # Combine all strata into one final dataframe
    final_df <- bind_rows(all_strata)
    
    # Create the simplified dataset with only 'rando_bras' and stratification variables
    stratification_columns <- grep("rando_|redcap_data_access_group", colnames(final_df), value = TRUE)
    simplified_df <- final_df[, stratification_columns, drop = FALSE]
  }
  
  # Convert all columns in simplified_df to factors
  simplified_df <- simplified_df %>% mutate(across(everything(), as.factor))
  
  rownames(final_df) <- NULL
  rownames(simplified_df) <- NULL
  
  # Create the output structure
  randomization_object <- list(
    settings = list(
      R_version = R.version.string,
      date = Sys.time(),
      initial_seed = seed,
      treatment_arms = arms,
      ratios = ratio,
      strat_vars = if (!is.null(strat_vars)) names(strat_vars) else NULL,
      strat_labels = strat_labels,
      strat_levels = strat_vars,  # Includes levels of stratification variables
      block_size_input = block_size,
      stratum_seeds = if (exists("stratum_seeds")) stratum_seeds else NULL, # Seeds for strata, if any
      project_acronym=project_acronym
    ),
    tables = list(
      full_dataset = final_df,
      simplified_dataset = simplified_df
    )
  )
  
  # Return the structured object
  return(randomization_object)
}


# Example usage:
stratification_vars <- list(terme_strate = c(1, 2), 
                            Centre = c(832,
                                       884,
                                       885,
                                       886,
                                       887,
                                       888,
                                       889,
                                       890),
                            Sex=c(1,2),
                            Poids=c(1,2,3),
                            Diabete=c(1,0),
                            HTA=c(1,0))
stratification_labels <- list(terme_strate = "term", Centre = "c",
                              Sex="sex",
                              Poids="pds",
                              Diabete="diab",
                              HTA="HTA")

# 
stratification_vars <- list(terme_strate = c(1, 2))
stratification_labels <- list(terme_strate = "term")


stratification_vars <- list(Centre = c(832, 884,885,886, 887, 888,889, 890))
stratification_labels <- list(Centre = "c")

# VARIBLE STRATA
results <- RandCapGen(
  n = 3600,
  id_prefix  = 'P',
  block_prefix = 'B',
  block_size = c(24,48),  # New argument bloc_size
  arms = c("A", "B","C","D","E","F",'G',"H","I"),
  ratio = c(1,3,2,1,1,2,4,9,1),
  strat_vars = stratification_vars,
  strat_labels = stratification_labels,
  seed = 22,
  project_acronym = "ANTEIPA"
)


#SANS STRATA
results <- RandCapGen(
  n = 144,
  id_prefix = 'P_',
  block_prefix = 'B_',
  block_size = c(14,28),  # New argument bloc_size
  arms = c("A", "B","C","D","E","F",'G'),
  ratio = c(1,3,2,1,1,2,4),
  strat_vars = NULL,
  strat_labels = NULL,
  seed = 22,
  project_acronym = "VOLEM"
)


results <- RandCapGen(
  n = 144,
  id_prefix = 'P_',
  block_prefix = 'B_',
  block_size = c(4,8),  # New argument bloc_size
  arms = c("A", "B"),
  ratio = c(1,3),
  seed = 22,
  project_acronym = "VOLEM"
)

# Display the results
full_dataset <- results$tables$full_dataset
simplified_dataset <- results$tables$full_dataset





library(grid)
library(gridExtra)
library(reshape2)

RandCapBalance <- function(randomization_object, output_path = "Randomization_Balance.pdf") {
  # Extraction du jeu de données complet et calcul du nombre de strates
  full_dataset <- randomization_object$tables$full_dataset
  num_strata <- length(unique(full_dataset$stratum))
  
  
  # Calcul de l’équilibre par bras de traitement au sein de chaque stratum
  if(num_strata==0){
    full_dataset$stratum<-"No stratum"
    balance_table <- table(Stratum = full_dataset$stratum, Treatment = full_dataset$rando_bras )
  }else{
    balance_table <- table(Stratum = full_dataset$stratum, Treatment = full_dataset$rando_bras )
  }
  
  
  # Convertir en data frame pour afficher les strates dans une colonne distincte
  balance_df <- as.data.frame(balance_table)
  balance_df <- reshape2::dcast(balance_df, Stratum ~ Treatment, value.var = "Freq")
  
  # Informations de configuration de randomisation
  strat_vars <- randomization_object$settings$strat_vars
  strat_vars_text <- if (is.null(strat_vars) || length(strat_vars) == 0) {
    "There is no stratification variable."
  } else {
    paste(strat_vars, collapse = ", ")
  }
  
  ratios <- paste(randomization_object$settings$ratio, collapse = ", ")
  treatment_arms <- paste(randomization_object$settings$treatment_arms, collapse = ", ")
  
  # Sauvegarde en PDF avec une mise en page améliorée
  pdf(output_path, width = 8.5, height = 11)
  
  # Titre et informations résumées
  grid.newpage()
  y_pos <- 0.9
  grid.text("Randomization Balance Summary", x = 0.5, y = y_pos, just = "center", 
            gp = gpar(fontsize = 18, fontface = "bold", col = "#2E86C1"))
  y_pos <- y_pos - 0.05
  
  # Texte pour le nombre de strates et la présence de stratification
  if (num_strata > 1) {
    grid.text(paste("Number of Strata:", num_strata), x = 0.1, 
              y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", 
                        col = "#2874A6"))
  } else {
    grid.text("Number of Strata:", x = 0.1, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", 
                        col = "#2874A6"))
    grid.text("The randomization is not stratified.", x = 0.3, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))  # Texte normal et noir
  }
  y_pos <- y_pos - 0.03
  
  # Texte pour les variables de stratification
  grid.text(paste("Stratification Variables:", strat_vars_text), x = 0.1, y = y_pos, just = "left", 
            gp = gpar(fontsize = 12, col = "#34495E"))
  y_pos <- y_pos - 0.03
  
  grid.text(paste("Treatment Arms:", treatment_arms), x = 0.1, y = y_pos, just = "left", 
            gp = gpar(fontsize = 12, col = "#34495E"))
  y_pos <- y_pos - 0.03
  grid.text(paste("Ratios:", ratios), x = 0.1, y = y_pos, just = "left", 
            gp = gpar(fontsize = 12, col = "#34495E"))
  
  # Ajustement pour la position de la table sur la première page
  first_page_y_pos <- y_pos - 0.05
  rows_per_page <- 22
  start_row <- 1
  end_row <- rows_per_page
  num_rows <- nrow(balance_df)
  
  # Thème pour la table avec des couleurs pour les colonnes
  table_theme <- ttheme_minimal(
    core = list(
      fg_params = list(col = "black"),
      bg_params = list(
        # Application des couleurs : première colonne et colonnes de traitement
        fill = c(
          rep("#D1E8FF", nrow(balance_df)),  # Couleur pour la colonne des strates
          rep(c("#EBF5FB", "#EBF5FB"), each = nrow(balance_df) * (ncol(balance_df) - 1) / 2) # Couleur pour les colonnes de traitement
        ),
        col = NA
      )
    ),
    colhead = list(
      fg_params = list(col = "white", fontface = "bold"),
      bg_params = list(fill = "#2E86C1")  # Couleur pour la ligne de titre
    )
  )
  
  # Pagination et position dynamique
  first_page <- TRUE
  while (start_row <= num_rows) {
    # Nouvelle page si continuation
    if (!first_page) {
      grid.newpage()
      y_pos <- 0.9  # Positionner au début de la page pour les pages suivantes
    } else {
      y_pos <- first_page_y_pos  # Position de départ après "Ratios" pour la première page
    }
    
    # Sélection des lignes pour la page courante
    table_page <- balance_df[start_row:min(end_row, num_rows), ]
    
    # Création de la table avec l'en-tête répétée
    table_grob <- tableGrob(table_page, rows = NULL, theme = table_theme)
    
    # Dessiner la table en commençant à `y_pos`
    pushViewport(viewport(y = y_pos, height = unit(0.7, "npc"), just = "top"))
    grid.draw(table_grob)
    popViewport()
    
    # Mise à jour des indices de ligne pour la page suivante
    start_row <- end_row + 1
    end_row <- start_row + rows_per_page - 1
    first_page <- FALSE  # Indiquer que les pages suivantes sont des continuations
  }
  
  # Fermer le PDF
  dev.off()
  
  message("Randomization balance summary saved to ", output_path)
  return(balance_df)
}
RandCapBalance(results)




RandCapSettings <- function(rand_obj, output_path = "Randomization_settings.pdf", 
                            display_block_size = FALSE, display_generated_seeds = FALSE) {
  # Démarrage du fichier PDF
  pdf(output_path, width = 8.5, height = 11, family = "serif")
  
  # Titre
  grid.newpage()
  grid.text("Randomization Settings", x = 0.05, y = 0.9, just = "left", 
            gp = gpar(fontsize = 20, fontface = "bold", col = "#2E86C1"))
  
  # Variables
  r_version <- rand_obj$settings$R_version
  date_randomization <- rand_obj$settings$date
  seed <- rand_obj$settings$initial_seed
  treatment_arms <- paste(rand_obj$settings$treatment_arms, collapse = ", ")
  ratios <- paste(rand_obj$settings$ratio, collapse = ", ")
  strat_vars <- rand_obj$settings$strat_vars
  strat_levels <- rand_obj$settings$strat_levels
  block_size <- paste(rand_obj$settings$block_size_input, collapse = ", ")
  generated_seeds <- rand_obj$settings$stratum_seeds
  
  # Fonction d'affichage avec ajustement de l'alignement
  display_setting <- function(label, value, y_pos) {
    # Label en gras et en bleu
    grid.text(label, x = 0.05, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    # Valeur en noir avec espacement ajusté
    grid.text(value, x = 0.4, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, col = "#34495E"))
  }
  
  # Fonction d'affichage avec retour à la ligne des seeds
  display_wrapped_seeds <- function(label, seeds, y_pos, wrap_length = 5) {
    # Label pour les generated seeds
    grid.text(label, x = 0.05, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    # Affichage des seeds avec retour à la ligne
    wrapped_seeds <- paste(seeds, collapse = ", ")
    seed_lines <- strwrap(wrapped_seeds, width = wrap_length)
    for (line in seed_lines) {
      y_pos <- y_pos - 0.03
      grid.text(line, x = 0.4, y = y_pos, just = "left", 
                gp = gpar(fontsize = 12, col = "#34495E"))
    }
    return(y_pos)
  }
  
  # Affichage des paramètres avec ajustements d'alignement
  y_pos <- 0.85
  display_setting("R Version:", r_version, y_pos); y_pos <- y_pos - 0.03
  display_setting("Date of Randomization:", date_randomization, y_pos); y_pos <- y_pos - 0.03
  display_setting("Initial Seed:", seed, y_pos); y_pos <- y_pos - 0.03
  
  # Affichage conditionnel des seeds générés
  if (display_generated_seeds) {
    y_pos <- display_wrapped_seeds("Generated Seeds:", generated_seeds, y_pos, wrap_length = 40)
    y_pos <- y_pos - 0.03
  }
  
  display_setting("Treatment Arms:", treatment_arms, y_pos); y_pos <- y_pos - 0.03
  display_setting("Ratios:", ratios, y_pos); y_pos <- y_pos - 0.03
  
  # Affichage des variables de stratification et niveaux avec retour à la ligne rapproché
  if (is.null(strat_vars) || length(strat_vars) == 0) {
    # Pas de variable de stratification
    grid.text("Stratification:", x = 0.05, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid.text("There is no stratification variable.", x = 0.4, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    # Afficher les variables de stratification
    grid.text("Stratification Variables:", x = 0.05, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    y_pos <- y_pos - 0.03
    for (i in seq_along(strat_vars)) {
      grid.text(paste0("− ", strat_vars[i], ": ", paste(strat_levels[[i]], collapse = ", ")),
                x = 0.1, y = y_pos, just = "left",
                gp = gpar(fontsize = 12, col = "#34495E"))
      y_pos <- y_pos - 0.025
    }
  }
  
  # Affichage conditionnel de la taille des blocs
  y_pos <- y_pos - 0.03
  grid.text("Block Size per Stratum:", x = 0.05, y = y_pos, just = "left", 
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  
  if (display_block_size) {
    grid.text(block_size, x = 0.4, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    grid.text("Not allowed to be displayed", x = 0.4, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, col = "#E74C3C"))
  }
  
  # Fermeture du PDF
  dev.off()
  
  message("Randomization settings saved to ", output_path)
}

RandCapSettings(results, display_block_size = F)


RandCapTable <- function(randomization_object, save_for_REDCap = TRUE, save_random_table = TRUE) {
  # Sauvegarde la table simplifiée si save_for_REDCap est TRUE
  if (save_for_REDCap) {
    redcap_table_name <- paste0(randomization_object$settings$project_acronym, "_REDCap_table","_",Sys.Date(),".csv")
    write.csv(randomization_object$tables$simplified_dataset, redcap_table_name, row.names = FALSE)
  }
  
  # Sauvegarde la table complète si save_random_table est TRUE
  if (save_random_table) {
    randomization_table_name <- paste0(randomization_object$settings$project_acronym, "_randomization_table","_",Sys.Date(),".csv")
    write.csv(randomization_object$tables$full_dataset, randomization_table_name, row.names = FALSE)
  }
  
  # Retourne les deux tables dans une liste
  return(list(
    simplified_table = randomization_object$tables$simplified_dataset,
    full_table = randomization_object$tables$full_dataset
  ))
}




RandCapTable(results,save_for_REDCap = T, save_random_table = T)

