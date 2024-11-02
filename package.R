
library(blockrand)
library(grid)
library(gridExtra)
library(tidyverse)

## RANDCAPGEN---------
# Define the custom function
RandCapGen <- function(n, id_prefix=NULL, block_prefix=NULL, block_sizes,
                       arms, ratio = NULL, strat_vars = NULL,
                       strat_vars_prefix = NULL, seed=NULL,
                       uneq_beg = FALSE, uneq_mid = FALSE, uneq_min = 1, 
                       uneq_maxit = 1000,project_acronym=NULL) {
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
  
  # Calculate block sizes by dividing each element in block_sizes by 
  #the number of balanced levels
  block.sizes <- block_sizes / length(balanced_levels)
  
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
  
  # Check if stratification variables are provided
  if (is.null(strat_vars)) {
    # No stratification: Perform a single global block randomization
    rand_df <- blockrand(
      n = n,
      num.levels = length(balanced_levels),
      id.prefix = id_prefix,
      block.prefix = block_prefix,
      block.sizes = block.sizes,
      levels = balanced_levels,
      uneq.beg = uneq_beg,
      uneq.mid = uneq_mid,
      uneq.min = uneq_min,
      uneq.maxit = uneq_maxit
    )
    
    # Rename the treatment column to "rando_bras"
    colnames(rand_df)[colnames(rand_df) == "treatment"] <- "rando_bras"
    
    rand_df <- melanger_blocs(rand_df)
    
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
        stratum_label <- paste0(strat_vars_prefix[[strat_var_names]], "",
                                strat_comb[[1]])
      } else {
        stratum_label <- paste(
          mapply(function(var_name, value) {
            paste0(strat_vars_prefix[[var_name]], "", value)
          }, strat_var_names, strat_comb),
          collapse = "|"
        )
      }
      
      # Define a custom block prefix with stratum order 
      #(S1, S2, ...) if stratification is present
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
        levels = balanced_levels,
        uneq.beg = uneq_beg,
        uneq.mid = uneq_mid,
        uneq.min = uneq_min,
        uneq.maxit = uneq_maxit
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
    stratification_columns <- grep("rando_|redcap_data_access_group", 
                                   colnames(final_df), value = TRUE)
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
      strat_vars_prefix = strat_vars_prefix,
      strat_levels = strat_vars,  # Includes levels of stratification variables
      block_sizes_input = block_sizes,
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


#RANDCAP BALANCE--------
RandCapBalance <- function(randomization_object, 
                           output_path = "Randomization_Balance.pdf") {
  # Extract the full dataset and calculate the number of strata
  full_dataset <- randomization_object$tables$full_dataset
  num_strata <- length(unique(full_dataset$stratum))
  
  # Handle the case where no strata are defined
  if (num_strata == 0) {
    full_dataset <- full_dataset %>% mutate(stratum = "No stratum")
  }
  
  # Create the balance table with dplyr and tidyr
  balance_df <- full_dataset %>%
    count(stratum, rando_bras) %>%
    pivot_wider(names_from = rando_bras, values_from = n, values_fill = 0) %>%
    mutate(Total = rowSums(select(., -stratum)))
  
  # Get randomization settings
  strat_vars <- randomization_object$settings$strat_vars
  strat_vars_text <- if (is.null(strat_vars) || length(strat_vars) == 0) {
    "There is no stratification variable."
  } else {
    paste(strat_vars, collapse = ", ")
  }
  
  ratios <- if (is.null(randomization_object$settings$ratio)) {
    "Equal size ratio"
  } else {
    paste(randomization_object$settings$ratio, collapse = ", ")
  }
  
  treatment_arms <- paste(randomization_object$settings$treatment_arms, collapse = ", ")
  
  # Start the PDF
  pdf(output_path, width = 8.5, height = 11)
  
  # Determine the layout for the first page based on the number of rows in the table
  grid.newpage()
  
  if (nrow(balance_df) >= 20) {
    # Case with at least 20 rows in the table
    text_layout_height <- 0.3
    table_layout_height <- 0.7
    table_rows <- 20
  } else {
    # Case with fewer than 20 rows in the table
    text_layout_height <- 0.3
    table_layout_height <- min(1, 0.025 * nrow(balance_df) + 0.1) # 2.5% * nrow + offset for spacing
    remaining_height <- 1 - (text_layout_height + table_layout_height)
    table_rows <- nrow(balance_df)
  }
  
  # Display text in the upper section (30% of the page)
  pushViewport(viewport(layout = grid.layout(10, 1)))  # Divide into 10 for proportional layout
  pushViewport(viewport(layout.pos.row = 1:3, layout.pos.col = 1))  # Upper section for text
  
  # Display headers and values with blue and bold style for headers and black for values
  grid.text("Randomization Balance Summary", x = 0.25, y = 0.8, just = "left", 
            gp = gpar(fontsize = 18, fontface = "bold", col = "#2E86C1"))
  
  # Display number of strata
  if (num_strata == 0) {
    grid.text("Number of Strata:", x = 0.15, y = 0.55, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid.text("The randomization is not stratified", x = 0.4, y = 0.55, just = "left", 
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    grid.text("Number of Strata:", x = 0.15, y = 0.55, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid.text(num_strata, x = 0.4, y = 0.55, just = "left", 
              gp = gpar(fontsize = 12, col = "#34495E"))
  }
  
  grid.text("Stratification variables:", x = 0.15, y = 0.4, just = "left", 
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  grid.text(strat_vars_text, x = 0.4, y = 0.4, just = "left", 
            gp = gpar(fontsize = 12, col = "#34495E"))
  
  grid.text("Treatment Arms:", x = 0.15, y = 0.25, just = "left", 
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  grid.text(treatment_arms, x = 0.4, y = 0.25, just = "left", 
            gp = gpar(fontsize = 12, col = "#34495E"))
  
  grid.text("Ratios:", x = 0.15, y = 0.1, just = "left", 
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  grid.text(ratios, x = 0.4, y = 0.1, just = "left", 
            gp = gpar(fontsize = 12, col = "#34495E"))
  
  popViewport()
  
  # Display the table in the lower section
  table_grob1 <- tableGrob(balance_df[1:table_rows, ], rows = NULL)
  table_theme_page1 <- ttheme_minimal(
    core = list(
      fg_params = list(col = "black"),
      bg_params = list(
        fill = c(
          rep("#D1E8FF", table_rows),                             # Color for strata column
          rep("#EBF5FB", table_rows * (ncol(balance_df) - 2)),    # Color for treatment columns
          rep("#D1E8FF", table_rows)                              # Color for Total column
        ),
        col = NA
      )
    ),
    colhead = list(
      fg_params = list(col = "white", fontface = "bold"),
      bg_params = list(fill = "#2E86C1")
    )
  )
  
  pushViewport(viewport(y = 1 - text_layout_height, 
                        height = unit(table_layout_height, "npc"), 
                        just = "top"))
  grid.draw(tableGrob(balance_df[1:table_rows, ], rows = NULL, 
                      theme = table_theme_page1))
  popViewport()
  
  # Add an empty layout below if there are fewer than 20 rows in the table
  if (nrow(balance_df) < 20) {
    pushViewport(viewport(y = text_layout_height, height = 
                            unit(remaining_height, "npc"), just = "top"))
    popViewport()
  }
  
  # Pages for remaining rows
  start_row <- table_rows + 1
  num_rows <- nrow(balance_df)
  rows_per_page <- 31
  
  while (start_row <= num_rows) {
    grid.newpage()
    
    # Check if this is the last page
    is_last_page <- (start_row + rows_per_page - 1 >= num_rows)
    end_row <- if (is_last_page) num_rows else (start_row + rows_per_page - 1)
    table_page <- balance_df[start_row:end_row, ]
    
    # Colors for each page
    page_rows <- nrow(table_page)
    page_fill_colors <- c(
      rep("#D1E8FF", page_rows),                       # Color for strata column
      rep("#EBF5FB", page_rows * (ncol(balance_df) - 2)), # Color for treatment columns
      rep("#D1E8FF", page_rows)                        # Color for Total column
    )
    
    table_theme <- ttheme_minimal(
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = page_fill_colors, col = NA)
      ),
      colhead = list(
        fg_params = list(col = "white", fontface = "bold"),
        bg_params = list(fill = "#2E86C1")
      )
    )
    table_grob <- tableGrob(table_page, rows = NULL, theme = table_theme)
    
    if (is_last_page) {
      # Last page: Adjust the height based on remaining rows with a factor of 3%
      remaining_rows <- end_row - start_row + 1
      last_page_height <- min(0.9, 0.025 * remaining_rows)
      
      # Position table in the upper section of the last page
      pushViewport(viewport(y = 0.9, height = unit(last_page_height, "npc"), just = "top"))
    } else {
      # Intermediate pages with table at the top
      pushViewport(viewport(y = 0.87, height = unit(0.75, "npc"), just = "top"))
    }
    
    grid.draw(table_grob)
    popViewport()
    
    # Move to the next set of rows
    start_row <- end_row + 1
  }
  
  # Close the PDF
  dev.off()
  
  message("Randomization balance summary saved to ", output_path)
  return(balance_df)
}


# RANDCAP SETTINGS FUNCTION-------


RandCapSettings <- function(rand_obj, output_path = "Randomization_settings.pdf", 
                            display_block_sizes = FALSE) {
  # Start PDF file
  pdf(output_path, width = 8.5, height = 11, family = "serif")
  
  # Title
  grid.newpage()
  grid.text("Randomization Settings", x = 0.1, y = 0.9, just = "left", 
            gp = gpar(fontsize = 20, fontface = "bold", col = "#2E86C1"))
  
  # Variables
  r_version <- rand_obj$settings$R_version
  date_randomization <- rand_obj$settings$date
  seed <- if (is.null(rand_obj$settings$initial_seed)
  ) "No seed was specified" else rand_obj$settings$initial_seed
  treatment_arms <- paste(rand_obj$settings$treatment_arms, collapse = ", ")
  ratios <- if (is.null(rand_obj$settings$ratio)
  ) "Equal size ratio" else paste(rand_obj$settings$ratio, 
                                  collapse = ", ")
  strat_vars <- rand_obj$settings$strat_vars
  strat_levels <- rand_obj$settings$strat_levels
  block_sizes <- paste(rand_obj$settings$block_sizes_input, collapse = ", ")
  
  # Display function with alignment adjustments
  display_setting <- function(label, value, y_pos) {
    # Bold blue label
    grid.text(label, x = 0.1, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    # Black value with adjusted spacing
    grid.text(value, x = 0.45, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, col = "#34495E"))
  }
  
  # Display settings with alignment adjustments
  y_pos <- 0.85
  display_setting("R Version:", r_version, y_pos); y_pos <- y_pos - 0.03
  display_setting("Date of Randomization:", date_randomization, 
                  y_pos); y_pos <- y_pos - 0.03
  display_setting("Seed:", seed, y_pos); y_pos <- y_pos - 0.03
  display_setting("Treatment Arms:", treatment_arms, y_pos);
  y_pos <- y_pos - 0.03
  display_setting("Ratios:", ratios, y_pos); y_pos <- y_pos - 0.03
  
  # Stratification variables and levels with line wrapping for Centers
  if (is.null(strat_vars) || length(strat_vars) == 0) {
    # No stratification variable
    grid.text("Stratification:", x = 0.1, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid.text("There is no stratification variable.", x = 0.45,
              y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    # Display stratification variables
    grid.text("Stratification Variables:", x = 0.1, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    y_pos <- y_pos - 0.03
    
    for (i in seq_along(strat_vars)) {
      # Check for Center-related variables
      if (tolower(strat_vars[i]) %in% c("center", "centre", 
                                        "redcap_data_access_group")) {
        # Wrap levels for "Centers"
        levels_text <- paste(strat_levels[[i]], collapse = ", ")
        wrapped_levels <- strwrap(levels_text, width = 80)
        
        # Print "Centers:" as the label and wrap the levels across lines
        grid.text("- Centers:", x = 0.15, y = y_pos, just = "left",
                  gp = gpar(fontsize = 12))
        y_pos <- y_pos - 0.03
        
        for (line in wrapped_levels) {
          grid.text(line, x = 0.2, y = y_pos, just = "left",
                    gp = gpar(fontsize = 12, col = "#34495E"))
          y_pos <- y_pos - 0.03
        }
        
      } else {
        # For other stratification variables, display without line wrapping
        grid.text(paste0("− ", strat_vars[i], ": ", paste(strat_levels[[i]], 
                                                          collapse = ", ")),
                  x = 0.15, y = y_pos, just = "left",
                  gp = gpar(fontsize = 12, col = "#34495E"))
        y_pos <- y_pos - 0.025
      }
    }
  }
  
  # Conditional display of block sizes
  y_pos <- y_pos - 0.03
  grid.text("Block Size per Stratum:", x = 0.1, y = y_pos, just = "left", 
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  
  if (display_block_sizes) {
    grid.text(block_sizes, x = 0.45, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    grid.text("Not allowed to be displayed", x = 0.45, y = y_pos, just = "left", 
              gp = gpar(fontsize = 12, col = "#E74C3C"))
  }
  
  # Close the PDF
  dev.off()
  
  message("Randomization settings saved to ", output_path)
}



#RANDCAP TABLE--------
RandCapTable <- function(randomization_object, save_for_REDCap = TRUE, 
                         save_random_table = TRUE) {
  # Save the simplified table if save_for_REDCap is TRUE
  if (save_for_REDCap) {
    redcap_table_name <- paste0(randomization_object$settings$project_acronym,
                                "_REDCap_table", "_", Sys.Date(), ".csv")
    write.csv(randomization_object$tables$simplified_dataset,
              redcap_table_name, row.names = FALSE)
  }
  
  # Save the full table if save_random_table is TRUE
  if (save_random_table) {
    randomization_table_name <- paste0(
      randomization_object$settings$project_acronym, 
      "_randomization_table", "_", Sys.Date(), ".csv")
    write.csv(randomization_object$tables$full_dataset, 
              randomization_table_name, row.names = FALSE)
  }
  
  # Return both tables in a list
  return(list(
    simplified_table = randomization_object$tables$simplified_dataset,
    full_table = randomization_object$tables$full_dataset
  ))
}


