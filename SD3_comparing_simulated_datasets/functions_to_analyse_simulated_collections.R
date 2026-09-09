# Code to perform simulated collection analysis (SD3 analysis)

# Set seed for consistency
set.seed(678999821)


# NORMAL RANGES --------------------------------------------------- #
# ellipsoid_x_radii_range <- c(min = 75, max = 125)
# ellipsoid_y_radii_range <- c(min = 75, max = 125)
# ellipsoid_z_radii_range <- c(min = 75, max = 125)
# network_width_range <- c(min = 25, max = 35)
# 
# mixed_cell_type_A_proportion_range <- c(min = 0.3, max = 0.7)
# ringed_ring_width_factor_range <- c(min = 0.1, max = 0.2)
# separated_cluster1_x_coordinate_range <- c(min = 125, max = 175)


# OUTSIDE OF NORMAL RANGES ---------------------------------------- #
# ellipsoid_x_radii_range_small <- c(min = 50, max = 100)
# ellipsoid_x_radii_range_large <- c(min = 100, max = 150)
# network_width_range_small <- c(min = 20, max = 30)
# network_width_range_large <- c(min = 30, max = 40)
#
# mixed_cell_type_A_proportion_range_small <- c(min = 0.1, max = 0.5)
# mixed_cell_type_A_proportion_range_large <- c(min = 0.5, max = 0.9)
# ringed_ring_width_factor_range_small <- c(min = 0.05, max = 0.15)
# ringed_ring_width_factor_range_large <- c(min = 0.15, max = 0.25)
# separated_cluster1_x_coordinate_range_small <- c(min = 100, max = 150)
# separated_cluster1_x_coordinate_range_large <- c(min = 150, max = 200)


# OTHER PARAMETERS ------------------------------------------------ #

# BACKGROUND ***
# number_of_cells = 30000
# length = 600
# width = 600
# height = 300
# minimum_distance_between_cells = 10\

# SEPARATED SPHERE CLUSTER ***
# radius = 100
# centre_loc <- c(450, 300, 150)

# NETWORK CLUSTER ***
# edges = 20
# containment sphere radius = 125


# FUNCTIONS ------------------------------------------------------- #
# Output is a 'parameters' data frame with 'n_simulations' rows, and columns 
# containing parameters specific to the chosen shape and arrangement.
# Values are randomly chosen from within the 'normal ranges' specified above.
create_random_parameters_for_simulated_group <- function(n_simulations, 
                                                         shape, 
                                                         arrangement) {
  
  # Create data frame for parameters
  parameters_df <- data.frame(number = seq_len(n_simulations))
  
  # Set random parameter ranges here
  ellipsoid_x_radii_range <- c(min = 75, max = 125)
  ellipsoid_y_radii_range <- c(min = 75, max = 125)
  ellipsoid_z_radii_range <- c(min = 75, max = 125)
  network_width_range <- c(min = 25, max = 35)
  
  mixed_cell_type_A_proportion_range <- c(min = 0.3, max = 0.7)
  ringed_ring_width_factor_range <- c(min = 0.1, max = 0.2)
  separated_cluster1_x_coordinate_range <- c(min = 125, max = 175)
  
  
  ### Shapes
  if (shape == "ellipsoid") {
    ellipsoid_x_radii <- runif(n_simulations, 
                               min = ellipsoid_x_radii_range["min"], 
                               max = ellipsoid_x_radii_range["max"])
    
    ellipsoid_y_radii <- runif(n_simulations, 
                               min = ellipsoid_y_radii_range["min"], 
                               max = ellipsoid_y_radii_range["max"])
    
    ellipsoid_z_radii <- runif(n_simulations, 
                               min = ellipsoid_z_radii_range["min"], 
                               max = ellipsoid_z_radii_range["max"])
    
    parameters_df$ellipsoid_x_radius <- ellipsoid_x_radii
    parameters_df$ellipsoid_y_radius <- ellipsoid_y_radii
    parameters_df$ellipsoid_z_radius <- ellipsoid_z_radii
    
  }
  else if (shape == "network") {
    network_widths <- runif(n_simulations,
                            min = network_width_range["min"],
                            max = network_width_range["max"])
    
    parameters_df$network_width <- network_widths
  }
  else {
    stop(paste(shape, "is not a valid shape", 
               "Should be 'ellipsoid', 'network'."))
  }
  
  
  
  ### Arrangements
  if (arrangement == "mixed") {
    mixed_cell_type_A_proportions <- runif(n_simulations, 
                                           min = mixed_cell_type_A_proportion_range["min"], 
                                           max = mixed_cell_type_A_proportion_range["max"])
    
    parameters_df$mixed_cell_type_A_proportion <- mixed_cell_type_A_proportions
    
  }
  else if (arrangement == "ringed") {
    ringed_ring_width_factors <- runif(n_simulations,
                                       min = ringed_ring_width_factor_range["min"],
                                       max = ringed_ring_width_factor_range["max"])
    
    # Get ring width as well, depending on if it is an ellipsoid or network shape
    if (shape == "ellipsoid") {
      ringed_ring_widths <- ringed_ring_width_factors * 
        (ellipsoid_x_radii + ellipsoid_y_radii + ellipsoid_z_radii) / 3
    }
    else if (shape == "network") {
      ringed_ring_widths <- ringed_ring_width_factors * network_widths
    }
    parameters_df$ringed_ring_width_factor <- ringed_ring_width_factors
    parameters_df$ringed_ring_width <- ringed_ring_widths
    
  }
  else if (arrangement == "separated") {
    separated_cluster1_x_coordinates <- runif(n_simulations,
                                              min = separated_cluster1_x_coordinate_range["min"],
                                              max = separated_cluster1_x_coordinate_range["max"])
    parameters_df$separated_cluster1_x_coordinate <- separated_cluster1_x_coordinates
  }
  else {
    stop(paste(arrangement, "is not a valid arrangement.", 
               "Should be 'mixed', 'ringed', 'separated'."))
  }  
  
  return(parameters_df)
}

# Output is a 'parameters' data frame
# The 'chosen_parameter' column is altered to include values randomly chosen from the 'range_for_chosen_parameter'
alter_chosen_parameter_for_simulated_group <- function(parameters, 
                                                       chosen_parameter, 
                                                       range_for_chosen_parameter) {
  
  if (is.null(parameters[[chosen_parameter]])) {
    stop(paste(chosen_parameter, "is not a parameter in 'parameters'."))
  }
  
  n_simulations <- nrow(parameters)
  
  parameters[[chosen_parameter]] <- runif(n_simulations, 
                                          min = range_for_chosen_parameter["min"], 
                                          max = range_for_chosen_parameter["max"])
  
  if (chosen_parameter == "ringed_ring_width_factor") {
    if ("network_width" %in% colnames(parameters)) {
      parameters[["ringed_ring_width"]] <- parameters[["ringed_ring_width_factor"]] * parameters[["network_width"]]
    }
    else {
      parameters[["ringed_ring_width"]] <- parameters[["ringed_ring_width_factor"]] * 
        (parameters[["ellipsoid_x_radius"]] + parameters[["ellipsoid_y_radius"]] + parameters[["ellipsoid_z_radius"]]) / 3
    }
  }
  
  return(parameters)
}

# Output is a list of simulation metadata
# Simulation metadata is made using 'parameters' data frame as input
create_simulation_metadata_for_simulated_group <- function(parameters, 
                                                           shape, 
                                                           arrangement) {
  
  # Define fixed background metadata parameters
  background_parameters <- list(
    number_of_cells = 30000,
    length = 600, # Units: micrometers (um)
    width = 600,
    height = 300,
    minimum_distance_between_cells = 10,
    cell_types = c("A", "B", "O"),
    cell_proportions = c(0, 0, 1)
  )
  
  spe_metadata_background <- spe_metadata_background_template("random", NULL)
  spe_metadata_background$background$length <- background_parameters$length
  spe_metadata_background$background$width <- background_parameters$width
  spe_metadata_background$background$height <- background_parameters$height
  spe_metadata_background$background$n_cells <- background_parameters$number_of_cells
  spe_metadata_background$background$minimum_distance_between_cells <- background_parameters$minimum_distance_between_cells
  spe_metadata_background$background$cell_types <- background_parameters$cell_types
  spe_metadata_background$background$cell_proportions <- background_parameters$cell_proportions
  
  # Determine cluster_type
  if (arrangement == "mixed") {
    cluster_type <- "regular"
  }
  else if (arrangement == "ringed") {
    cluster_type <- "ring"
  }
  else if (arrangement == "separated") {
    cluster_type <- "regular"
  }
  else {
    stop(paste(arrangement, "is not a valid arrangement.", 
               "Should be 'mixed', 'ringed', 'separated'."))
  }
  
  # Add cluster metadata
  spe_metadata_cluster <- spe_metadata_cluster_template(cluster_type, shape, spe_metadata_background)
  spe_metadata_cluster$cluster_1$cluster_cell_types <- c('A', 'B')
  spe_metadata_cluster$cluster_1$centre_loc <- c(300, 300, 150)
  
  # Add sphere if arrangement is separated
  if (arrangement == "separated") {
    spe_metadata_cluster <- spe_metadata_cluster_template(cluster_type, "sphere", spe_metadata_cluster)
    spe_metadata_cluster$cluster_2$cluster_cell_types <- c('A', 'B')
    spe_metadata_cluster$cluster_2$cluster_cell_proportions <- c(0, 1)
    spe_metadata_cluster$cluster_2$radius <- 100
    spe_metadata_cluster$cluster_2$centre_loc <- c(450, 300, 150)
    
  }
  
  simulation_metadata <- list()
  
  for (index in seq_len(nrow(parameters))) {
    simulation_metadata[[index]] <- spe_metadata_cluster
    
    if (shape == "ellipsoid") {
      simulation_metadata[[index]][["cluster_1"]][["x_radius"]] <- parameters$ellipsoid_x_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["y_radius"]] <- parameters$ellipsoid_y_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["z_radius"]] <- parameters$ellipsoid_z_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["y_z_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
      simulation_metadata[[index]][["cluster_1"]][["x_z_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
      simulation_metadata[[index]][["cluster_1"]][["x_y_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
    }
    else if (shape == "network") {
      simulation_metadata[[index]][["cluster_1"]][["n_edges"]] <- 20
      simulation_metadata[[index]][["cluster_1"]][["width"]] <- parameters$network_width[index]
      simulation_metadata[[index]][["cluster_1"]][["radius"]] <- 125
    }
    else {
      stop(paste(shape, "is not a valid shape", 
                 "Should be 'ellipsoid', 'network'."))
    }
    
    
    ### Arrangements
    if (arrangement == "mixed") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- 
        c(parameters$mixed_cell_type_A_proportion[index], 1 - parameters$mixed_cell_type_A_proportion[index])
      
    }
    else if (arrangement == "ringed") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- c(1, 0)
      simulation_metadata[[index]][["cluster_1"]][["ring_cell_types"]] <- c('A', 'B')
      simulation_metadata[[index]][["cluster_1"]][["ring_cell_proportions"]] <- c(0, 1)
      simulation_metadata[[index]][["cluster_1"]][["ring_width"]] <- parameters$ringed_ring_width[index]
    }
    else if (arrangement == "separated") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- c(1, 0)
      simulation_metadata[[index]][["cluster_1"]][["centre_loc"]] <- 
        c(parameters$separated_cluster1_x_coordinate[index], 300, 150)
    }
    else {
      stop(paste(arrangement, "is not a valid arrangement.", 
                 "Should be 'mixed', 'ringed', 'separated'."))
    }  
  }
  return(simulation_metadata)
}


# Input are two lists of simulation metadata for two groups
# Output is complete analysis of these two groups (i.e. one set)
analyse_simulated_set <- function(simulation_metadata_for_simulated_group1,
                                  simulation_metadata_for_simulated_group2) {
  
  # Check
  if (length(simulation_metadata_for_simulated_group1) != length(simulation_metadata_for_simulated_group2)) {
    stop("Number of simulations in set 1 and set 2 are different.")
  }
  
  # Set defined parameters/values
  cell_types <- c('A', 'B')
  n_simulations <- length(simulation_metadata_for_simulated_group1)
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  n_splits <- 10
  thresholds <- seq(0.01, 1, 0.01)
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  # z-coords for slices
  bottom_z_coord_of_slices <- c(85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205)
  top_z_coord_of_slices <- bottom_z_coord_of_slices + 10
  n_slices <- length(bottom_z_coord_of_slices)
  
  # Function to get a random slice from spe
  get_random_slice_from_spe <- function(spe, 
                                        bottom_z_coord_of_slices, 
                                        top_z_coord_of_slices) {
    
    slices_from_spe <- list()
    
    if (length(bottom_z_coord_of_slices) != length(top_z_coord_of_slices)) stop("Lengths of bottom_z_coords_of_slices and top_z_coords_of_slices should be equal.")
    number_of_slices <- length(bottom_z_coord_of_slices)
    
    i <- sample(length(bottom_z_coord_of_slices), 1) # random index to choose from slice coordinates
    bottom_z_coord_of_random_slice <- bottom_z_coord_of_slices[i]
    top_z_coord_of_random_slice <- top_z_coord_of_slices[i]
    z_coords_of_cells_in_spe <- spatialCoords(spe)[ , "Cell.Z.Position"]
    random_slice_from_spe <- spe[ , bottom_z_coord_of_random_slice < z_coords_of_cells_in_spe & z_coords_of_cells_in_spe < top_z_coord_of_random_slice]
    spatialCoords(random_slice_from_spe) <- spatialCoords(random_slice_from_spe)[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    return(random_slice_from_spe)
  }
  
  # Function to create empty metric df list
  create_empty_metric_df_list <- function(
    cell_types,
    n_simulations,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("simulation", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           COO = COO_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  # Function analyse spe in 3D
  analyse_simulation3D <- function(spe, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   i, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (i - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- i
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for gradient-based data which increments after each target cell type
    pair_index <- n_cell_type_combinations * (i - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F,
                                                           feature_colname = "Cell.Type")
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                       paste(reference_cell_type, target_cell_type, sep = ","))
        

        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
          metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
        }
        else {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
          metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
          if (reference_cell_type != target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
            metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
          }
        }
        if (reference_cell_type == target_cell_type) {
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
        }
        
        # Spatial heterogeneity metrics
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        
        if (reference_cell_type != target_cell_type) {
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(spe, 
                                                                              n_splits,
                                                                              reference_cell_type, 
                                                                              target_cell_type,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
          
          if (is.null(proportion_grid_metrics)) {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            PBSAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                         "proportion",
                                                         weight_method = "queen")
            
            PBP_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                      "proportion",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
          } 
          
          entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe, 
                                                                   n_splits,
                                                                   c(reference_cell_type, target_cell_type), 
                                                                   plot_image = F,
                                                                   feature_colname = "Cell.Type")
          
          if (is.null(entropy_grid_metrics)) {
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            EBSAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                         "entropy",
                                                         weight_method = "queen")
            
            EBP_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                      "entropy",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
          }    
        }
        else {
          metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
          metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
          metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
          metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
        }
        
        pair_index <- pair_index + 1
      }
    }
    return(metric_df_list)
  }
  
  # Function analyse spe in 2D
  analyse_simulation2D <- function(spe, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   i, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (i - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- i
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for gradient-based data which increments after each target cell type
    pair_index <- n_cell_type_combinations * (i - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F,
                                                           feature_colname = "Cell.Type")
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                       paste(reference_cell_type, target_cell_type, sep = ","))
    
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
          metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
        }
        else {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
          metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
          if (reference_cell_type != target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
            metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
          }
        }
        if (reference_cell_type == target_cell_type) {
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
        }
        
        # Spatial heterogeneity metrics
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
        
        if (reference_cell_type != target_cell_type) {
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(spe, 
                                                                              n_splits,
                                                                              reference_cell_type, 
                                                                              target_cell_type,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
          
          if (is.null(proportion_grid_metrics)) {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                         "proportion",
                                                         weight_method = "queen")
            
            PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                      "proportion",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
          } 
          
          entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe, 
                                                                   n_splits,
                                                                   c(reference_cell_type, target_cell_type), 
                                                                   plot_image = F,
                                                                   feature_colname = "Cell.Type")
          
          if (is.null(entropy_grid_metrics)) {
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                         "entropy",
                                                         weight_method = "queen")
            
            EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                      "entropy",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
          }    
        }
        else {
          metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
          metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
          metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
          metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
        }
        
        pair_index <- pair_index + 1
      }
    }
    return(metric_df_list)
  }
  
  
  # Define metric df lists
  group1_metric_df_list3D <- group1_metric_df_list2D <- group2_metric_df_list3D <- group2_metric_df_list2D <-
    create_empty_metric_df_list(cell_types, n_simulations, radii_colnames, thresholds_colnames)
  
  for (i in seq(n_simulations)) {
    print(i)
    # Simulate spes in 3D, and get 2D slices
    group1_spe3D <- simulate_spe_metadata3D(simulation_metadata_for_simulated_group1[[i]], plot_image = F)
    group2_spe3D <- simulate_spe_metadata3D(simulation_metadata_for_simulated_group2[[i]], plot_image = F)
    
    group1_spe2D <- get_random_slice_from_spe(group1_spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)
    group2_spe2D <- get_random_slice_from_spe(group2_spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)

    # Analyse spes
    group1_metric_df_list3D <- analyse_simulation3D(group1_spe3D,
    cell_types,
    radii,
    thresholds,
    n_splits,
    i,
    group1_metric_df_list3D)
    
    group2_metric_df_list3D <- analyse_simulation3D(group2_spe3D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group2_metric_df_list3D)
    
    group1_metric_df_list2D <- analyse_simulation2D(group1_spe2D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group1_metric_df_list2D)
    
    group2_metric_df_list2D <- analyse_simulation2D(group2_spe2D,
                                                    cell_types,
                                                    radii,
                                                    thresholds,
                                                    n_splits,
                                                    i,
                                                    group2_metric_df_list2D)
  }
  
  simulated_set_analysis <- list(
    group1_metric_df_list3D = group1_metric_df_list3D,
    group2_metric_df_list3D = group2_metric_df_list3D,
    group1_metric_df_list2D = group1_metric_df_list2D,
    group2_metric_df_list2D = group2_metric_df_list2D
  )
  
  ## Turn gradient metrics into AUC and add to metric_df list
  get_AUC_for_radii_gradient_metrics <- function(y) {
    x <- radii
    h <- diff(x)[1]
    n <- length(x)
    
    AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
    
    return(AUC)
  }
  
  gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
  
  for (simulated_group_analysis in names(simulated_set_analysis)) {
    metric_df_list <- simulated_set_analysis[[simulated_group_analysis]]
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "reference", "target", metric_AUC_name)
      }
      else {
        subset_colnames <- c("simulation", "reference", metric_AUC_name)
      }
      
      df <- metric_df_list[[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      df <- df[ , subset_colnames]
      metric_df_list[[metric_AUC_name]] <- df
    }
    
    # PBP_AUC 3D
    PBP_df <- metric_df_list[["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    PBP_AUC_df <- PBP_df[ , c("simulation", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    simulated_set_analysis[[simulated_group_analysis]] <- metric_df_list
  }
  
  return(simulated_set_analysis)
}


# Output is p-values obtained from analysis of sets
# group1_chosen_parameters is a list of the form: list("mixed_cell_type_A_proportion" = c(min = 0.1, max = 0.5), ...)
analyse_simulated_collection <- function(n_simulations_in_group,
                                         n_sets,
                                         group1_shape,
                                         group1_arrangement,
                                         group1_chosen_parameters = NULL,
                                         group2_shape,
                                         group2_arrangement,
                                         group2_chosen_parameters)  {
  
  # Define constants
  metrics <- c("AMD",
               "MS_AUC", "NMS_AUC",
               "ACIN_AUC", "ANE_AUC", "ANC_AUC", "COO_AUC", "CK_AUC", "CL_AUC", "CG_AUC",
               "PBSAC", "EBSAC", "PBP_AUC", "EBP_AUC")
  
  cell_types <- c("A", "B")
  
  # Define result
  simulated_collection_p_values <- data.frame()
  
  # Subset metric_df
  subset_metric_df <- function(metric_df, reference_cell_type, target_cell_type, metric) {
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      metric_df <- metric_df[
        metric_df[["cell_types"]] == paste(reference_cell_type, target_cell_type, sep = ","), 
      ]
    }
    else if (metric %in% c("ANE_AUC")) {
      metric_df <- metric_df[
        metric_df[["target"]] == paste(reference_cell_type, target_cell_type, sep = ","), 
      ]
    }
    else {
      metric_df <- metric_df[
        metric_df[["reference"]] == reference_cell_type & metric_df[["target"]] == target_cell_type, 
      ]
    }
    return(metric_df)
  }
  
  # Calculate p-values from a single set
  calculate_p_values_from_simulated_set <- function(simulated_set_analysis,
                                                    metrics,
                                                    cell_types) {
    
    group1_metric_df_list3D <- simulated_set_analysis[["group1_metric_df_list3D"]]
    group2_metric_df_list3D <- simulated_set_analysis[["group2_metric_df_list3D"]]
    group1_metric_df_list2D <- simulated_set_analysis[["group1_metric_df_list2D"]]
    group2_metric_df_list2D <- simulated_set_analysis[["group2_metric_df_list2D"]]
    
    combined_p_value_df <- data.frame()
    
    for (metric in metrics) {
      metric_p_value_df <- data.frame()
      
      group1_metric_df3D <- group1_metric_df_list3D[[metric]]
      group2_metric_df3D <- group2_metric_df_list3D[[metric]]
      group1_metric_df2D <- group1_metric_df_list2D[[metric]]
      group2_metric_df2D <- group2_metric_df_list2D[[metric]]

      for (reference_cell_type in cell_types) {
        for (target_cell_type in cell_types) {
          
          group1_metric_df3D_subset <- subset_metric_df(group1_metric_df3D, reference_cell_type, target_cell_type, metric)
          group2_metric_df3D_subset <- subset_metric_df(group2_metric_df3D, reference_cell_type, target_cell_type, metric)
          group1_metric_df2D_subset <- subset_metric_df(group1_metric_df2D, reference_cell_type, target_cell_type, metric)
          group2_metric_df2D_subset <- subset_metric_df(group2_metric_df2D, reference_cell_type, target_cell_type, metric)
          
          p_value3D <- wilcox.test(group1_metric_df3D_subset[[metric]], group2_metric_df3D_subset[[metric]])$p.value
          p_value2D <- wilcox.test(group1_metric_df2D_subset[[metric]], group2_metric_df2D_subset[[metric]])$p.value
          
          # Add rows for 3D and 2D results
          metric_p_value_df <- rbind(
            metric_p_value_df,
            data.frame(
              reference = reference_cell_type,
              target    = target_cell_type,
              metric    = metric,
              dimension = "3D",
              p_value   = p_value3D
            ),
            data.frame(
              reference = reference_cell_type,
              target    = target_cell_type,
              metric    = metric,
              dimension = "2D",
              p_value   = p_value2D
            )
          )
        }
      }
      combined_p_value_df <- rbind(combined_p_value_df, metric_p_value_df)
    }
    return(combined_p_value_df)
  }
  
  
  # Iterate for each set
  for (set_index in seq_len(n_sets)) {
    print(paste("Set:", set_index))
    group1_parameters <- create_random_parameters_for_simulated_group(n_simulations = n_simulations_in_group,
                                                                      shape = group1_shape,
                                                                      arrangement = group1_arrangement)
    
    group2_parameters <- create_random_parameters_for_simulated_group(n_simulations = n_simulations_in_group,
                                                                      shape = group2_shape,
                                                                      arrangement = group2_arrangement)
    
    if (!is.null(group1_chosen_parameters)) {
      for (i in seq_len(length(group1_chosen_parameters))) {
        group1_parameters <- alter_chosen_parameter_for_simulated_group(parameters = group1_parameters,
                                                                        chosen_parameter = names(group1_chosen_parameters)[i],
                                                                        range_for_chosen_parameter = group1_chosen_parameters[i])
      }
    }
    
    if (!is.null(group2_chosen_parameters)) {
      for (i in seq_len(length(group2_chosen_parameters))) {
        group2_parameters <- alter_chosen_parameter_for_simulated_group(parameters = group2_parameters,
                                                                        chosen_parameter = names(group2_chosen_parameters)[i],
                                                                        range_for_chosen_parameter = group2_chosen_parameters[[i]])
      }
    }
    
    group1_simulation_metadata <- create_simulation_metadata_for_simulated_group(parameters = group1_parameters,
                                                                                 shape = group1_shape,
                                                                                 arrangement = group1_arrangement)
    
    group2_simulation_metadata <- create_simulation_metadata_for_simulated_group(parameters = group2_parameters,
                                                                                 shape = group2_shape,
                                                                                 arrangement = group2_arrangement)  
    
    simulated_set_analysis <- analyse_simulated_set(group1_simulation_metadata,
                                                    group2_simulation_metadata)
    
    simulated_set_p_values <- calculate_p_values_from_simulated_set(simulated_set_analysis,
                                                                    metrics,
                                                                    cell_types)
    
    simulated_set_p_values$set <- set_index
    
    simulated_collection_p_values <- rbind(simulated_collection_p_values, simulated_set_p_values)
    
  }
  return(simulated_collection_p_values)
}


