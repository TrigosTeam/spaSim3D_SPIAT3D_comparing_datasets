# Code to show comparison of stomics fly (drosophila) datasets.
# Requires analysis of each fly sample using function in 'spaSim-3D_SPIAT-3D_analysing_individual_datasets' repo.

### Libraries -----
library(SpatialExperiment)
library(plotly)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(scales)
library(DescTools)
### Read analysis data -------
setwd("~/R/public_data_analysis")

stomics_dataset_names <- c(
  "E14-16h",
  "E16-18h",
  "L1",
  "L2",
  "L3"
)


stomics_analysis_list <- list()

for (stomics_dataset_name in stomics_dataset_names) {
  stomics_analysis_list[[stomics_dataset_name]] <- readRDS(paste("stomics_", stomics_dataset_name, "_metric_df_list.RDS", sep = ""))
}


### Format analysis data by adding AUC --------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACIN", "ANE", "ANC", "FCO", "CK", "CL", "CG")) {
    return("radius")
  }
  else if (metric %in% c("PBP", "EBP")) {
    return("threshold")  
  }
  else {
    stop("Invalid metric. Must be gradient-based")
  }
}


## Turn gradient radii metrics into AUC and add to metric_df list
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}

radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "FCO", "CK", "CL", "CG")

## Turn threshold radii metrics into AUC and add to metric_df list
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")


for (stomics_dataset_name in stomics_dataset_names) {
  metric_df_list <- stomics_analysis_list[[stomics_dataset_name]]
  
  for (metric in gradient_radii_metrics) {
    metric_AUC_name <- paste(metric, "AUC", sep = "_")
    
    if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "FCO", "CK", "CL", "CG")) {
      subset_colnames <- c("slice", "reference", "target", metric_AUC_name)
    }
    else {
      subset_colnames <- c("slice", "reference", metric_AUC_name)
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
  PBP_AUC_df <- PBP_df[ , c("slice", "reference", "target", "PBP_AUC")]
  metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
  
  # EBP_AUC 3D
  EBP_df <- metric_df_list[["EBP"]]
  EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
    sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
  })
  EBP_AUC_df <- EBP_df[ , c("slice", "cell_types", "EBP_AUC")]
  metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
  
  stomics_analysis_list[[stomics_dataset_name]] <- metric_df_list
}



### Format analysis data by combining metric_df_lists from different datasets, and changing 3D slice index to 0 --------

metrics <- c("AMD",
             "MS_AUC", "NMS_AUC",
             "ACIN_AUC", "ANE_AUC", "ANC_AUC", "FCO_AUC", "CK_AUC", "CL_AUC", "CG_AUC",
             "PBSAC", "EBSAC", "PBP_AUC", "EBP_AUC")


metric_df_list_combined <- list()

for (metric in metrics) {
  
  metric_df_list_combined[[metric]] <- data.frame()
  
  for (stomics_dataset_name in stomics_dataset_names) {
    
    metric_df_list <- stomics_analysis_list[[stomics_dataset_name]]
    
    df <- metric_df_list[[metric]]
    
    # Ensure slice column is numeric, not character
    df$slice <- as.numeric(df$slice) 
    
    # Change 3D slice index to 0
    df$slice[df$slice == max(df$slice)] <- 0 
    
    # Add stomics dataset column to each metric df
    df$dataset <- stomics_dataset_name 
    
    # Change column name containing metric values to 'values' for all metrics
    colnames(df)[colnames(df) == metric] <- "value"
    
    metric_df_list_combined[[metric]] <- rbind(metric_df_list_combined[[metric]], df)
  }
}



### Compare 3D and 2D functions (all metrics) ----

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

plot_metric_3D_vs_2D_comparison <- function(metric_3D_vs_2D_comparison_results,
                                            cell_types_of_interest,
                                            metrics) {
  
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  figs <- list()
  
  for (metric in metrics) {
    metric_data <- metric_3D_vs_2D_comparison_results[[metric]]
    
    metric_data2D <- metric_data[metric_data$slice != 0, ]
    metric_data3D <- metric_data[metric_data$slice == 0, ]
    
    figs[[metric]] <- list()
    
    for (reference_cell_type in cell_types_of_interest) {
      for (target_cell_type in cell_types_of_interest) {
        
        metric_data2D_subset <- subset_metric_df(metric_data2D, reference_cell_type, target_cell_type, metric)
        metric_data3D_subset <- subset_metric_df(metric_data3D, reference_cell_type, target_cell_type, metric)
        
        metric_data2D_subset$value <- (metric_data2D_subset$value) - (median(metric_data2D_subset[metric_data2D_subset[["dataset"]] == "E14-16h", "value"], na.rm = TRUE))
        metric_data3D_subset$value <- (metric_data3D_subset$value) - (metric_data3D_subset[metric_data3D_subset[["dataset"]] == "E14-16h", "value"])
        
        fig <- ggplot(metric_data2D_subset, aes(x = dataset, y = value)) +
          geom_boxplot(outlier.shape = NA, fill = "lightgray") +
          geom_jitter(width = 0.2, color = "#0062c5", alpha = 0.5) +
          geom_point(data = metric_data3D_subset, aes(x = dataset, y = value),
                     color = "#bb0036", shape = 8, size = 3) +
          labs(x = "Fly stage", y = sub("_AUC$", "", metric), title = paste(reference_cell_type, target_cell_type, sep = "/")) +
          theme_minimal() +
          theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            
            # Center and size title 
            plot.title = element_text(size = 14, hjust = 0.5),
            
            # Axis title sizes
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            
            # Axis tick label sizes
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 14)
          ) +
          scale_y_continuous(
            breaks = pretty_breaks(n = 3),
            labels = sci_clean_threshold
          )
        
        
        # Add line that connects each box plot's median value
        tryCatch({
          median_df <- aggregate(value ~ dataset, data = metric_data2D_subset, FUN = median)
          
          fig <- fig + geom_line(data = median_df, aes(x = dataset, y = value, group = 1),
                                 color = "black", linetype = "solid", linewidth = 0.5)
        }, error = function(e) {
          message("Skipping median line due to error: ", e$message)
        })
        
        # Add line that connects each 3D red star value
        fig <- fig + geom_line(data = metric_data3D_subset, aes(x = dataset, y = value, group = 1),
                               color = "#bb0036", linetype = "solid", linewidth = 0.5)
        
        pair <- paste(reference_cell_type, target_cell_type, sep = '/')
        figs[[metric]][[pair]] <- fig
        
      }
    }
  }
  return(figs)
}

plot_metric_3D_vs_2D_comparison_no_box_plot <- function(metric_3D_vs_2D_comparison_results,
                                                        cell_types_of_interest,
                                                        metrics) {
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  figs <- list()
  
  for (metric in metrics) {
    metric_data <- metric_3D_vs_2D_comparison_results[[metric]]
    
    metric_data2D <- metric_data[metric_data$slice != 0, ]
    metric_data3D <- metric_data[metric_data$slice == 0, ]
    
    figs[[metric]] <- list()
    
    for (reference_cell_type in cell_types_of_interest) {
      for (target_cell_type in cell_types_of_interest) {
        
        metric_data2D_subset <- subset_metric_df(metric_data2D, reference_cell_type, target_cell_type, metric)
        metric_data3D_subset <- subset_metric_df(metric_data3D, reference_cell_type, target_cell_type, metric)
        
        metric_data2D_subset$value <- (metric_data2D_subset$value) - (median(metric_data2D_subset[metric_data2D_subset[["dataset"]] == "E14-16h", "value"], na.rm = TRUE))
        metric_data3D_subset$value <- (metric_data3D_subset$value) - (metric_data3D_subset[metric_data3D_subset[["dataset"]] == "E14-16h", "value"])
        
        # Compute medians
        medians <- metric_data2D_subset %>%
          group_by(dataset) %>%
          summarise(median_value = median(value, na.rm = TRUE)) %>%
          filter(!is.na(median_value))
        
        fig <- ggplot(metric_data2D_subset, aes(x = dataset, y = value)) +
          geom_jitter(width = 0.2, color = "#0062c5", alpha = 0.5) +
          geom_point(data = metric_data3D_subset, aes(x = dataset, y = value),
                     color = "#bb0036", shape = 8, size = 3) +
          geom_crossbar(data = medians, aes(x = dataset, y = median_value, ymin = median_value, ymax = median_value),
                        width = 0.5, color = "black") +
          labs(x = "Fly stage", y = sub("_AUC$", "", metric), title = paste(reference_cell_type, target_cell_type, sep = "/")) +
          theme_minimal() +
          theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            
            # Center and size title 
            plot.title = element_text(size = 14, hjust = 0.5),
            
            # Axis title sizes
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            
            # Axis tick label sizes
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 14)
          ) +
          scale_y_continuous(
            breaks = pretty_breaks(n = 3),
            labels = sci_clean_threshold
          )
        
        
        # Add line that connects each box plot's median value
        tryCatch({
          median_df <- aggregate(value ~ dataset, data = metric_data2D_subset, FUN = median)
          
          fig <- fig + geom_line(data = median_df, aes(x = dataset, y = value, group = 1),
                                 color = "black", linetype = "solid", linewidth = 0.5)
        }, error = function(e) {
          message("Skipping median line due to error: ", e$message)
        })
        
        # Add line that connects each 3D red star value
        fig <- fig + geom_line(data = metric_data3D_subset, aes(x = dataset, y = value, group = 1),
                               color = "#bb0036", linetype = "solid", linewidth = 0.5)
        
        pair <- paste(reference_cell_type, target_cell_type, sep = '/')
        figs[[metric]][[pair]] <- fig
        
      }
    }
  }
  return(figs)
}

tally_number_of_3D_values_in_and_above_and_below_2D_box_plot <- function(metric_3D_vs_2D_comparison_results,
                                                                         cell_types_of_interest,
                                                                         metrics) {
  
  check_if_value_is_between_q1_and_q3 <- function(value, num_vector) {
    
    q1 <- quantile(num_vector, 0.25, na.rm = TRUE)
    q3 <- quantile(num_vector, 0.75, na.rm = TRUE)
    
    if (length(value) == 0 || is.na(q1) || is.na(q3)) {
      return(F)
    }
    
    if (q1 < value && value < q3) {
      return(T)
    }
    else {
      return(F)
    }
  }
  
  check_if_value_is_above_q3 <- function(value, num_vector) {
    
    q3 <- quantile(num_vector, 0.75, na.rm = TRUE)
    
    if (length(value) == 0 || is.na(q3)) {
      return(F)
    }
    
    if (value > q3) {
      return(T)
    }
    else {
      return(F)
    }
  }
  
  check_if_value_is_below_q1 <- function(value, num_vector) {
    
    q1 <- quantile(num_vector, 0.25, na.rm = TRUE)
    
    if (length(value) == 0 || is.na(q1)) {
      return(F)
    }
    
    if (value < q1) {
      return(T)
    }
    else {
      return(F)
    }
  }
  
  tallies <- list()
  
  fly_stages <- c("E14-16h", "E16-18h", "L1", "L2", "L3" )
  
  for (metric in metrics) {
    metric_data <- metric_3D_vs_2D_comparison_results[[metric]]
    
    metric_data2D <- metric_data[metric_data$slice != 0, ]
    metric_data3D <- metric_data[metric_data$slice == 0, ]
    
    tallies[[metric]] <- list()
    
    for (reference_cell_type in cell_types_of_interest) {
      for (target_cell_type in cell_types_of_interest) {
        
        metric_data2D_subset <- subset_metric_df(metric_data2D, reference_cell_type, target_cell_type, metric)
        metric_data3D_subset <- subset_metric_df(metric_data3D, reference_cell_type, target_cell_type, metric)
        
        metric_data2D_subset$value <- (metric_data2D_subset$value) - (median(metric_data2D_subset[metric_data2D_subset[["dataset"]] == "E14-16h", "value"], na.rm = TRUE))
        metric_data3D_subset$value <- (metric_data3D_subset$value) - (metric_data3D_subset[["value"]])[1]
        
        tally <- 0
        
        for (fly_stage in fly_stages) {
          
          # Ignore the first fly stage which is always correct.
          if (fly_stage == "E14-16h") {
            next
          }
          
          # Add to 'ones' if value is below box plot
          if (check_if_value_is_below_q1(
            metric_data3D_subset[metric_data3D_subset[["dataset"]] == fly_stage, "value"],
            metric_data2D_subset[metric_data2D_subset[["dataset"]] == fly_stage, "value"])) {
            tally <- tally + 1
          }
          # Add to 'tens' if value is in box plot
          if (check_if_value_is_between_q1_and_q3(
            metric_data3D_subset[metric_data3D_subset[["dataset"]] == fly_stage, "value"],
            metric_data2D_subset[metric_data2D_subset[["dataset"]] == fly_stage, "value"])) {
            tally <- tally + 10
          }
          # Add to 'hundreds' if value is above box plot
          if (check_if_value_is_above_q3(
            metric_data3D_subset[metric_data3D_subset[["dataset"]] == fly_stage, "value"],
            metric_data2D_subset[metric_data2D_subset[["dataset"]] == fly_stage, "value"])) {
            tally <- tally + 100
          }
        }
        pair <- paste(reference_cell_type, target_cell_type, sep = '/')
        tallies[[metric]][[pair]] <- tally
      }
    }
  }
  return(tallies)
}

plot_tallies_with_categories <- function(tallies, metrics) {
  
  figs <- list()
  
  for (metric in metrics) {
    
    tallies_subset <- tallies[[metric]]
    
    plot_df <- stack(tallies_subset)
    
    split_vec <- strsplit(as.character(plot_df$ind), '/')
    plot_df$reference <- sapply(split_vec, `[`, 1)
    plot_df$target <- sapply(split_vec, `[`, 2)
    mat <- matrix(nrow = 6, ncol = 6)
    rownames(mat) <- c("CNS", "epidermis", "carcass", "muscle", "midgut", "fat body")
    colnames(mat) <- c("CNS", "epidermis", "carcass", "muscle", "midgut", "fat body")
    
    for (reference in unique(plot_df$reference)) {
      for (target in unique(plot_df$target)) {
        
        mat[reference, target] <- plot_df$values[plot_df$ind == paste(reference, target, sep = "/")]
      }
    }
    
    # Function to classify each number
    classify_value <- function(x) {
      
      if (x == 0) {
        return("N/A")
      }
      
      hundreds <- x %/% 100
      tens <- (x %% 100) %/% 10
      ones <- x %% 10
      
      if (ones >= 3) {
        return("over")
      } else if (tens >= 3) {
        return("consistent")
      } else if (hundreds >= 3) {
        return("under")
      } else {
        return("inconsistent")
      }
    }
    # Apply the function to every element in the matrix
    classified_mat <- matrix(
      data = sapply(mat, classify_value),
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
    
    # Define custom colors for each label
    label_colors <- c(
      "over" = "#73ec81",
      "under" = "#b8db50",
      "consistent" = "#004a07",
      "inconsistent" = "#bb0036",
      "N/A" = "gray"
    )
    
    fig <- Heatmap(classified_mat, 
                   column_title = sub("_AUC$", "", metric),
                   name = "label",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   col = label_colors)
    
    figs[[metric]] <- fig
  }
  return(figs)
}

plot_tallies_with_categories_bar <- function(tallies, metrics) {
  
  classified_df <- data.frame()
  
  for (metric in metrics) {
    tallies_subset <- tallies[[metric]]
    
    plot_df <- stack(tallies_subset)
    
    split_vec <- strsplit(as.character(plot_df$ind), '/')
    plot_df$reference <- sapply(split_vec, `[`, 1)
    plot_df$target <- sapply(split_vec, `[`, 2)
    
    mat <- matrix(nrow = 6, ncol = 6)
    rownames(mat) <- c("CNS", "epidermis", "carcass", "muscle", "midgut", "fat body")
    colnames(mat) <- c("CNS", "epidermis", "carcass", "muscle", "midgut", "fat body")
    
    for (reference in unique(plot_df$reference)) {
      for (target in unique(plot_df$target)) {
        mat[reference, target] <- plot_df$values[plot_df$ind == paste(reference, target, sep = "/")]
      }
    }
    
    # Function to classify each number
    classify_value <- function(x) {
      
      if (x == 0) {
        return("N/A")
      }
      
      hundreds <- x %/% 100
      tens <- (x %% 100) %/% 10
      ones <- x %% 10
      
      if (ones >= 3) {
        return("over")
      } else if (tens >= 3) {
        return("consistent")
      } else if (hundreds >= 3) {
        return("under")
      } else {
        return("inconsistent")
      }
    }
    # Build the classification dataframe
    classified_df_subset <- data.frame(
      metric = character(),
      category = character(),
      pair = character(),
      stringsAsFactors = FALSE
    )
    
    for (i in rownames(mat)) {
      for (j in colnames(mat)) {
        value <- mat[i, j]
        if (!is.na(value)) {
          classified_df_subset <- rbind(classified_df_subset, data.frame(
            metric = metric, # Change later
            category = classify_value(value),
            pair = paste(i, j, sep = "/")
          ))
        }
      }
    }
    
    classified_df <- rbind(classified_df, classified_df_subset)
    
    # Define custom colors for each label
    label_colors <- c(
      "over" = "#73ec81",
      "under" = "#b8db50",
      "consistent" = "#004a07",
      "inconsistent" = "#bb0036"
    )
    
  }
  # Plot proportions
  classified_df <- subset(classified_df, category != "N/A")
  
  # Remove AUC from metric names
  classified_df$metric <- sub("_AUC$", "", classified_df$metric)
  
  # Create a proportions_df to save values
  proportions_df <- classified_df %>%
    group_by(metric, category) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = category,
      values_from = n,
      values_fill = 0
    ) %>%
    rowwise()
  
  proportions_df$total <- 
    as.numeric(proportions_df$consistent) + 
    as.numeric(proportions_df$over) + 
    as.numeric(proportions_df$under) + 
    as.numeric(proportions_df$inconsistent)
  
  proportions_df <- proportions_df %>%
    mutate(
      over = over / total,
      under = under / total,
      consistent = consistent / total,
      inconsistent = inconsistent / total
    ) %>%
    ungroup() %>%
    select(-total) %>%
    arrange(
      desc(inconsistent), 
      desc(over),
      desc(under),
      desc(consistent)
    )
  write.csv(proportions_df, "~/R/values_from_figures/fly_proportion_values.csv")
  
  # Get metric order from here
  metric_order <- proportions_df %>% 
    mutate(
      over_under = over + under
    ) %>% 
    arrange(
      desc(inconsistent), 
      desc(over_under),
      desc(consistent)
    ) %>%
    pull(metric)
  
  # Factor metric
  classified_df$metric <- factor(classified_df$metric, levels = metric_order)
  
  # Factor categories
  classified_df$category <- factor(
    classified_df$category,
    levels = c("inconsistent", "over", "under", "consistent")
  )
  
  fig <- ggplot(classified_df, aes(x = metric, fill = category)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = label_colors) +
    labs(
      x = "Metric",
      y = "Proportion",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(
      title = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 16, angle = 90, hjust = 1),
      axis.text.y  = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  
  
  return(fig)
}

plot_metric_3D_vs_2D_comparison_categorised <- function(metric_3D_vs_2D_comparison_results,
                                                        cell_types_of_interest,
                                                        metrics,
                                                        tallies) {
  
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  figs <- list()
  
  for (metric in metrics) {
    metric_data <- metric_3D_vs_2D_comparison_results[[metric]]
    
    metric_data2D <- metric_data[metric_data$slice != 0, ]
    metric_data3D <- metric_data[metric_data$slice == 0, ]
    
    figs[[metric]] <- list()
    
    # Subset tallies for current metric
    tallies_subset <- tallies[[metric]]
    
    tallies_subset <- stack(tallies_subset)
    
    split_vec <- strsplit(as.character(tallies_subset$ind), '/')
    tallies_subset$reference <- sapply(split_vec, `[`, 1)
    tallies_subset$target <- sapply(split_vec, `[`, 2)
    
    # Function to classify each number
    classify_value <- function(x) {
      
      if (x == 0) {
        return("N/A")
      }
      
      hundreds <- x %/% 100
      tens <- (x %% 100) %/% 10
      ones <- x %% 10
      
      if (ones >= 3) {
        return("over")
      } else if (tens >= 3) {
        return("consistent")
      } else if (hundreds >= 3) {
        return("under")
      } else {
        return("inconsistent")
      }
    }
    
    tallies_subset$category <- Vectorize(classify_value)(tallies_subset$values)
    
    for (category in c("consistent", "inconsistent", "over", "under")) {
      
      figs[[metric]][[category]] <- list()
      
      for (reference_cell_type in cell_types_of_interest) {
        for (target_cell_type in cell_types_of_interest) {
          
          if (tallies_subset[tallies_subset$reference == reference_cell_type & tallies_subset$target == target_cell_type, "category"] != category) {
            next
          }
          
          metric_data2D_subset <- subset_metric_df(metric_data2D, reference_cell_type, target_cell_type, metric)
          metric_data3D_subset <- subset_metric_df(metric_data3D, reference_cell_type, target_cell_type, metric)
          
          metric_data2D_subset$value <- (metric_data2D_subset$value) - (median(metric_data2D_subset[metric_data2D_subset[["dataset"]] == "E14-16h", "value"], na.rm = TRUE))
          metric_data3D_subset$value <- (metric_data3D_subset$value) - (metric_data3D_subset[metric_data3D_subset[["dataset"]] == "E14-16h", "value"])
          
          fig <- ggplot(metric_data2D_subset, aes(x = dataset, y = value)) +
            geom_boxplot(outlier.shape = NA, fill = "lightgray") +
            geom_jitter(width = 0.2, color = "#0062c5", alpha = 0.5) +
            geom_point(data = metric_data3D_subset, aes(x = dataset, y = value),
                       color = "#bb0036", shape = 8, size = 3) +
            labs(x = "Fly stage", y = metric, title = paste(reference_cell_type, target_cell_type, sep = "/")) +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              
              # Center and size title 
              plot.title = element_text(size = 14, hjust = 0.5),
              
              # Axis title sizes
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              
              # Axis tick label sizes
              axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 14)
            ) +
            scale_y_continuous(
              breaks = pretty_breaks(n = 3),
              labels = sci_clean_threshold
            )
          
          
          # Add line that connects each box plot's median value
          tryCatch({
            median_df <- aggregate(value ~ dataset, data = metric_data2D_subset, FUN = median)
            
            fig <- fig + geom_line(data = median_df, aes(x = dataset, y = value, group = 1),
                                   color = "black", linetype = "solid", linewidth = 0.5)
          }, error = function(e) {
            message("Skipping median line due to error: ", e$message)
          })
          
          # Add line that connects each 3D red star value
          fig <- fig + geom_line(data = metric_data3D_subset, aes(x = dataset, y = value, group = 1),
                                 color = "#bb0036", linetype = "solid", linewidth = 0.5)
          
          pair <- paste(reference_cell_type, target_cell_type, sep = '/')
          figs[[metric]][[category]][[pair]] <- fig
          
        }
      }
    }
  }
  return(figs)
}

### Compare 3D and 2D analysis -----
metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "FCO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


stomics_fly_cells <- c("CNS", "epidermis", "carcass", "muscle", "midgut", "fat body")


# 3D vs 2D comparison
fig_metric_3D_vs_2D_comparison <- plot_metric_3D_vs_2D_comparison(metric_df_list_combined, stomics_fly_cells, metrics)
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison.pdf", width = 25, height = 15)

for (metric in metrics) {
  grid.arrange(do.call(arrangeGrob, c(fig_metric_3D_vs_2D_comparison[[metric]], nrow = 6, ncol = 6)))
}

dev.off()


# 3D vs 2D comparison (without box plot)
fig_metric_3D_vs_2D_comparison_no_box_plot <- plot_metric_3D_vs_2D_comparison_no_box_plot(metric_df_list_combined, stomics_fly_cells, metrics)
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_no_box_plot.pdf", width = 25, height = 15)

for (metric in metrics) {
  grid.arrange(do.call(arrangeGrob, c(fig_metric_3D_vs_2D_comparison_no_box_plot[[metric]], nrow = 6, ncol = 6)))
}

dev.off()


# Get tallies for categorical trends
tallies <- tally_number_of_3D_values_in_and_above_and_below_2D_box_plot(metric_df_list_combined,
                                                                        stomics_fly_cells,
                                                                        metrics)


# Heatmaps of tallies
tallies_plot1 <- plot_tallies_with_categories(tallies, metrics)
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_heatmaps.pdf", width = 6, height = 5)

for (metric in metrics) {
  fig <- tallies_plot1[[metric]]
  plot(fig)
}

dev.off()


# Bar plots of tallies
tallies_plot2 <- plot_tallies_with_categories_bar(tallies, metrics)
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_bar.pdf", width = 8, height = 3.6)
print(tallies_plot2)
dev.off()


# Box plots of examples of different cateogrical trends
# Consistent
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_consistent_example_MS.pdf", width = 4.2, height = 2.8)
print(fig_metric_3D_vs_2D_comparison$MS_AUC$`muscle/CNS`)
dev.off()

# Over
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_over_example_ACIN.pdf", width = 4.2, height = 2.8)
print(fig_metric_3D_vs_2D_comparison$ACIN_AUC$`fat body/carcass`)
dev.off()

# Under
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_under_example_AMD.pdf", width = 4.2, height = 2.8)
print(fig_metric_3D_vs_2D_comparison$AMD$`epidermis/carcass`)
dev.off()

# Inconsistent
setwd("~/R/plots/public_data/stomics_comparison")
pdf("stomics_3D_vs_2D_comparison_inconsistent_example_CL.pdf", width = 4.2, height = 2.8)
print(fig_metric_3D_vs_2D_comparison$CL_AUC$`midgut/fat body`)
dev.off()




# 3D vs 2D comparison, separated into categorical trends
# fig_metric_3D_vs_2D_comparison_categorised <- plot_metric_3D_vs_2D_comparison_categorised(metric_df_list_combined, stomics_fly_cells, metrics, tallies)
# 
# setwd("~/R/plots/public_data/stomics_comparison")
# 
# for (category in c("consistent", "inconsistent", "over", "under")) {
#   
#   pdf(paste("stomics_3D_vs_2D_comparison_", category, ".pdf", sep = ""), width = 25, height = 15)
#   
#   for (metric in metrics) {
#     
#     if (length(fig_metric_3D_vs_2D_comparison_categorised[[metric]][[category]]) == 0) {
#       next 
#     }
#     
#     grid.arrange(do.call(arrangeGrob, c(fig_metric_3D_vs_2D_comparison_categorised[[metric]][[category]])))
#   }
#   
#   dev.off() 
# }

