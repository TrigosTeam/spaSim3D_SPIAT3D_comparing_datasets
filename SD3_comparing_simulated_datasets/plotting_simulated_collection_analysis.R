# Code to plot simulated collection analysis

# Libraries ----
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


# Functions -----
# *** Function 1. ***
# Show two box plots, one for 2D and other for 3D, containing p-values from all metrics and cell combination pairs
box_plot_p_values_3D_vs_2D <- function(simulated_collection_analysis) {
  
  # Create the plot
  fig <- ggplot(data = simulated_collection_analysis, 
                aes(x = dimension, y = p_value)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +
    geom_jitter(width = 0.2, color = "#0062c5", alpha = 0.5) +
    geom_hline(yintercept = (0.05), color = "#bb0036", linetype = "dashed") + # horizontal line at p = 0.05
    labs(title = "All p-values vs dimension",
         x = "Dimension",
         y = "p-value") +
    ylim(0, 1) +
    stat_summary(
      fun = median,
      geom = "text",
      aes(label = round(..y.., 3)), # Display median values rounded to 3 decimals
      color = "black",
      vjust = -0.5 # Adjust vertical position of text
    ) +
    theme_minimal()
  
  return(fig)
}

# *** Function 2. ***
# Each metric gets a single main graph, with sub-graphs for each reference-target pair.
# Each sub-graph has two box plots, one for 2D and the other for 3D, similar to above box plot function.
box_plot_p_values_for_each_metric_and_pair_3D_vs_2D <- function(simulated_collection_analysis, metrics) {
  
  # Add a combined 'pair' column for reference and target
  simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                              simulated_collection_analysis$target,
                                              sep = "/")
  
  # Empty list to carry all the figs for each metric
  figs <- list()
  
  for (metric in metrics) {
    # Filter table for current metric
    simulated_collection_analysis_metric_subset <- 
      simulated_collection_analysis[simulated_collection_analysis$metric == metric, ]
    
    # Create the plot
    fig <- ggplot(data = simulated_collection_analysis_metric_subset, 
                  aes(x = dimension, y = p_value)) +
      geom_boxplot(outlier.shape = NA, fill = "lightgray") +
      geom_jitter(width = 0.2, color = "#0062c5", alpha = 0.5) +
      geom_hline(yintercept = (0.05), color = "#bb0036", linetype = "dashed") + # horizontal line at p = 0.05
      facet_wrap(~pair, scales = "free") +
      labs(title = sub("_AUC$", "", metric),     # Remove AUC from metrics
           x = "Dimension",
           y = "p-value") +
      ylim(0, 1) +
      # stat_summary(
      #   fun = median,
      #   geom = "text",
      #   aes(label = round(..y.., 3)), # Display median values rounded to 3 decimals
      #   color = "#0062c5",
      #   vjust = -0.5 # Adjust vertical position of text
      # ) +
      theme_minimal()
    
    figs[[metric]] <- fig
  }
  return(figs)
}

# *** Function 3. ***
# Same as function 2, but log of p_values
box_plot_log_p_values_for_each_metric_and_pair_3D_vs_2D <- function(simulated_collection_analysis, metrics) {
  
  # Add a combined 'pair' column for reference and target
  simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                              simulated_collection_analysis$target,
                                              sep = "/")
  
  # Empty list to carry all the figs for each metric
  figs <- list()
  
  for (metric in metrics) {
    # Filter table for current metric
    simulated_collection_analysis_metric_subset <- 
      simulated_collection_analysis[simulated_collection_analysis$metric == metric, ]
    
    # Create the plot
    fig <- ggplot(data = simulated_collection_analysis_metric_subset, 
                  aes(x = dimension, y = log(p_value))) +
      geom_boxplot(outlier.shape = NA, fill = "lightgray") +
      geom_jitter(width = 0.2, color = "#0062c5", alpha = 0.5) +
      geom_hline(yintercept = log(0.05), color = "#bb0036", linetype = "dashed") + # horizontal line at p = 0.05
      facet_wrap(~pair, scales = "free") +
      labs(title = sub("_AUC$", "", metric),     # Remove AUC from metrics
           x = "Dimension",
           y = "log(p-value)") +
      ylim(log(0.001), log(1)) +
      # stat_summary(
      #   fun = median,
      #   geom = "text",
      #   aes(label = round(..y.., 3)), # Display median values rounded to 3 decimals
      #   color = "#0062c5",
      #   vjust = -0.5 # Adjust vertical position of text
      # ) +
      theme_minimal()
    
    figs[[metric]] <- fig
  }
  return(figs)
}

# *** Function 4. ***
# Create a heatmap with metrics on the rows and cell type combination on the columns.
# p-value false positive or false negative rate shown
heatmap_p_values_for_each_metric_and_pair_3D_vs_2D <- function(simulated_collection_analysis, 
                                                               metrics, 
                                                               threshold) {
  
  # Add a combined 'pair' column for reference and target
  simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                              simulated_collection_analysis$target,
                                              sep = "/")
  
  # Get cell type combinations pairs
  pairs <- unique(simulated_collection_analysis$pair)
  
  # Set up metrics vs cell type combinations matrix
  metrics_vs_pairs_matrix3D <- matrix(nrow = length(pairs), ncol = length(metrics))
  rownames(metrics_vs_pairs_matrix3D) <- pairs
  colnames(metrics_vs_pairs_matrix3D) <- metrics
  
  # 2D matrix has same format as 3D matrix
  metrics_vs_pairs_matrix2D <- metrics_vs_pairs_matrix3D
  
  # Iterate through each metric and cell type combination
  for (pair in pairs) {
    for (metric in metrics) {
      
      simulated_collection_analysis_filtered <- 
        simulated_collection_analysis[simulated_collection_analysis$metric == metric &
                                        simulated_collection_analysis$pair == pair, ]
      
      # Get all the p-values for the current cell type pair and metric
      p_values3D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "3D", "p_value"]
      metrics_vs_pairs_matrix3D[pair, metric] <- sum(p_values3D[!is.na(p_values3D)] < threshold) / sum(!is.na(p_values3D))
      
      p_values2D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "2D", "p_value"]
      metrics_vs_pairs_matrix2D[pair, metric] <- sum(p_values2D[!is.na(p_values2D)] < threshold) / sum(!is.na(p_values2D))
    }
  }
  
  # Remove AUC from metrics
  metrics <- sub("_AUC$", "", metrics)
  colnames(metrics_vs_pairs_matrix3D) <- metrics
  colnames(metrics_vs_pairs_matrix2D) <- metrics
  
  heatmap3D <- Heatmap(metrics_vs_pairs_matrix3D, 
                       name = "sig",
                       col = colorRamp2(c(0, 0.5, 1),
                                        c("lightgrey", "#b8db50", "#007128")),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       na_col = "white",
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if (!is.nan(metrics_vs_pairs_matrix3D[i, j])) { # Check if the value is not NaN
                           grid.text(sprintf("%.2f", metrics_vs_pairs_matrix3D[i, j]), x, y, gp = gpar(fontsize = 10))
                         }
                       },
                       row_title = "3D",
                       column_labels = metrics)
  
  heatmap2D <- Heatmap(metrics_vs_pairs_matrix2D, 
                       name = "sig",
                       col = colorRamp2(c(0, 0.5, 1),
                                        c("lightgrey", "#b8db50", "#007128")),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       na_col = "white",
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if (!is.nan(metrics_vs_pairs_matrix2D[i, j])) { # Check if the value is not NaN
                           grid.text(sprintf("%.2f", metrics_vs_pairs_matrix2D[i, j]), x, y, gp = gpar(fontsize = 10))
                         }
                       },
                       row_title = "2D",
                       column_labels = metrics)
  
  fig <- heatmap3D %v% heatmap2D
  fig <- draw(fig, gap = unit(1, "cm"))
  
  return(fig)
}

# *** Function 5. ***
# Same as function 4, but taking the difference between 3D and 2D
heatmap_p_values_for_each_metric_and_pair_3D_subtract_2D <- function(simulated_collection_analysis, 
                                                                     metrics, 
                                                                     threshold) {
  
  # Add a combined 'pair' column for reference and target
  simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                              simulated_collection_analysis$target,
                                              sep = "/")
  
  # Get cell type combinations pairs
  pairs <- unique(simulated_collection_analysis$pair)
  
  # Set up metrics vs cell type combinations matrix
  metrics_vs_pairs_matrix3D <- matrix(nrow = length(pairs), ncol = length(metrics))
  rownames(metrics_vs_pairs_matrix3D) <- pairs
  colnames(metrics_vs_pairs_matrix3D) <- metrics
  
  # 2D matrix has same format as 3D matrix
  metrics_vs_pairs_matrix2D <- metrics_vs_pairs_matrix3D
  
  # Iterate through each metric and cell type combination
  for (pair in pairs) {
    for (metric in metrics) {
      
      simulated_collection_analysis_filtered <- 
        simulated_collection_analysis[simulated_collection_analysis$metric == metric &
                                        simulated_collection_analysis$pair == pair, ]
      
      # Get all the p-values for the current cell type pair and metric
      p_values3D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "3D", "p_value"]
      metrics_vs_pairs_matrix3D[pair, metric] <- sum(p_values3D[!is.na(p_values3D)] < threshold) / sum(!is.na(p_values3D))
      
      p_values2D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "2D", "p_value"]
      metrics_vs_pairs_matrix2D[pair, metric] <- sum(p_values2D[!is.na(p_values2D)] < threshold) / sum(!is.na(p_values2D))
    }
  }

  metrics_vs_pairs_matrix_diff <- metrics_vs_pairs_matrix3D - metrics_vs_pairs_matrix2D

  # Remove AUC from metrics
  metrics <- sub("_AUC$", "", metrics)
  colnames(metrics_vs_pairs_matrix_diff) <- metrics
  
  heatmap_diff <- Heatmap(metrics_vs_pairs_matrix_diff, 
                          name = "sig_dff",
                          col = colorRamp2(c(-1, 0, 1),
                                           c("#0062c5", "white", "#bb0036")),
                          rect_gp = gpar(col = "black", lwd = 1),
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          na_col = "white",
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if (!is.nan(metrics_vs_pairs_matrix_diff[i, j])) { # Check if the value is not NaN
                              grid.text(sprintf("%.2f", metrics_vs_pairs_matrix_diff[i, j]), x, y, gp = gpar(fontsize = 10))
                            }
                          },
                          column_labels = metrics,
                          column_names_rot = 90,
                          column_names_centered = FALSE,
                          column_names_gp = gpar(col = "#4D4D4D"),
                          
                          column_title = "Metric", 
                          column_title_side = "bottom")

  return(heatmap_diff)
}


# Plotting setup -- the only part you need to change **** ----- 
setwd("~/R/SD3_data")
simulated_collection_analysis <- readRDS("ME_vs_ME_simulated_collection_analysis.RDS")
group_comparison_name <- "ME_vs_ME"

# Plotting - don't need to change below, just run code below ****  -----

metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "FCO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


# Plot 1. Box plot for p-values
fig <- box_plot_p_values_3D_vs_2D(simulated_collection_analysis)

setwd(paste("~/R/plots/SD3/", group_comparison_name, sep = ""))
pdf(paste(group_comparison_name, "box_plot_p_values_3D_vs_2D.pdf", sep = "_"), width = 10, height = 8)
print(fig)
dev.off()


# Plot 2. Box plot for p-values for each metric and pair
figs <- box_plot_p_values_for_each_metric_and_pair_3D_vs_2D(simulated_collection_analysis,
                                                            metrics)

setwd(paste("~/R/plots/SD3/", group_comparison_name, sep = ""))
pdf(paste(group_comparison_name, "box_plot_p_values_for_each_metric_and_pair_3D_vs_2D.pdf", sep = "_"))
for (metric in metrics) {
  print(figs[[metric]])
}
dev.off()

# Plot 3. Box plot for log of p-values for each metric and pair
figs <- box_plot_log_p_values_for_each_metric_and_pair_3D_vs_2D(simulated_collection_analysis,
                                                                metrics)

setwd(paste("~/R/plots/SD3/", group_comparison_name, sep = ""))
pdf(paste(group_comparison_name, "box_plot_log_p_values_for_each_metric_and_pair_3D_vs_2D.pdf", sep = "_"))
for (metric in metrics) {
  print(figs[[metric]])
}
dev.off()


# Plot 4. Heatmap for p-values for each metric and pair
fig <- heatmap_p_values_for_each_metric_and_pair_3D_vs_2D(simulated_collection_analysis,
                                                          metrics,
                                                          threshold = 0.05)

setwd(paste("~/R/plots/SD3/", group_comparison_name, sep = ""))
pdf(paste(group_comparison_name, "heatmap_p_values_for_each_metric_and_pair_3D_vs_2D.pdf", sep = "_"), width = 10, height = 8)
print(fig)
dev.off()

# Plot 5. Heatmap difference for p-values for each metric and pair
fig <- heatmap_p_values_for_each_metric_and_pair_3D_subtract_2D(simulated_collection_analysis,
                                                                metrics,
                                                                threshold = 0.05)

setwd(paste("~/R/plots/SD3/", group_comparison_name, sep = ""))
pdf(paste(group_comparison_name, "heatmap_p_values_for_each_metric_and_pair_3D_subtract_2D.pdf", sep = "_"), width = 10, height = 4)
print(fig)
dev.off()

