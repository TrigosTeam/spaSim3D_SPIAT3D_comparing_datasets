# Code to plot simulated collection analysis, showing only cell pair A/B (ignoring cell pairs A/A, B/A, B/B)

# Libraries ----
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


# Functions -----
# *** Function 4. ***
# Create a heatmap with metrics on the rows and cell type combination on the columns.
# p-value false positive or false negative rate shown
heatmap_p_values_for_each_metric_and_for_pair_A_B_3D_vs_2D <- function(simulated_collection_analysis, 
                                                                       metrics, 
                                                                       threshold) {
  
  # Add a combined 'pair' column for reference and target
  simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                              simulated_collection_analysis$target,
                                              sep = "|")
  
  simulated_collection_analysis <- simulated_collection_analysis[simulated_collection_analysis$pair == "A|B", ]
  
  # Set up metrics vs cell type combinations matrix
  metrics_vs_pairs_matrix3D <- matrix(nrow = 1, ncol = length(metrics))
  colnames(metrics_vs_pairs_matrix3D) <- metrics
  
  # 2D matrix has same format as 3D matrix
  metrics_vs_pairs_matrix2D <- metrics_vs_pairs_matrix3D
  
  # Iterate through each metric and cell type combination
  for (metric in metrics) {
    
    simulated_collection_analysis_filtered <- 
      simulated_collection_analysis[simulated_collection_analysis$metric == metric, ]
    
    # Get all the p-values for the current  metric
    p_values3D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "3D", "p_value"]
    metrics_vs_pairs_matrix3D[1, metric] <- sum(p_values3D[!is.na(p_values3D)] < threshold) / sum(!is.na(p_values3D))
    
    p_values2D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "2D", "p_value"]
    metrics_vs_pairs_matrix2D[1, metric] <- sum(p_values2D[!is.na(p_values2D)] < threshold) / sum(!is.na(p_values2D))
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
                       column_labels = metrics,
                       show_row_names = FALSE)
  
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
                       column_labels = metrics,
                       show_row_names = FALSE)
  
  fig <- heatmap3D %v% heatmap2D
  fig <- draw(fig, gap = unit(1, "cm"))
  
  return(fig)
}

# *** Function 5. ***
# Same as function 4, but taking the difference between 3D and 2D
heatmap_p_values_for_each_metric_and_for_pair_A_B_3D_subtract_2D <- function(simulated_collection_analysis, 
                                                                             metrics, 
                                                                             threshold) {
  
  # Add a combined 'pair' column for reference and target
  simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                              simulated_collection_analysis$target,
                                              sep = "|")
  
  simulated_collection_analysis <- simulated_collection_analysis[simulated_collection_analysis$pair == "A|B", ]
  
  # Set up metrics vs cell type combinations matrix
  metrics_vs_pairs_matrix3D <- matrix(nrow = 1, ncol = length(metrics))
  colnames(metrics_vs_pairs_matrix3D) <- metrics
  
  # 2D matrix has same format as 3D matrix
  metrics_vs_pairs_matrix2D <- metrics_vs_pairs_matrix3D
  
  # Iterate through each metric and cell type combination
  for (metric in metrics) {
    
    simulated_collection_analysis_filtered <- 
      simulated_collection_analysis[simulated_collection_analysis$metric == metric, ]
    
    # Get all the p-values for the current  metric
    p_values3D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "3D", "p_value"]
    metrics_vs_pairs_matrix3D[1, metric] <- sum(p_values3D[!is.na(p_values3D)] < threshold) / sum(!is.na(p_values3D))
    
    p_values2D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "2D", "p_value"]
    metrics_vs_pairs_matrix2D[1, metric] <- sum(p_values2D[!is.na(p_values2D)] < threshold) / sum(!is.na(p_values2D))
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
                          show_row_names = FALSE,
                          
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
group_comparison_name <- "ME_vs_ME_pair_A_B"

# Plotting - don't need to change below, just run code below ****  -----

# Define constants
metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "FCO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


# Plot 4. Heatmap for p-values for each metric and pair
fig <- heatmap_p_values_for_each_metric_and_for_pair_A_B_3D_vs_2D(simulated_collection_analysis,
                                                                  metrics,
                                                                  threshold = 0.05)

setwd("~/R/plots/SD3")
pdf(paste(group_comparison_name, "heatmap_p_values_for_each_metric_and_for_pair_A_B_3D_vs_2D.pdf", sep = "_"), width = 10, height = 3)
print(fig)
dev.off()

# Plot 5. Heatmap difference for p-values for each metric and pair
fig <- heatmap_p_values_for_each_metric_and_for_pair_A_B_3D_subtract_2D(simulated_collection_analysis,
                                                                        metrics,
                                                                        threshold = 0.05)

setwd("~/R/plots/SD3")
pdf(paste(group_comparison_name, "heatmap_p_values_for_each_metric_and_for_pair_A_B_3D_subtract_2D.pdf", sep = "_"), width = 10, height = 2)
print(fig)
dev.off()

# Alternative plotting - combining all collections? -----
setwd("~/R/SD3_data")
ME_vs_ME_simulated_collection_analysis <- readRDS("ME_vs_ME_simulated_collection_analysis.RDS")
MN_vs_MN_simulated_collection_analysis <- readRDS("MN_vs_MN_simulated_collection_analysis.RDS")
ME_vs_lME_simulated_collection_analysis <- readRDS("ME_vs_lME_simulated_collection_analysis.RDS")
MN_vs_lMN_simulated_collection_analysis <- readRDS("MN_vs_lMN_simulated_collection_analysis.RDS")


simulated_collection_analyses <- list("ME_vs_ME" = ME_vs_ME_simulated_collection_analysis,
                                      "MN_vs_MN" = MN_vs_MN_simulated_collection_analysis,
                                      "ME_vs_lME" = ME_vs_lME_simulated_collection_analysis,
                                      "MN_vs_lMN" = MN_vs_lMN_simulated_collection_analysis)

# Create a heatmap with metrics on the columns and collection on the rows
# p-value false positive or false negative rate shown
heatmap_p_values_for_each_metric_and_for_pair_A_B_for_many_collections_3D_vs_2D <- function(simulated_collection_analyses, 
                                                                                            metrics, 
                                                                                            threshold) {
  
  
  # Setup the matrices
  metrics_vs_collections_matrix3D <- matrix(nrow = length(simulated_collection_analyses), ncol = length(metrics))
  rownames(metrics_vs_collections_matrix3D) <- names(simulated_collection_analyses)
  colnames(metrics_vs_collections_matrix3D) <- metrics
  metrics_vs_collections_matrix2D <- metrics_vs_collections_matrix3D
  
  for (i in seq_len(length(simulated_collection_analyses))) {
    simulated_collection_analysis <- simulated_collection_analyses[[i]]
    
    # Add a combined 'pair' column for reference and target
    simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                                simulated_collection_analysis$target,
                                                sep = "|")
    
    simulated_collection_analysis <- simulated_collection_analysis[simulated_collection_analysis$pair == "A|B", ]
    
    # Iterate through each metric and cell type combination
    for (metric in metrics) {
      
      simulated_collection_analysis_filtered <- 
        simulated_collection_analysis[simulated_collection_analysis$metric == metric, ]
      
      # Get all the p-values for the current  metric
      p_values3D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "3D", "p_value"]
      metrics_vs_collections_matrix3D[i, metric] <- sum(p_values3D[!is.na(p_values3D)] < threshold) / sum(!is.na(p_values3D))
      
      p_values2D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "2D", "p_value"]
      metrics_vs_collections_matrix2D[i, metric] <- sum(p_values2D[!is.na(p_values2D)] < threshold) / sum(!is.na(p_values2D))
    } 
  }
  
  # Remove AUC from metrics
  metrics <- sub("_AUC$", "", metrics)
  colnames(metrics_vs_collections_matrix3D) <- metrics
  colnames(metrics_vs_collections_matrix2D) <- metrics
  
  heatmap3D <- Heatmap(metrics_vs_collections_matrix3D, 
                       name = "sig",
                       col = colorRamp2(c(0, 0.5, 1),
                                        c("lightgrey", "#b8db50", "#007128")),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       na_col = "white",
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if (!is.nan(metrics_vs_collections_matrix3D[i, j])) { # Check if the value is not NaN
                           grid.text(sprintf("%.2f", metrics_vs_collections_matrix3D[i, j]), x, y, gp = gpar(fontsize = 10))
                         }
                       },
                       row_title = "3D",
                       column_names_rot = 90,
                       column_names_centered = TRUE)
  
  heatmap2D <- Heatmap(metrics_vs_collections_matrix2D, 
                       name = "sig",
                       col = colorRamp2(c(0, 0.5, 1),
                                        c("lightgrey", "#b8db50", "#007128")),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       na_col = "white",
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if (!is.nan(metrics_vs_collections_matrix2D[i, j])) { # Check if the value is not NaN
                           grid.text(sprintf("%.2f", metrics_vs_collections_matrix2D[i, j]), x, y, gp = gpar(fontsize = 10))
                         }
                       },
                       row_title = "2D",
                       column_names_rot = 90,
                       column_names_centered = TRUE)
  
  fig <- heatmap3D %v% heatmap2D
  fig <- draw(fig, gap = unit(1, "cm"))
  
  return(fig)
}

# Create a heatmap with metrics on the columns and collection on the rows
# p-value false positive or false negative rate shown
heatmap_p_values_for_each_metric_and_for_pair_A_B_for_many_collections_3D_subtract_2D <- function(simulated_collection_analyses, 
                                                                                                  metrics, 
                                                                                                  threshold) {
  
  
  # Setup the matrices
  metrics_vs_collections_matrix3D <- matrix(nrow = length(simulated_collection_analyses), ncol = length(metrics))
  rownames(metrics_vs_collections_matrix3D) <- names(simulated_collection_analyses)
  colnames(metrics_vs_collections_matrix3D) <- metrics
  metrics_vs_collections_matrix2D <- metrics_vs_collections_matrix3D
  
  for (i in seq_len(length(simulated_collection_analyses))) {
    simulated_collection_analysis <- simulated_collection_analyses[[i]]
    
    # Add a combined 'pair' column for reference and target
    simulated_collection_analysis$pair <- paste(simulated_collection_analysis$reference, 
                                                simulated_collection_analysis$target,
                                                sep = "|")
    
    simulated_collection_analysis <- simulated_collection_analysis[simulated_collection_analysis$pair == "A|B", ]
    
    # Iterate through each metric and cell type combination
    for (metric in metrics) {
      
      simulated_collection_analysis_filtered <- 
        simulated_collection_analysis[simulated_collection_analysis$metric == metric, ]
      
      # Get all the p-values for the current  metric
      p_values3D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "3D", "p_value"]
      metrics_vs_collections_matrix3D[i, metric] <- sum(p_values3D[!is.na(p_values3D)] < threshold) / sum(!is.na(p_values3D))
      
      p_values2D <- simulated_collection_analysis_filtered[simulated_collection_analysis_filtered$dimension == "2D", "p_value"]
      metrics_vs_collections_matrix2D[i, metric] <- sum(p_values2D[!is.na(p_values2D)] < threshold) / sum(!is.na(p_values2D))
    } 
  }
  
  metrics_vs_collections_matrix <- metrics_vs_collections_matrix3D - metrics_vs_collections_matrix2D
  
  # Remove AUC from metrics
  metrics <- sub("_AUC$", "", metrics)
  colnames(metrics_vs_collections_matrix) <- metrics
  
  heatmap_diff <- Heatmap(
    metrics_vs_collections_matrix, 
    name = "sig_dff",
    col = colorRamp2(c(-1, 0, 1),
                     c("#0062c5", "white", "#bb0036")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    na_col = "white",
    rect_gp = gpar(col = "black", lwd = 1),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.nan(metrics_vs_collections_matrix[i, j])) {
        grid.text(sprintf("%.2f", metrics_vs_collections_matrix[i, j]),
                  x, y, gp = gpar(fontsize = 12))
      }
    },
      
    show_row_names = F,
    # row_names_rot = 90,
    # row_names_centered = TRUE,
    # row_names_side = "left",
    # row_names_gp = gpar(fontsize = 16),
    
    column_names_rot = 90,
    column_names_centered = FALSE,
    column_names_gp = gpar(fontsize = 16, col = "#4D4D4D"),
    
    column_title = "Metric", 
    column_title_gp = gpar(fontsize = 16),
    column_title_side = "bottom",
    
    show_heatmap_legend = FALSE
    
  )
  
  
  return(heatmap_diff)
}

# Define constants
metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


fig <- heatmap_p_values_for_each_metric_and_for_pair_A_B_for_many_collections_3D_vs_2D(simulated_collection_analyses,
                                                                                       metrics,
                                                                                       0.05)
setwd("~/R/plots/SD3")
pdf("heatmap_p_values_for_each_metric_and_for_pair_A_B_for_many_collections_3D_vs_2D.pdf", width = 9, height = 8.5)
print(fig)
dev.off()




fig <- heatmap_p_values_for_each_metric_and_for_pair_A_B_for_many_collections_3D_subtract_2D(simulated_collection_analyses,
                                                                                             metrics,
                                                                                             0.05)
setwd("~/R/plots/SD3")
pdf("heatmap_p_values_for_each_metric_and_for_pair_A_B_for_many_collections_3D_subtract_2D.pdf", width = 6.5, height = 7.5)
print(fig)
dev.off()



