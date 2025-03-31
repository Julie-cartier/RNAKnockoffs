
source("/home/julie/Documents/Paper_codes/New_cohorts_EXPERIMENTS/AEGIS/Settings.R")  # get X_AEGIS_exp and X_CRUKPAP
load("/home/julie/Documents/Paper_codes/Data/BC/X_BC_exp.R") # get X_BC_exp

# Function to create correlation heatmaps for the different data sets

create_heatmap <- function(cor_matrix, dataset_name, save_name){
  
  # Perform hierarchical clustering for better visualization
  hc <- hclust(dist(1 - cor_matrix))  
  
  # re-order features
  ordered_features <- colnames(cor_matrix)[hc$order]
  cor_matrix <- cor_matrix[ordered_features, ordered_features]
  
  # Melt correlation matrix to create heatmaps 
  df_Sigma <- melt(cor_matrix)
  colnames(df_Sigma) <- c("Feature1", "Feature2", "Correlation")

  
  HM <- ggplot(df_Sigma, aes(x = Feature1, y = Feature2, fill = Correlation)) +
    geom_tile() +  
    scale_fill_gradient2(low = "#FFFFFF", mid = "#4D9DDA", high = "#081AD1", midpoint = 0.5) + 
    theme_minimal() +  
    labs(x = dataset_name, y = "", fill = "Pearson Correlation : ") +
    theme(
      axis.text.x = element_blank(),  
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),  
      legend.title = element_text(size = 9, face = "bold") 
    ) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5))
  
  #ggsave(file = str_glue("/home/julie/Documents/Paper_codes/Paper_graphs/Correlation/Figures/",save_name,".png"), width = 18, height = 20, units = "cm")
  return(HM)
}

# Create heatmaps

HM_BC <- create_heatmap(abs(cor(X_BC)), "Breast Cancer Features", "X_BC_exp")

HM_AEGIS <- create_heatmap(abs(cor(X_AEGIS_exp)), "AEGIS Features", "X_AEGIS_exp")

HM_CRUKPAP <- create_heatmap(abs(cor(X_CRUKPAP)), "CRUKPAP Features", "X_CRUKPAP")



# Combine plots 
final_plot <- ggarrange(HM_CRUKPAP, HM_AEGIS, HM_BC, ncol = 3, labels = c("(a) ", "(b)", "(c)"), common.legend = TRUE, font.label = list(size = 8))


# Save the final merged plot
ggsave("/home/julie/Documents/Paper_codes/Paper_graphs/Correlation/Figures/SM_figure_17.png", final_plot, width = 20, height = 7, units = "cm")
