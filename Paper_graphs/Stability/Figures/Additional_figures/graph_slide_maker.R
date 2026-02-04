y.rep <- 10
beta.tot <- c()
table.selected.genes.LPLR.tot <- data.frame()
table.selected.genes.LPLR.oracle.tot <- data.frame()

data_analysis <- data.frame()

for (rep in 1:y.rep){
  print(rep)
  # load selected features for all methods 
  
  
  load(str_glue("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table.selected_genes_cv_KOPI_",rep,".R"))
  load(str_glue("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table_selected_genes_LPLR_",rep,".R"))
  load(str_glue("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/y_rep/table_selected_genes_LPLR_oracle_",rep,".R"))
  load(str_glue("/home/julie/Documents/Paper_codes/CRUKPAP_EXPERIMENTS/Stability/R_files/y_rep/SS_selected_features_rep_", rep, ".Rdata"))
  
  
  # Feature selection metrics (used in the paper table)
  
  list_non_null = which(beta!=0)
  true_genes = rep("False discoveries ", ncol(X_real))
  true_genes[ which(beta!=0)] <- "True genes "
  bin_width <- 0.01
  
  # prop_selection.genes_LPLR <- rowMeans(table.selected.genes.LPLR)
  # prop_selection.genes_LPLR_0 <- prop_selection.genes_LPLR[ !prop_selection.genes_LPLR == 0]
  # true_genes_LPLR_0 = true_genes[ !prop_selection.genes_LPLR == 0]
  # 
  # df.genes <- data.frame(freq = prop_selection.genes_LPLR_0, status = as.factor(true_genes_LPLR_0))
  # 
  # prop_selection.genes_LPLR_oracle <- rowMeans(table.selected.genes.LPLR.oracle)
  # prop_selection.genes_LPLR_oracle_0 <- prop_selection.genes_LPLR_oracle[ !prop_selection.genes_LPLR_oracle == 0]
  # true_genes_LPLR_oracle_0 <- true_genes[ !prop_selection.genes_LPLR_oracle == 0]
  # 
  # df.genes.oracle <- data.frame(freq = prop_selection.genes_LPLR_oracle_0, status = as.factor(true_genes_LPLR_oracle_0))
  # 
  # 
  # prop_selection_0.2 <- rowSums(table.selected_genes.KOPI_0.2)
  # prop_selection_0.3 <- rowSums(table.selected_genes.KOPI_0.3)
  # 
  # prop_selection_0.2_0 <- prop_selection_0.2[ !prop_selection_0.2 == 0]
  # prop_selection_0.3_0 <- prop_selection_0.3[ !prop_selection_0.3 == 0]
  # 
  # data_0.2_KOPI <- data.frame(freq = prop_selection_0.2_0/100, status = as.factor(true_genes[ !prop_selection_0.2 == 0]))
  # data_0.3_KOPI <- data.frame(freq = prop_selection_0.3_0/100, status = as.factor(true_genes[ !prop_selection_0.3 == 0]))
  
  
  
  
  # get the non null beta
  
  list_non_null <- which(beta!=0)
  true_genes <- rep("False discoveries ", ncol(X_real))
  true_genes[which(beta!=0)] <- "True genes "
  bin_width <- 0.01
  
  # get selection frequency 
  
  # for the LASSO
  table.selected.genes.LPLR
  prop_selection.genes_LPLR <- rowMeans(table.selected.genes.LPLR)
  print(length(which(colMeans(table.selected.genes.LPLR) == 0)))
  prop_selection.genes_LPLR_0 <- prop_selection.genes_LPLR[ !prop_selection.genes_LPLR == 0]# Remove features that have never been selected
  true_genes_LPLR_0 <- true_genes[!prop_selection.genes_LPLR == 0]
  
  df.genes <- data.frame(freq = prop_selection.genes_LPLR_0, status = as.factor(true_genes_LPLR_0))
  
  max_bin_count_LPLR <- df.genes$freq %>% cut(breaks = seq(0, 1, by = bin_width), include.lowest = TRUE) %>% table() %>% max() # To adjust y-axis size
  
  # for the oracle LASSO
  
  prop_selection.genes_LPLR_oracle <- rowMeans(selected_genes.MB)
  print(length(which(colMeans(table.selected.genes.LPLR.oracle) == 0)))
  prop_selection.genes_LPLR_oracle_0 <- prop_selection.genes_LPLR_oracle[ !prop_selection.genes_LPLR_oracle == 0]# Remove features that have never been selected
  true_genes_LPLR_oracle_0 <- true_genes[ !prop_selection.genes_LPLR_oracle == 0]
  
  df.genes.oracle <- data.frame(freq = prop_selection.genes_LPLR_oracle_0, status = as.factor(true_genes_LPLR_oracle_0))
  
  max_bin_count_LPLR_oracle <- df.genes.oracle$freq %>% cut(breaks = seq(0, 1, by = bin_width), include.lowest = TRUE) %>% table() %>% max() # To adjust y-axis size
  
  # For KOPI
  
  prop_selection_0.2 <- rowMeans(table.selected_genes.KOPI_0.2)
  prop_selection_0.3 = rowMeans(table.selected_genes.KOPI_0.3)
  print(length(which(colMeans(table.selected_genes.KOPI_0.2) == 0)))
  print(length(which(colMeans(table.selected_genes.KOPI_0.3) == 0)))
  
  prop_selection_0.2_0 <- prop_selection_0.2[ !prop_selection_0.2 == 0]# Remove features that have never been selected
  prop_selection_0.3_0 <- prop_selection_0.3[ !prop_selection_0.3 == 0]# Remove features that have never been selected
  
  data_0.2_KOPI <- data.frame(freq = prop_selection_0.2_0, status = as.factor(true_genes[ !prop_selection_0.2 == 0]))
  data_0.3_KOPI <- data.frame(freq = prop_selection_0.3_0, status = as.factor(true_genes[ !prop_selection_0.3 == 0]))
  
  max_bin_count_KOPI_0.2 <- data_0.2_KOPI$freq %>% cut(breaks = seq(0, 1, by = bin_width), include.lowest = TRUE) %>% table() %>% max() # To adjust y-axis size
  max_bin_count_KOPI_0.3 <- data_0.3_KOPI$freq %>% cut(breaks = seq(0, 1, by = bin_width), include.lowest = TRUE) %>% table() %>% max()# To adjust y-axis size
  
  
  
  ylim.value <- max(max_bin_count_KOPI_0.3, max_bin_count_KOPI_0.2)
  
  hist.genes <- ggplot(df.genes, aes(x = freq, fill = status)) +
    geom_histogram(binwidth = 0.01, color = "black") +
    coord_cartesian(xlim = c(0, 1)) +
    #ylim(0, ylim.value) +
    labs(title = "Selection proportion of genes selected more than once",
         x = "Selection frequency",
         y = "Number of genes") +
    theme(axis.text.x = element_text(size = 26, face = "bold"),
          axis.title.x = element_text(size = 28),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 30, face = "bold"),
          plot.subtitle = element_text(hjust = 0),
          plot.caption = element_text(size = 22, face = "bold"),
          panel.background = element_rect(fill = "#F6F6F6"),
          panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          plot.title = element_blank(),
          legend.position = c(0.3, 0.8),
          legend.text = element_text(size = 28),    
          legend.title = element_blank(),  
          plot.margin = unit(c(0.6, 0.2, 0.2, 0.6), "cm")
    ) +
    scale_fill_manual(
      values = c("#7c6fb0", "#FF5733", "#FFFF66")
      #name = expression("LASSO (" * lambda[min] * "), n = 2150")
    )
  
  
  
  
  hist.genes.oracle <- ggplot(df.genes.oracle, aes(x = freq, fill = status)) +
    geom_histogram(binwidth = 0.01, color = "black") +
    #ylim(0, ylim.value) +
    coord_cartesian(xlim = c(0, 1)) + 
    labs(title = "Selection proportion of genes selected more than once",
         x = "Selection frequency",
         y = "Number of genes") +
    theme(axis.text.x = element_text(size = 26, face = "bold"),
          axis.title.x = element_text(size = 28),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 30, face = "bold"),
          plot.subtitle = element_text(hjust = 0),
          plot.caption = element_text(size = 22, face = "bold"),
          panel.background = element_rect(fill = "#F6F6F6"),
          panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          plot.title = element_blank(),
          legend.position = c(0.3, 0.8),
          legend.text = element_text(size = 28),    
          legend.title = element_blank(),  
          plot.margin = unit(c(0.6, 0.2, 0.2, 0.6), "cm")
    ) +
    scale_fill_manual(
      values = c("#b2a3e8", "#FF5733", "#FFFF66")
      #name = expression("LASSO (" * lambda[oracle] * "), n = 381")
    )
  
  
  
  hist_KOPI_0.2 <- ggplot(data_0.2_KOPI, aes(x = freq, fill = status)) +
    geom_histogram(binwidth = 0.01, color = "black" ) +
    coord_cartesian(xlim = c(0, 1)) + 
    ylim(0, ylim.value) +
    labs(title = "Selection proportion of genes selected more than once",
         x = "Selection frequency",
         y = "Number of genes") +
    theme(axis.text.x = element_text(size = 26, face = "bold"),
          axis.title.x = element_text(size = 28),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 30, face = "bold"),
          plot.subtitle = element_text(hjust = 0),
          plot.caption = element_text(size = 22, face = "bold"),
          panel.background = element_rect(fill = "#F6F6F6"),
          panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          plot.title = element_blank(),
          legend.position = c(0.3, 0.8),
          legend.text = element_text(size = 28),    
          legend.title = element_blank(),  
          plot.margin = unit(c(0.6, 0.2, 0.2, 0.6), "cm")
    ) +
    scale_fill_manual(values = c("#17becf", "#FF5733", "#FFFF66"))
  #name = "KOPI (q = 0.2), n = 113") # Define colors for different statuses
  
  hist_KOPI_0.2
  
  
  
  
  hist_KOPI_0.3 <- ggplot(data_0.3_KOPI, aes(x = freq, fill = status)) +
    geom_histogram(binwidth = 0.01, color = "black" ) +
    coord_cartesian(xlim = c(0, 1)) + 
    ylim(0, ylim.value) +
    labs(title = "Selection proportion of genes selected more than once",
         x = "Selection frequency",
         y = "Number of genes") +
    theme(axis.text.x = element_text(size = 26, face = "bold"),
          axis.title.x = element_text(size = 28),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 30, face = "bold"),
          plot.subtitle = element_text(hjust = 0),
          plot.caption = element_text(size = 24, face = "bold"),
          panel.background = element_rect(fill = "#F6F6F6"),
          panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
          plot.title = element_blank(),
          legend.position = c(0.3, 0.8),
          legend.text = element_text(size = 28),    
          legend.title = element_blank(),
          plot.margin = unit(c(0.6, 0.2, 0.2, 0.6), "cm")
    ) +
    scale_fill_manual(values =  c("#008899", "#FF5733", "#FFFF66"))
  #name = "KOPI (q = 0.3), n = 520") # Define colors for different status
  
  hist_KOPI_0.3
  
  
  
  
  
  
  graph_min <- ggarrange(hist.genes, hist_KOPI_0.2, hist_KOPI_0.3, ncol = 3, common.legend = FALSE, font.label = list(size = 22)) 
  
  ggsave(file = str_glue("/home/julie/Documents/Paper_codes/Paper_graphs/Stability/Figures/q3/min_y_",rep,".png"), plot = graph_min, height = 9, width = 22)
  
  graph_oracle <- ggarrange(hist.genes.oracle, hist_KOPI_0.2, hist_KOPI_0.3, ncol = 3, common.legend = FALSE, font.label = list(size = 22)) 
  
  ggsave(file = str_glue("/home/julie/Documents/Paper_codes/Paper_graphs/Stability/Figures/q3/SS_y_",rep,".png"), plot = graph_oracle, height = 9, width = 22)
}