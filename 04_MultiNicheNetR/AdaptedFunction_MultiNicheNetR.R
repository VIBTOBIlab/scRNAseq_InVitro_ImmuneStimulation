#### make_group_lr_prod_plots ####
# OBJECTIVE
#     Create an similar plot as obtained using the make_sample_lr_prod_plot per group instead of per sample
#     This implies calculating the median over samples of the scaled LR pseudobulk expression product
# ADAPTATIONS
#     1. Instead of the scaled_LR_pb_prod that is plotted in the dots (color) for all the samples separately, a median value over these samples (per condition) is calculated and used


make_group_lr_prod_plots <- function(prioritization_tables, prioritized_tbl_oi){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  pb_exprs_data_subset = prioritization_tables$sample_prioritization_tbl %>% 
    dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
    dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                  lr_interaction = paste(ligand, receptor, sep = " - "))  %>%  
    dplyr::arrange(receiver) %>% 
    dplyr::group_by(receiver) %>%  
    dplyr::arrange(sender, .by_group = TRUE)
  
  pb_exprs_data_subset = pb_exprs_data_subset %>% 
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = pb_exprs_data_subset$sender_receiver %>% unique()))%>%
    group_by(group, sender_receiver, lr_interaction)%>%
    mutate(median_scaled_LR_pb_prod = median(scaled_LR_pb_prod))
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4.25)
  names(keep_sender_receiver_values) = levels(pb_exprs_data_subset$keep_sender_receiver)
  
  p1 = pb_exprs_data_subset %>%
    ggplot(aes(group, lr_interaction, color = median_scaled_LR_pb_prod)) +
    geom_point() +
    facet_grid(sender_receiver~., scales = "free", space = "free", switch = "y")+
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(2.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      # strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Median scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + xlab("") + ylab("") +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(pb_exprs_data_subset$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  return(p1)
  
}


#### make_group_lr_prod_activity_plots ####
# OBJECTIVE
#     Create an similar plot as obtained using the make_sample_lr_prod_activity_plots per group instead of per sample
#     This implies calculating the median over samples of the scaled LR pseudobulk expression product
# ADAPTATIONS
#     1. Instead of the scaled_LR_pb_prod that is plotted in the dots (color) for all the samples separately, a median value over these samples (per condition) is calculated and used
#     2. A point graph was added to indicate in which condition the ligand-receptor pair was prioritized


make_group_lr_prod_activity_plots <- function(prioritization_tables, 
                                              prioritized_tbl_oi, widths = NULL,
                                              prioritized_tbl_oi_perCondition = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  sample_data <- prioritization_tables$sample_prioritization_tbl %>% 
    dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
    dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                  lr_interaction = paste(ligand, receptor, sep = " - "))   %>%  
    dplyr::arrange(receiver) %>% 
    dplyr::group_by(receiver) %>%  
    dplyr::arrange(sender, .by_group = TRUE)
  
  sample_data = sample_data %>% 
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sample_data$sender_receiver %>% unique()))%>%
    group_by(group, sender_receiver, lr_interaction)%>%
    mutate(median_scaled_LR_pb_prod = median(scaled_LR_pb_prod))
  
  condition_data <- prioritized_tbl_oi_perCondition%>%
    # Perform similar preprocessing as the other tables
    mutate(lr_interaction = paste(ligand, receptor, sep = " - "),
           sender_receiver = paste(sender, receiver, sep = " --> "))%>%
    arrange(receiver)%>%
    group_by(receiver)%>%
    arrange(sender, .by_group = T)%>%
    # Factorization needed to ensure the same order used in the other graphs
    mutate(sender_receiver = factor(sender_receiver, levels = sender_receiver%>%unique))
  
  group_data = prioritization_tables$group_prioritization_table_source  %>% 
    dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                  lr_interaction = paste(ligand, receptor, sep = " - "))  %>% 
    dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, activity, activity_scaled, direction_regulation, prioritization_score) %>% 
    dplyr::filter(id %in% sample_data$id) %>%  
    dplyr::arrange(receiver) %>% 
    dplyr::group_by(receiver) %>%  
    dplyr::arrange(sender, .by_group = TRUE)
  group_data = group_data %>% 
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% 
                                             unique()))
  
  group_data_celltype_specificity = prioritization_tables$group_prioritization_tbl  %>% 
    dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))%>% 
    dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% 
    dplyr::filter(id %in% sample_data$id) %>%  
    dplyr::arrange(receiver) %>% 
    dplyr::group_by(receiver) %>%  
    dplyr::arrange(sender, .by_group = TRUE)
  group_data_celltype_specificity = group_data_celltype_specificity %>% 
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data_celltype_specificity$sender_receiver %>% 
                                             unique()))
  
  group_data_frac_expression = prioritization_tables$group_prioritization_table_source  %>% 
    dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% 
    dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_frac_expression = group_data_frac_expression %>% 
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% 
                                             unique()))
  
  group_data = group_data %>% 
    inner_join(group_data_celltype_specificity) %>% 
    inner_join(group_data_frac_expression)
  group_data = group_data %>% 
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% 
                                             unique()))
  rm(group_data_celltype_specificity)
  rm(group_data_frac_expression)
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
  
  p1 = sample_data %>%
    ggplot(aes(group, lr_interaction, color = median_scaled_LR_pb_prod)) +
    geom_point() +
    facet_grid(sender_receiver~., scales = "free", space = "free", switch = "y")+
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Median scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + 
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  p_conditions <- ggplot(condition_data,
                         mapping = aes(x = group, y = lr_interaction))+
    geom_point(shape = 4)+
    facet_grid(sender_receiver~., scales = "free", space = "free", switch = "y")+
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      legend.position = 'none'
    )
    
  
  p2 = group_data %>%
    ggplot(aes(direction_regulation , lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 90),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold", angle = 90),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.placement = "outside"
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver")
  max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
  
  p3 = group_data %>%
    ggplot(aes(direction_regulation , lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 90),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.placement = "outside"
    ) + labs(fill = "Ligand\nActivity in Receiver")
  max_activity = (group_data$activity) %>% max()
  min_activity = (group_data$activity) %>% min()
  custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white",high = "darkorange",midpoint = 0)
  
  p3 = p3 + custom_scale_fill
  
  # add the plot visualizing cell-type specificity
  # cs_data = group_data %>% filter(group %in% prioritized_tbl_oi$group) %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data = group_data %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% tidyr::gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data$LR[cs_data$LR == "scaled_pb_ligand"] = "ligand"
  cs_data$LR[cs_data$LR == "scaled_pb_receptor"] = "receptor"
  frac_data = group_data %>% distinct(sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% tidyr::gather(LR, fraction_expression, fraction_ligand_group:fraction_receptor_group)
  frac_data$LR[frac_data$LR == "fraction_ligand_group"] = "ligand"
  frac_data$LR[frac_data$LR == "fraction_receptor_group"] = "receptor"
  
  cs_data = cs_data %>% inner_join(frac_data)
  
  p_cs = cs_data %>% 
    ggplot(aes(LR , lr_interaction, color = celltype_specificity, size = fraction_expression)) +
    geom_point() +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    viridis::scale_color_viridis() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 90),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.placement = "outside"
    ) + labs(color = "Scaled celltype specificity") + labs(size = "Fraction of expression")
  
  
  if(!is.null(widths)){
    p = patchwork::wrap_plots(
      p1, p_conditions, p2, p3, p_cs,
      nrow = 1,guides = "collect",
      widths = widths
    )
  } else {
    p = patchwork::wrap_plots(
      p1, p_conditions, p2, p3, p_cs,
      nrow = 1,guides = "collect",
      widths = c(sample_data$group %>% unique() %>% length(), 
                ( condition_data$group %>% unique() %>% length())/2, 
                 2*(sample_data$group %>% unique() %>% length()), 
                 2*(sample_data$group %>% unique() %>% length()),
                 2*(sample_data$group %>% unique() %>% length()))
    )
  }
  
  return(p)
  
}





