# Plot predictions for a given patient and different data splits (C values)
plot_states <- function(r, title_) {
  colsr <- colnames(r)
  colnames(r)[1:length(colsr)-1] <- strtrim(colsr[1:length(colsr)-1], 5)
  
  r[r == 9] <- "None"
  r[r == 0] <- "Euthymia"
  r[r == 1] <- "Depression"
  r[r == 2] <- "Mania"
  r[r == 3] <- "Mixed"
  
  
  df_long <- r %>%
    mutate(index = row_number()) %>%
    pivot_longer(cols = -index)
  
  cols <- c('white', 'chartreuse4', 'darkorchid3', 'red', 'darkorange')
  names(cols) <- c("None", "Euthymia", "Depression", "Mania", "Mixed")
  q <- ggplot(df_long, aes(x = index, y = name, fill = factor(value))) +
    geom_tile(height = .3) +
    scale_fill_manual(values = cols, limits = force, breaks = c("Euthymia", 
                                                                "Depression", 
                                                                "Mania", 
                                                                "Mixed"), name = "State") +
    
    
    scale_x_continuous(limits = c(0, dim(r)[1]),
                       breaks = seq(from = 0, 
                                    to = dim(r)[1], 
                                    by = floor(dim(r)[1]/10)), 
                       guide = guide_axis(check.overlap = TRUE),
                       expand = c(0,0))+
    labs(title = title_, x = "Observation index", y = "C") +
    theme(legend.position = "top", 
          legend.title = element_blank(),
          legend.key.height =  unit(0.2, "cm"),
          legend.key.width =  unit(3, "cm"),
          legend.key.spacing.x = unit(1, "cm"),
          text = element_text(size = 25)) +
    guides(fill  = guide_legend(label.position = "top", 
                                nrow = 1, reverse = F))
  return(q)
}

# Visualize EMR and CM for different data splits (1-C)
visualize_performance <- function(dat) {
  
  vcol <- viridis(4)
  
  measures <- c("solid", "longdash") ; names(measures) <- c("EMR", "CM")
  models <- c(vcol[1], vcol[2], vcol[3], vcol[4]) ; names(models) <- c("HMM", "HMMe", "RF", "GNB")
  methods <- c(3,2,1) ; names(methods) <- c("CCA", "PCA", "mRMR")
  
  ggplot(dat, aes(x = TestSize/100, y = value, 
                  color = Model, 
                  linetype = Measure,
                  shape = Method)) + 
    geom_line() +
    theme_bw() +
    scale_linetype_manual(values = measures) +
    scale_shape_manual(values = methods) + 
    scale_color_manual( values = models) +
    labs(x = "1 - C", y = "Value", color  = "Model", linetype = "Measure", size = "Model") +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "top",
          text = element_text(size = 25),
          legend.key.height =  unit(0.5, "cm"),
          legend.key.width =  unit(3, "cm"),
          legend.key.spacing.x = unit(1, "cm"),
          legend.title = element_blank()) +
    guides(linetype = guide_legend(label.position = "top", 
                                   nrow = 1, reverse = F),
           color  = guide_legend(label.position = "top", 
                                 nrow = 1, reverse = F),
           size = guide_legend(label.position = "top", 
                               nrow = 1, reverse = F))
  
}

# Visualize the (in)consistency and stability of the prediction 
visualize_pi_of_k <- function(dat) {
  dat_melt <- cbind(melt(dat$mean), 
                    "lower" = melt(dat$lower)$value, 
                    "upper" = melt(dat$upper)$value,
                    "K" = dat$K)
  
  Ks <- dat_melt$K
  dat_melt <- data.frame(separate_wider_delim(dat_melt, cols = variable, 
                                              delim = "_", names = c("Model", "Method")))
  
  vcol <- viridis(4)
  models <- c(vcol[1], vcol[2], vcol[3], vcol[4]) ; names(models) <- c("HMM", "HMMe", "RF", "GNB")
  methods <- c(3,2,1) ; names(methods) <- c("CCA", "PCA", "mRMR")
  
  ggplot(dat_melt, 
         aes(x = K, y = value, colour = Model, shape = Method)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = lower,
                    ymax = upper, 
                    fill = Model),
                alpha = 0.15, linetype = 0) +
    theme_bw() +
    scale_shape_manual(values = methods) + 
    scale_color_manual( values = models) +
    labs(x = "K", y = "Performance Index") +
    scale_x_continuous(limits = c(min(Ks), max(Ks)),
                       breaks = seq(from = min(Ks), 
                                    to = max(Ks), 
                                    by = 0.5)) +
    scale_y_continuous(limits = c(0, 1))  + 
    theme(legend.position = "top", 
          plot.title = element_text(hjust = 0.5), 
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
          text = element_text(size = 25),
          legend.key.height =  unit(0.5, "cm"),
          legend.key.width =  unit(3, "cm"),
          legend.key.spacing.x = unit(1, "cm")) +
    guides(colour  = guide_legend(label.position = "top", 
                                  nrow = 1, reverse = F),
           fill = guide_legend(label.position = "top", 
                               nrow = 1, reverse = F))
  
}
