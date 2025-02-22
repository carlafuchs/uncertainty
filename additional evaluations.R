x1 <- c(-2, 0, 0)
x2 <- c(-2, -2, 0)

variation_v <- c(1, 1.5, 2.5)

list_results_new <- list()

for (i in 2:3) {

  list_results_new[[i]] <- list()

  for(j in 1:3) {
    list_results_new[[i]][[j]] <- data.frame(index = seq.int(from = 20, to = 500, by = 10), u_a = NA, u_e = NA, p_pos = NA, p_neg = NA)
    for (k in seq.int(from = 20, to = 500, by = 10)) {
      df <- generate_data(n = k, separation = 2, noise = 0.1, variation = variation_v[j], ood_fraction = 0)
      list_results_new[[i]][[j]][list_results_new[[i]][[j]]$index == k, c(2:5)] <- evaluate_a(df, x1[i], x2[i])
    }
  }

}

df <- generate_data(n = 100, separation = 2, noise = 0.1, variation = 2.5, ood_fraction = 0)
ggplot(df, aes(x = X1, y = X2, color = factor(y))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "black"), na.value = "gray") +
  labs(title = "Synthetic Data for Logistic Regression",
       color = "Class")  +
  ggtitle(NULL) +
  theme_minimal()

plots_list_al_s <- list()

for (i in 1:3) {
  al_data_1_s <- cbind(list_results_new[[i]][[1]], list_results_new[[i]][[2]], list_results_new[[i]][[3]])
  colnames(al_data_1_s) <- c("step1", "a1", "e1", "pp1", "pn1", "step2", "a2", "e2", "pp2", "pn2", "step3", "a3", "e3", "pp3", "pn3")
  
  al_data_1_s <- al_data_1_s %>%
    select(step1, a1, a2, a3) %>%
    rename("1" = a1, "1.5" = a2, "2.5" = a3) %>%
    pivot_longer(cols = c("1", "1.5", "2.5"))
  
  plots_list_al_s[[i]] <- ggplot(al_data_1_s, aes(x = step1, y = value, color = name)) +
    geom_line(size = 0.6) +  # Thinner lines
    geom_point(size = 1.8, alpha = 0.8) +  # Adjust point size & transparency
    scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +  # Better colors
    labs(title = "Aleatoric Uncertainty",
         x = "Instances",
         y = "Value",
         color = "Standard Deviation =") +
    theme_minimal() +
    ylim(0,1)
}

combined_plot_al_s <- (plots_list_al_s[[1]] | plots_list_al_s[[2]] | plots_list_al_s[[3]]) + 
  plot_layout(guides = "collect")

plots_list_ep_s <- list()

for (i in 1:3) {
  ep_data_1_s <- cbind(list_results_new[[i]][[1]], list_results_new[[i]][[2]], list_results_new[[i]][[3]])
  colnames(ep_data_1_s) <- c("step1", "a1", "e1", "pp1", "pn1", "step2", "a2", "e2", "pp2", "pn2", "step3", "a3", "e3", "pp3", "pn3")
  
  ep_data_1_s <- ep_data_1_s %>%
    select(step1, e1, e2, e3) %>%
    rename("1" = e1, "1.5" = e2, "2.5" = e3) %>%
    pivot_longer(cols = c("1", "1.5", "2.5"))
  
  plots_list_ep_s[[i]] <- ggplot(ep_data_1_s, aes(x = step1, y = value, color = name)) +
    geom_line(size = 0.6) +  # Thinner lines
    geom_point(size = 1.8, alpha = 0.8) +  # Adjust point size & transparency
    scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +  # Better colors
    labs(title = "Epistemic Uncertainty",
         x = "Instances",
         y = "Value",
         color = "Standard Deviation =") +
    theme_minimal() +
    ylim(0,1)
}

combined_plot_ep_s <- (plots_list_ep_s[[1]] | plots_list_ep_s[[2]] | plots_list_ep_s[[3]]) + 
  plot_layout(guides = "collect")

