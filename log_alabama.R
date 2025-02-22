library(alabama)
library(ggplot2)
library(patchwork)

evaluate_a <- function(df, x1, x2) {
  
  logit_model <- glm(y ~ X1 + X2, data = df, family = binomial)
  mle_coef <- coef(logit_model)
  
  logit <- function(x) log(x / (1 - x))
  sigmoid <- function(x) 1 / (1 + exp(-x))
  
  prediction <- function(coeff, X1, X2) {
    sigmoid(coeff[1] + coeff[2] * X1 + coeff[3] * X2)
  }
  
  log_likelihood <- function(coeff, X1, X2, y) {
    P <- prediction(coeff, X1, X2)
    sum(y * log(P) + (1 - y) * log(1 - P))
  }
  
  grad_objective_function <- function(coeff) {
    P <- prediction(coeff, df$X1, df$X2)
    residuals <- df$y - P
    c(sum(residuals), sum(residuals * df$X1), sum(residuals * df$X2))
  }
  
  constraint_function <- function(coeff, target) {
    coeff[1] + coeff[2] * x1 + coeff[3] * x2 - logit(target)
  }
  
  grad_constraint_function <- function(coeff) {
    matrix(c(1, x1, x2), nrow = 1)  # Gradient vector for equality constraint
  }
  
  objective_function <- function(coeff) {
    -log_likelihood(coeff, df$X1, df$X2, df$y)  # Negate for minimization
  }
  
  rhs_pos <- function(coeff) max(2 * prediction(coeff, x1, x2) - 1, 0)
  rhs_neg <- function(coeff) max(1 - 2 * prediction(coeff, x1, x2), 0)
  
  Q_p <- seq(0.5, 1, length.out = 52)[-c(1, 52)]
  Q_n <- seq(0, 0.5, length.out = 52)[-c(1, 52)]
  
  phi_pos <- rhs_pos(mle_coef)
  phi_neg <- rhs_neg(mle_coef)
  
  alpha_pos_values <- list(mle_coef)    
  alpha_neg_values <- list(mle_coef)
  
  for (i in 1:50) {
    
    alpha_p <- max(Q_p)
    alpha_n <- min(Q_n)
    
    if (2 * alpha_p - 1 > phi_pos) {
      
      con <- function(coeff) constraint_function(coeff, alpha_p)
      grad_con <- function(coeff) grad_constraint_function(coeff)
      
      result <- auglag(
        par = mle_coef,
        fn = objective_function,
        gr = grad_objective_function,
        heq = con,
        heq.jac = grad_con
      )
      
      phi_pos <- max(phi_pos, min(exp(-result$value - log_likelihood(mle_coef, df$X1, df$X2, df$y)), 2 * alpha_p - 1))
      alpha_pos_values[[i + 1]] <- result$par
    }
    
    if (1 - 2 * alpha_n > phi_neg) {
      
      con <- function(coeff) constraint_function(coeff, alpha_n)
      grad_con <- function(coeff) grad_constraint_function(coeff)
      
      result <- auglag(
        par = mle_coef,
        fn = objective_function,
        gr = grad_objective_function,
        heq = con,
        heq.jac = grad_con
      )
      
      phi_neg <- max(phi_neg, min(exp(-result$value - log_likelihood(mle_coef, df$X1, df$X2, df$y)), 1 - 2 * alpha_n))
      alpha_neg_values[[i + 1]] <- result$par
    }
    
    Q_p <- Q_p[Q_p != alpha_p]
    Q_n <- Q_n[Q_n != alpha_n]
  }
  
  if (is.na(phi_pos) || is.na(phi_neg)) return(c(NA, NA, NA, NA))
  
  u_a <- min(1 - phi_pos, 1 - phi_neg)
  u_e <- min(phi_pos, phi_neg)
  

  # Compute probabilities p_pos and p_neg
  if (phi_pos > phi_neg) {
    p_pos <- 1 - (u_a + u_e)
  } else if (phi_pos == phi_neg) {
    p_pos <- (1 - (u_a + u_e)) / 2
  } else {
    p_pos <- 0
  }
  
  p_neg <- 1 - (p_pos + u_a + u_e)
  
  return(c(u_a, u_e, p_pos, p_neg))
}


# Function to generate synthetic data
generate_data <- function(n = 200, separation = 2, noise = 0.1, ood_fraction = 0.1) {
  
  # Generate two class centers
  mu_1 <- c(-separation, -separation)  # Class 0 center
  mu_2 <- c(separation, separation)    # Class 1 center
  
  # Generate class 0 samples
  X0 <- matrix(rnorm(n, mean = rep(mu_1, each = n/2), sd = 1), ncol = 2)
  y0 <- rep(0, n/2)  # Class 0 labels
  
  # Generate class 1 samples
  X1 <- matrix(rnorm(n, mean = rep(mu_2, each = n/2), sd = 1), ncol = 2)
  y1 <- rep(1, n/2)  # Class 1 labels
  
  # Combine datasets
  X <- rbind(X0, X1)
  y <- c(y0, y1)
  
  # Add aleatoric uncertainty (flip labels with noise probability)
  flip_indices <- sample(1:n, size = round(noise * n), replace = FALSE)
  y[flip_indices] <- 1 - y[flip_indices]  # Flip labels
  
  # Generate Out-of-Distribution (OOD) samples
  num_ood <- round(n * ood_fraction)
  X_ood <- matrix(rnorm(num_ood * 2, mean = 0, sd = 4), ncol = 2)  # More spread out
  y_ood <- rep(NA, num_ood)  # Mark OOD samples with NA (unknown class)
  
  # Merge OOD samples with training data
  X <- rbind(X, X_ood)
  y <- c(y, y_ood)
  
  # Create dataframe
  df <- data.frame(X1 = X[,1], X2 = X[,2], y = y)
  
  return(df)
}



highlight_x <- c(-2, 0, 0)
highlight_y <- c(-2, -2, 0)
plots_list <- list()

for (i in 1:3) {
  plots_list[[i]] <- ggplot(df, aes(x = X1, y = X2, color = factor(y))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_point(x = highlight_x[i], y = highlight_y[i], color = "red", size = 4) +  # Highlighted point
  geom_text(x = highlight_x[i], y = highlight_y[i], label = expression(x[q]), 
            vjust = -1, hjust = 1, color = "red", fontface = "bold") +  # Annotation
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "black"), na.value = "gray") +
  labs(title = "Synthetic Data for Logistic Regression",
       color = "Class")  +
  ggtitle(NULL) +
  theme_minimal()
}

# Plot the data


combined_plot <- (plots_list[[1]] | plots_list[[2]] | plots_list[[3]]) + 
  plot_layout(guides = "collect")


results_1_03 <- data.frame(index = seq.int(from = 20, to = 500, by = 10), u_a = NA, u_e = NA, p_pos = NA, p_neg = NA)

for (i in seq.int(from = 20, to = 500, by = 10)) {
  df <- generate_data(n = i, separation = 2, noise = 0.3, ood_fraction = 0)
  results_1_03[results_1_03$index == i, c(2:5)] <- evaluate_a(df, 0, 0)
}

x1 <- c(-2, 0, 0)
x2 <- c(-2, -2, 0)
noise_v <- c(0.1, 0.3, 0.5)
separation_v <- c(0.4, 1.2, 2)

list_results <- list()

for (i in 1:3) {
  
  list_results[[i]] <- list()
  
  for(j in 1:3) {
    list_results[[i]][[j]] <- data.frame(index = seq.int(from = 20, to = 500, by = 10), u_a = NA, u_e = NA, p_pos = NA, p_neg = NA)
    for (k in seq.int(from = 20, to = 500, by = 10)) {
      df <- generate_data(n = k, separation = 2, noise = noise_v[j], ood_fraction = 0)
      list_results[[i]][[j]][list_results[[i]][[j]]$index == k, c(2:5)] <- evaluate_a(df, x1[i], x2[i])
    }
  }
  
  for(j in 1:3) {
    list_results[[i]][[j + 3]] <- data.frame(index = seq.int(from = 20, to = 500, by = 10), u_a = NA, u_e = NA, p_pos = NA, p_neg = NA)
    for (k in seq.int(from = 20, to = 500, by = 10)) {
      df <- generate_data(n = k, separation = separation_v[j], noise = 0.1, ood_fraction = 0)
      list_results[[i]][[j + 3]][list_results[[i]][[j + 3]]$index == k, c(2:5)] <- evaluate_a(df, x1[i], x2[i])
    }
  }
  
}

evaluate_a(df, 0, 0)


df <- generate_data(n = 100, separation = 0.4, noise = 0.1, ood_fraction = 0)
ggplot(df, aes(x = X1, y = X2, color = factor(y))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "black"), na.value = "gray") +
  labs(title = "Synthetic Data for Logistic Regression",
       color = "Class")  +
  ggtitle(NULL) +
  theme_minimal()
