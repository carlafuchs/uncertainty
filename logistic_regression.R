# Load necessary library
library(ggplot2)
library(cmaes)
library(nloptr)
library(alabama)

# Set seed for reproducibility
set.seed(123)

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


### MAIN PART

# # Function to evaluate previous draws
# evaluate <- function(df, x1, x2) {
#   
#   # Fit logistic regression model (MLE estimation)
#   logit_model <- glm(y ~ X1 + X2, data = df, family = binomial)
#   
#   # Extract MLE estimates
#   mle_coef <- coef(logit_model)
#   
#   print(df$y)
#   print(prod(df$y))
#   
#   prediction <- function(coeff) {
#     1 / (1 + exp(-(coeff[1] + coeff[2] * x1 + coeff[3] * x2)))
#   }
#   
#   likelihood <- function(coeff, y) {
#     P <- prediction(coeff)  # Compute probabilities
#     likelihood <- prod(P^y * (1 - P)^(1 - y))   # Compute likelihood
#     return(likelihood)
#   }
#   
#   # cat("MLE", mle_coefficients)
#   
#   # Define lhs and rhs
#   lhs <- function(coeff) {likelihood(coeff, df$y) / likelihood(mle_coef, df$y)}
#   
#   rhs_pos <- function(coeff) {max(2 * prediction(coeff) - 1, 0)}
#   
#   rhs_neg <- function(coeff) {max(1 - 2 * prediction(coeff), 0)}
#   
#   print(mle_coef)
#   
#   print(lhs(mle_coef))
#   print(rhs_pos(mle_coef))
#   print(rhs_neg(mle_coef))
#   
#   # Try to compute root_pos safely
#   # Define the objective function (negate for minimization)
#   objective_function_pos <- function(x) {
#     f_value <- -min(lhs(x), rhs_pos(x))  # We negate since CMA-ES minimizes
#     
#     # Define acceptable function range
#     if (f_value < 0 || f_value > 1) {
#       return(Inf)  # Reject solution by assigning a very high cost
#     }
#     
#     return(f_value)
#   }
#   
#   # Set dimension (number of variables)
#   dimension <- 3  # Change based on your problem
#   
#   # Run CMA-ES optimization
#   result_pos <- cma_es(
#     par = mle_coef,  # Initial point (e.g., zero vector)
#     fn = objective_function_pos,  # Function to optimize
#     lower = rep(-500, dimension),  # Lower bounds
#     upper = rep(500, dimension),   # Upper bounds
#     control = list(trace = TRUE, maxit = 500, sigma = 0.01)  # Maximum iterations
#   )
#   
#   root_pos <- result_pos$par
#   
#   phi_pos <- result_pos$value
#   
#   # Try to compute root_neg safely
#   
#   # Define the objective function (negate for minimization)
#   objective_function_neg <- function(x) {
#     f_value <- -min(lhs(x), rhs_neg(x))  # We negate since CMA-ES minimizes
#     
#     # Define acceptable function range
#     if (f_value < 0 || f_value > 1) {
#       return(Inf)  # Reject solution by assigning a very high cost
#     }
#     
#     return(f_value)  # Otherwise, return the normal function value
#   }
#   
#   # Set dimension (number of variables)
#   dimension <- 3  # Change based on your problem
#   
#   # Run CMA-ES optimization
#   result_neg <- cma_es(
#     par = mle_coef,  # Initial point (e.g., zero vector)
#     fn = objective_function_neg,  # Function to optimize
#     lower = rep(-500, dimension),  # Lower bounds
#     upper = rep(500, dimension),   # Upper bounds
#     control = list(maxit = 500)  # Maximum iterations
#   )
# 
#   root_neg <- result_neg$par
#   
#   phi_neg <- result_neg$value
#   
#   # If either phi_pos or phi_neg is NA, return default values
#   if (is.na(phi_pos) || is.na(phi_neg)) {
#     return(c(NA, NA, NA, NA))
#   }
#   
#   # Compute u_a and u_e
#   u_a <- min(1 - phi_pos, 1 - phi_neg)
#   u_e <- min(phi_pos, phi_neg)
#   
#   # Compute probabilities p_pos and p_neg
#   if (phi_pos > phi_neg) {
#     p_pos <- 1 - (u_a + u_e)
#   } else if (phi_pos == phi_neg) {
#     p_pos <- (1 - (u_a + u_e)) / 2
#   } else {
#     p_pos <- 0
#   }
#   
#   p_neg <- 1 - (p_pos + u_a + u_e)
#   
#   return(c(u_a, u_e, p_pos, p_neg))
# }


### MAIN PART

# Function to evaluate previous draws
evaluate <- function(df, x1, x2) {
  
  # Fit logistic regression model (MLE estimation)
  logit_model <- glm(y ~ X1 + X2, data = df, family = binomial)
  
  # Extract MLE estimates
  mle_coef <- coef(logit_model)
  
  cat("MLE: ", mle_coef)
  
  logit <- function(x) {
    log(x / (1 - x))
  }
  
  sigmoid <- function(x) {
    1 / (1 + exp(-x))
  }
  
  prediction <- function(coeff, X1, X2) {
   sigmoid(coeff[1] + coeff[2] * X1 + coeff[3] * X2)
  }
  
  log_likelihood <- function(coeff, X1, X2, y) {
    P <- prediction(coeff, X1, X2)  # Compute probabilities
    likelihood <- sum(y * log(P) + (1 - y) * log(1 - P))   # Compute likelihood
    return(likelihood)
  }
  
  mle_lik <- log_likelihood(mle_coef, df$X1, df$X2, df$y)
  
  rhs_pos <- function(coeff) {max(2 * prediction(coeff, x1, x2) - 1, 0)}
  rhs_neg <- function(coeff) {max(1 - 2 * prediction(coeff, x1, x2), 0)}
  
  Q_p <- seq(0.5, 1, length.out = 52)[-c(1, 52)]
  Q_n <- seq(0, 0.5, length.out = 52)[-c(1, 52)]
  
  phi_pos <- rhs_pos(mle_coef)
  phi_neg <- rhs_neg(mle_coef)
  
  alpha_pos_values <- list()      
  alpha_neg_values <- list()
  alpha_pos_values[[1]] <- mle_coef
  alpha_neg_values[[1]] <- mle_coef
  
  objective_function <- function(coeff) {
    return(log_likelihood(coeff, df$X1, df$X2, df$y))
  }
  
  grad_objective_function <- function(coeff) {
    # Compute predicted probabilities
    P <- prediction(coeff, df$X1, df$X2)
    
    # Compute residuals (difference between actual and predicted values)
    residuals <- df$y - P
    
    # Compute gradients
    grad_beta0 <- sum(residuals)  # Intercept term
    grad_beta1 <- sum(residuals * df$X1)
    grad_beta2 <- sum(residuals * df$X2)
    
    return(c(grad_beta0, grad_beta1, grad_beta2))
  }
  
  grad_constraint_gradient <- function(coeff) {
    return(c(1, x1, x2))  # Gradient vector
  }
  
  for (i in 1:50) {
    
    alpha_p <- max(Q_p)
    alpha_n <- min(Q_n)
    
    print(alpha_p)
    print(alpha_n)
    
    if(2 * alpha_p - 1 > phi_pos) {

      constraint_function <- function(coeff) {
        return(coeff[1] + coeff[2] * x1 + coeff[3] * x2 - logit(alpha_p))
      }
      
      cat("constraint", constraint_function(mle_coef))
      
      # Run optimization with constraints
      result <- nloptr(
        x0 = c(0,0,0), 
        eval_f = objective_function,  # Objective function
        eval_grad_f = grad_objective_function,
        eval_g_ineq = NULL,           # No inequality constraints
        eval_g_eq = constraint_function,  # Equality constraint
        eval_jac_g_eq = grad_constraint_gradient,
        opts = list(
          maximize = TRUE,
          "algorithm" = "NLOPT_LD_SLSQP",  # Sequential quadratic programming
          "maxeval" = 500,
          "xtol_rel" = 1e-20
        )
      )
      
      cat("POS: ", result$solution)
      
      phi_pos <- max(phi_pos, min(exp(result$objective - mle_lik), 2 * alpha_p - 1))
      
      cat("PHI_POS: ", phi_pos)
      
      alpha_pos_values[[i + 1]] <- result$solution
    }
    
    if(1 - 2 * alpha_n > phi_neg) {
      
      constraint_function <- function(coeff) {
        return(coeff[1] + coeff[2] * x1 + coeff[3] * x2 - logit(alpha_n))
      }
      
      # Run optimization with constraints
      result <- nloptr(
        x0 = c(0,0,0), 
        eval_f = objective_function,  # Objective function
        eval_grad_f = grad_objective_function,
        eval_g_ineq = NULL,           # No inequality constraints
        eval_g_eq = constraint_function,  # Equality constraint
        eval_jac_g_eq = grad_constraint_gradient,
        opts = list(
          maximize = TRUE,
          "algorithm" = "NLOPT_LD_SLSQP",  # Sequential quadratic programming
          "maxeval" = 500,
          "xtol_rel" = 1e-20
        )
      )
      
      cat("NEG: ", result$solution)
      
      phi_neg <- max(phi_neg, min(exp(result$objective - mle_lik), 1 - 2 * alpha_n))
      
      cat("PHI_NEG: ", phi_neg)
      
      alpha_neg_values[[i + 1]] <- result$solution
    }
    
    Q_p <- Q_p[Q_p != alpha_p]
    Q_n <- Q_n[Q_n != alpha_n]
  }
  
  # If either phi_pos or phi_neg is NA, return default values
  if (is.na(phi_pos) || is.na(phi_neg)) {
    return(c(NA, NA, NA, NA))
  }
  
  cat("pp", phi_pos)
  cat("pn", phi_neg)
  
  # Compute u_a and u_e
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


# try out for observation

evaluate(df, 0, 0)

# Generate dataset
df <- generate_data(n = 20, separation = 2, noise = 0.1, ood_fraction = 0)

# Plot the data
ggplot(df, aes(x = X1, y = X2, color = factor(y))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("blue", "red", "black"), na.value = "gray") +
  labs(title = "Synthetic Data for Logistic Regression",
       color = "Class (NA = OOD)") +
  theme_minimal()


al_u <- c()
ep_u <- c()
result <- list()
  
for (i in 1:5) {
  df <- generate_data(n = 100*i, separation = 2, noise = 0.3, ood_fraction = 0)
  result[[i]] <- evaluate(df, -2.5, -2.5)
  al_u[i] <- result[[i]][1]
  ep_u[i] <- result[[i]][2]
}

plot(al_u)
plot(ep_u)
