library(alabama)

evaluate <- function(df, x1, x2) {
  
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
  
  p_pos <- ifelse(phi_pos > phi_neg, 1 - (u_a + u_e), (1 - (u_a + u_e)) / 2)
  p_neg <- 1 - (p_pos + u_a + u_e)
  
  return(c(u_a, u_e, p_pos, p_neg))
}

df <- generate_data(n = 150, separation = 2, noise = 0.1, ood_fraction = 0)

evaluate(df, 0, 0)
