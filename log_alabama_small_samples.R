library(alabama)

evaluate_ga_new <- function(df, x_values) {
  
  logit <- function(x) {
    log(x / (1 - x))
  }
  
  sigmoid <- function(x) {
    1 / (1 + exp(-x))
  }
  
  # Prediction function (now works for any number of covariates)
  prediction <- function(coeff, X) {
    sigmoid(X %*% coeff)  # Matrix multiplication
  }
  
  # Log-likelihood function (also generalized)
  log_likelihood <- function(coeff, X, y) {
    P <- prediction(coeff, X)
    likelihood <- sum(y * log(P) + (1 - y) * log(1 - P))
    return(-likelihood)  # Minimize negative log-likelihood
  }
  
  X <- df[, !names(df) %in% "y"]
  X <- as.matrix(cbind(1, X))
  y <- df$y
  
  mle_result <- optim(rep(c(-49.2355550, 44.9594934, 4.5014763, -3.9887503, 5.8290367, 6.1459694, -1.1921229, 3.0849678, 1.5226629,  
                            0.4943313, -4.4001744, -0.7766716, -0.5352403, 2.0734418, 5.2412949, -3.6102580, -1.7302657, 3.5612977,  
                            -5.1410566, -1.1246404, 2.7921003, -3.6933781, 2.2664529, 0.6995009, 3.5709001, -2.3721866, -7.8948901,  
                            -0.3304857, 2.8131755, 11.1628471, 2.5352575, -0.7010891, 0.4996494, -6.7268639)
  ), log_likelihood, X = X, y = y, method = "BFGS")
  mle_coef <- mle_result$par
  
  # Prepare the model matrix

  mle_lik <- log_likelihood(mle_coef, X, y)
  
  rhs_pos <- function(coeff) {max(2 * prediction(coeff, c(1, x_values)) - 1, 0)}
  rhs_neg <- function(coeff) {max(1 - 2 * prediction(coeff, c(1, x_values)), 0)}
  
  Q_p <- seq(0.5, 1, length.out = 52)[-c(1, 52)]
  Q_n <- seq(0, 0.5, length.out = 52)[-c(1, 52)]
  
  cat("Prediction", prediction(mle_coef, c(1, x_values)))
  
  phi_pos <- rhs_pos(mle_coef)
  
  cat("Phi_pos", phi_pos)
  phi_neg <- rhs_neg(mle_coef)
  
  objective_function <- function(coeff) {
    return(log_likelihood(coeff, X, y))
  }
  
  grad_objective_function <- function(coeff) {
    P <- prediction(coeff, X)
    residuals <- y - P
    grad_beta <- t(X) %*% residuals  # Compute gradient dynamically
    return(grad_beta) 
  }
  
  grad_constraint_gradient <- function(coeff) {
    return(matrix(c(1, x_values), nrow = 1))
  }
  
  for (i in 1:50) {
    
    
    alpha_p <- max(Q_p)
    alpha_n <- min(Q_n)
    
    
    cat("phi_pos", phi_pos, "alpha", alpha_p)
    
    cat("equ", 2 * alpha_p - 1 > phi_pos)
    
    if (2 * alpha_p - 1 > phi_pos) {
      
      cat("check")
      
      constraint_function <- function(coeff) {
        return(sum(coeff * c(1, x_values)) - logit(alpha_p))  # Dynamic constraint
      }
      
      # Run optimization with auglag
      result <- auglag(
        par = mle_coef, 
        fn = objective_function,  
        gr = grad_objective_function,
        heq = constraint_function,  
        heq.jac = grad_constraint_gradient
      )
      
      phi_pos <- max(phi_pos, min(exp(-result$value - mle_lik), 2 * alpha_p - 1))
      
    }
    
    if (1 - 2 * alpha_n > phi_neg) {
      
      constraint_function <- function(coeff) {
        return(sum(coeff * c(1, x_values)) - logit(alpha_n))  # Dynamic constraint
      }
      
      # Run optimization with auglag
      result <- auglag(
        par = mle_coef, 
        fn = objective_function,  
        gr = grad_objective_function,
        heq = constraint_function,  
        heq.jac = grad_constraint_gradient
      )
      
      phi_neg <- max(phi_neg, min(exp(-result$value - mle_lik), 1 - 2 * alpha_n))
      
    }
    
    Q_p <- Q_p[Q_p != alpha_p]
    Q_n <- Q_n[Q_n != alpha_n]
  }
  
  if (is.na(phi_pos) || is.na(phi_neg)) {
    return(c(NA, NA, NA, NA))
  }
  
  u_a <- min(1 - phi_pos, 1 - phi_neg)
  u_e <- min(phi_pos, phi_neg)
  
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


ionosphere <- read.table("ionosphere.data", header = FALSE, sep = ",")
colnames(ionosphere)[colnames(ionosphere) == "V35"] <- "y"
ionosphere$y <- as.numeric(sub("g", 1, sub("b", 0, ionosphere$y)))
ionosphere$V2 <- NULL

# Run Evaluations
evaluate_ga_new(ionosphere[sample(351, 50),], c(1, 0.99539, -0.05889, 0.85243, 0.02306, 0.83398, -0.37708, 1.00000, 0.03760, 
                                               0.85243, -0.17755, 0.59755, -0.44945, 0.60536, -0.38223, 0.84356, -0.38542, 0.58212,
                                               -0.32192, 0.56971, -0.29674, 0.36946, -0.47357, 0.56811, -0.51171, 0.41078, -0.46168,
                                               0.21266, -0.34090, 0.42267, -0.54487, 0.18641, -0.45300))

