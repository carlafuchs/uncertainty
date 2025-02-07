
evaluate <- function(df, x_values) {
  
  # Fit logistic regression model (MLE estimation)
  logit_model <- glm(y ~ ., data = df, family = binomial)  # Use all covariates dynamically
  
  cat(logit_model$fitted.values)
  # Extract MLE estimates
  mle_coef <- coef(logit_model)
  
  cat("MLE Coefficients: ", mle_coef, "\n")
  
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
    return(likelihood)
  }
  
  # Prepare the model matrix
  X <- model.matrix(logit_model)  # Extract design matrix
  y <- df$y
  
  mle_lik <- log_likelihood(mle_coef, X, y)
  
  rhs_pos <- function(coeff) {max(2 * prediction(coeff, c(1, x_values)) - 1, 0)}
  rhs_neg <- function(coeff) {max(1 - 2 * prediction(coeff, c(1, x_values)), 0)}
  
  Q_p <- seq(0.5, 1, length.out = 52)[-c(1, 52)]
  Q_n <- seq(0, 0.5, length.out = 52)[-c(1, 52)]
  
  phi_pos <- rhs_pos(mle_coef)
  phi_neg <- rhs_neg(mle_coef)
  
  alpha_pos_values <- list()      
  alpha_neg_values <- list()
  alpha_pos_values[[1]] <- mle_coef
  alpha_neg_values[[1]] <- mle_coef
  
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
    return(c(1, x_values))  # Gradient vector for equality constraint
  }
  
  for (i in 1:50) {
    
    alpha_p <- max(Q_p)
    alpha_n <- min(Q_n)
    
    if(2 * alpha_p - 1 > phi_pos) {
      
      constraint_function <- function(coeff) {
        return(sum(coeff * c(1, x_values)) - logit(alpha_p))  # Dynamic constraint
      }
      
      # Run optimization with constraints
      result <- nloptr(
        x0 = alpha_pos_values[[i]], 
        eval_f = objective_function,  
        eval_grad_f = grad_objective_function,
        eval_g_ineq = NULL,           
        eval_g_eq = constraint_function,  
        eval_jac_g_eq = grad_constraint_gradient,
        opts = list(
          maximize = TRUE,
          "algorithm" = "NLOPT_LD_SLSQP",
          "maxeval" = 500,
          "xtol_rel" = 1e-6
        )
      )
      
      cat("POS: ", result$solution, "\n")
      
      phi_pos <- max(phi_pos, min(exp(result$objective - mle_lik), 2 * alpha_p - 1))
      
      cat("PHI_POS: ", phi_pos, "\n")
      
      alpha_pos_values[[i + 1]] <- result$solution
    }
    
    if(1 - 2 * alpha_n > phi_neg) {
      
      constraint_function <- function(coeff) {
        return(sum(coeff * c(1, x_values)) - logit(alpha_n))  # Dynamic constraint
      }
      
      # Run optimization with constraints
      result <- nloptr(
        x0 = alpha_neg_values[[i]], 
        eval_f = objective_function,  
        eval_grad_f = grad_objective_function,
        eval_g_ineq = NULL,           
        eval_g_eq = constraint_function,  
        eval_jac_g_eq = grad_constraint_gradient,
        opts = list(
          maximize = TRUE,
          "algorithm" = "NLOPT_LD_SLSQP",
          "maxeval" = 500,
          "xtol_rel" = 1e-6
        )
      )
      
      cat("NEG: ", result$solution, "\n")
      
      phi_neg <- max(phi_neg, min(exp(result$objective - mle_lik), 1 - 2 * alpha_n))
      
      cat("PHI_NEG: ", phi_neg, "\n")
      
      alpha_neg_values[[i + 1]] <- result$solution
    }
    
    Q_p <- Q_p[Q_p != alpha_p]
    Q_n <- Q_n[Q_n != alpha_n]
  }
  
  if (is.na(phi_pos) || is.na(phi_neg)) {
    return(c(NA, NA, NA, NA))
  }
  
  cat("pp", phi_pos, "\n")
  cat("pn", phi_neg, "\n")
  
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

df <- read.table("data_banknote_authentication.txt", header = FALSE, sep = ",")
colnames(df)[colnames(df) == "V5"] <- "y"

ionosphere <- read.table("ionosphere.data", header = FALSE, sep = ",")
colnames(ionosphere)[colnames(ionosphere) == "V35"] <- "y"
ionosphere$y <- sub("g", 1, ionosphere$y)
ionosphere$y <- sub("b", 0, ionosphere$y)
ionosphere$y <- as.numeric(ionosphere$y)
ionosphere$V2 <- NULL



evaluate(ionosphere[sample(1:351, 50),], c(1, 0.99539, -0.05889, 0.85243 ,0.02306 ,0.83398 ,-0.37708 ,1.00000 ,0.03760 , 
                       0.85243 ,-0.17755 ,0.59755 ,-0.44945 ,0.60536 ,-0.38223 ,0.84356 ,-0.38542 ,0.58212 ,
                       -0.32192 ,0.56971 ,-0.29674 ,0.36946 ,-0.47357 ,0.56811 ,-0.51171 ,0.41078 ,-0.46168 ,
                       0.21266 ,-0.34090 ,0.42267 ,-0.54487 ,0.18641 ,-0.45300))


evaluate(df, c(0,0,0,0))
evaluate(df[1:10,], c(4.5459000, 8.16740, -2.45860000, -1.462100))

