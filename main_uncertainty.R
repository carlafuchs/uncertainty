library(ggplot2)
library(dplyr)
library(tidyr)
library(viridisLite)
# Set seed for reproducibility
set.seed(123)

# Define the Bernoulli probability
p_0 <- 0.1  # Probability of success

# Number of draws
n <- 500

# Vector to store draws and results
draws <- 0
u_a_values <- numeric(n)  # Store u_a values
u_e_values <- numeric(n)  # Store u_e values

# Function to evaluate previous draws
evaluate <- function(draws, t) {

  mle <- mean(draws)  # Compute MLE from previous draws
  
  # cat("MLE", mle, "draws:", draws)
  
  # Function to solve LHS = RHS for positive phi
  equation <- function(p) {
    lhs <- dbinom(sum(draws), size = length(draws), prob = p) /
      dbinom(sum(draws), size = length(draws), prob = mle)
    rhs <- max(2 * p - 1, 0)
    return(lhs - rhs)
  }

  curve(equation(x), from = 0, to = 1, n = 100, ylim = c(-1, 1))

  # Function to solve LHS = RHS for negative phi
  equation_neg <- function(p) {
    lhs <- dbinom(sum(draws), size = length(draws), prob = p) /
      dbinom(sum(draws), size = length(draws), prob = mle)
    rhs <- max(1 - 2 * p, 0)
    return(lhs - rhs)
  }

  # Try to compute phi_pos safely
  root_pos <- tryCatch(
    uniroot(equation, c(mle, 1))$root,
    error = function(e) NA
  )

  phi_pos <- max(2 * root_pos - 1, 0)

  # Try to compute phi_neg safely
  root_neg <- tryCatch(
    uniroot(equation_neg, c(0, mle))$root,
    error = function(e) NA
  )

  phi_neg <- max(1 - 2 * root_neg, 0)
  
  # If either phi_pos or phi_neg is NA, return default values
  if (is.na(phi_pos) || is.na(phi_neg)) {
    return(c(NA, NA, NA, NA))
  }

  # # Initialize variables
  # root_pos <- mle   # Starting value
  # root_neg <- mle
  # step_size <- 0.0001  # Increment amount
  # tolerance <- 1e-3  # Convergence tolerance
  # max_iter <- 1000  # Safety limit to prevent infinite loops
  # iter <- 0
  #
  # # Define functions
  # f <- function(p) {dbinom(sum(draws), size = length(draws), prob = p) /
  #                   dbinom(sum(draws), size = length(draws), prob = mle)}  # Function we want to match
  # g <- function(p) {max(2 * p - 1, 0)} # Target value
  # h <- function(p) {max(1 - 2 * p, 0)}
  #
  #   # If either phi_pos or phi_neg is NA, return default values
  # if (is.na(f(root_pos)) || is.na(f(root_neg)) || is.na(f(mle))) {
  #   return(c(NA, NA, NA, NA))
  # }
  #
  #
  # # Increase p until f(p) is approximately equal to g(p)
  # while (abs(f(root_pos) - g(root_pos)) > tolerance && iter < max_iter) {
  #   root_pos <- root_pos + step_size  # Increase p
  #   iter <- iter + 1    # Count iterations
  #   # cat("F:", f(root_pos), "G:", g(root_pos), "Root:", root_pos)
  # }
  #
  # while (abs(f(root_neg) - h(root_neg)) > tolerance && iter < max_iter) {
  #   root_neg <- root_neg - step_size  # Increase p
  #   iter <- iter + 1    # Count iterations
  # }
  #
  # phi_pos <- max(2 * root_pos - 1, 0)
  # phi_neg <- max(1 - 2 * root_neg, 0)
  #
  #
  #
  
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
  
  return(c(u_a, u_e, p_pos, p_neg, root_pos, root_neg))
}

# Generate draws sequentially and store results
for (t in 1:n) {
  draws[t] <- rbinom(1, 1, p_0)  # Bernoulli draw
  eval_result <- evaluate(draws, t)
  
  u_a_values[t] <- eval_result[1]  # Store u_a
  u_e_values[t] <- eval_result[2]  # Store u_e
  
  cat("Step", t, "- Draw:", draws[t], "- u_a:", u_a_values[t], "- u_e:", u_e_values[t],
      "p_pos", eval_result[3], "p_neg", eval_result[4], 
      "root_pos", eval_result[5], "root_neg", eval_result[6], "\n")
}
# Create a data frame for plotting
df_01 <- data.frame(
  Step = 1:n,
  u_a = u_a_values,
  u_e = u_e_values
)


# # Create a data frame for plotting
# df_03 <- data.frame(
#   Step = 1:n,
#   u_a = u_a_values,
#   u_e = u_e_values
# )

# Create a data frame for plotting
df_05 <- data.frame(
  Step = 1:n,
  u_a = u_a_values,
  u_e = u_e_values
)

df_final <- cbind(df_01, df_03, df_05)

colnames(df_final) <- c("step", "a1", "e1", "step2", "a2", "e2","step3", "a3", "e3")

df_final_a <- df_final %>% 
  select(step, a1, a2, a3) %>%
  rename("0.1" = a1, "0.3" = a2, "0.5" = a3) %>%
  pivot_longer(cols = c("0.1", "0.3", "0.5"))


ggplot(df_final_a, aes(x = step, y = value, color = name)) +
  geom_line() +
  geom_point() +
  labs(title = "Bernoulli - Aleatoric Uncertainty",
       x = "Instances",
       y = "Value",
       color = "p =") +
  theme_minimal()

ggplot(df_final_a, aes(x = step, y = value, color = name)) +
  geom_line(size = 0.6) +  # Thinner lines
  geom_point(size = 1.8, alpha = 0.8) +  # Adjust point size & transparency
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +  # Better colors
  labs(title = "Bernoulli - Aleatoric Uncertainty",
       x = "Instances",
       y = "Value",
       color = "p =") +
  theme_minimal()

df_final_e <- df_final %>% 
  select(step, e1, e2, e3) %>%
  rename("0.1" = e1, "0.3" = e2, "0.5" = e3) %>%
  pivot_longer(cols = c("0.1", "0.3", "0.5"))


ggplot(df_final_e, aes(x = step, y = value, color = name)) +
  geom_line(size = 0.6) +  # Thinner lines
  geom_point(size = 1.8, alpha = 0.8) +  # Adjust point size & transparency
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +  # Better colors
  labs(title = "Bernoulli - Epistemic Uncertainty",
       x = "Instances",
       y = "Value",
       color = "p =") +
  theme_minimal()

# Plot u_a and u_e over time
ggplot(df, aes(x = Step)) +
  geom_line(aes(y = u_a, color = "u_a"), size = 1) +
  geom_line(aes(y = u_e, color = "u_e"), size = 1) +
  geom_point(aes(y = u_a, color = "u_a"), size = 2) +
  geom_point(aes(y = u_e, color = "u_e"), size = 2) +
  labs(title = "Evolution of u_a and u_e over Time",
       x = "Step",
       y = "Value",
       color = "Legend") +
  theme_minimal() +
  scale_color_manual(values = c("u_a" = "blue", "u_e" = "red"))



